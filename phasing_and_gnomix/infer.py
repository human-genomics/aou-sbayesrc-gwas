"""Bounded sample batches; Gnofix is mandatory for every input individual."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import gzip
from pathlib import Path
import random
import time
import numpy as np
import pandas as pd
from phasing_and_gnomix.common import configure_runtime, derived_seed, read_psam, complete, write_json
from phasing_and_gnomix.models import load_model

_STATE = {}


def reconstruct_haplotypes(X, swap, C, M, W):
    window = np.minimum(np.arange(C) // M, W-1)
    pairs = X.reshape(-1,2,C)
    return np.where(swap[:,None,window], pairs[:,::-1,:], pairs).reshape(X.shape)


def model_input(pgen, indices, n, v, C, match):
    import pgenlib
    X = np.zeros((len(indices)*2,C), dtype=np.int8)
    order = np.argsort(match['pgen_idx'], kind='stable')
    variants = match['pgen_idx'][order].astype(np.uint32)
    with pgenlib.PgenReader(str(pgen).encode(), raw_sample_ct=n, variant_ct=v,
                            sample_subset=np.asarray(indices,dtype=np.uint32)) as reader:
        for s in range(0,len(order),4096):
            e=min(s+4096,len(order))
            alleles=np.empty((e-s,len(indices)*2),np.int32)
            reader.read_alleles_list(variants[s:e],alleles)
            if (alleles < 0).any() or (alleles > 1).any():
                raise ValueError('Unexpected missing or non-binary Beagle allele')
            a=alleles.T.astype(np.int8)
            f=match['flip'][order[s:e]]
            a[:,f]=1-a[:,f]
            X[:,match['model_idx'][order[s:e]]]=a
    return X


def calls(model, X, iids, chrom, seed):
    from src.Gnofix.gnofix import gnofix
    B=model.base.predict_proba(X)
    raw=model.smooth.predict_proba(B).argmax(axis=-1)
    W=model.W
    fixed=np.empty((len(X),W),np.uint8)
    swaps=np.zeros((len(iids),W),bool)
    for k,iid in enumerate(iids):
        per_sample=derived_seed(seed, chrom, iid)
        random.seed(per_sample);np.random.seed(per_sample)
        _,_,ym,yp,_,tracker=gnofix(X[2*k],X[2*k+1],B=B[2*k:2*k+2],smoother=model.smooth)
        fixed[2*k],fixed[2*k+1]=ym,yp
        swaps[k]=np.asarray(tracker[0])==1
    corrected=reconstruct_haplotypes(X,swaps,model.C,model.M,model.W)
    probabilities=model.predict_proba(corrected)
    support=np.take_along_axis(probabilities,fixed[...,None].astype(np.int64),axis=2)[...,0]
    if not np.isfinite(support).all() or (support < 0).any() or (support > 1).any():
        raise ValueError('Invalid posterior support')
    agreement=(np.sort(raw.reshape(-1,2,W),axis=1)==np.sort(fixed.reshape(-1,2,W),axis=1)).all(axis=1).mean(axis=1)
    return fixed,support,swaps,agreement


def _init(runtime, model, pgen, match, samples, windows, chrom, seed):
    configure_runtime(runtime)
    _STATE.update(model=load_model(model),pgen=pgen,match=dict(np.load(match)),
                  samples=read_psam(samples),win=pd.read_csv(windows,sep='\t'),chrom=chrom,seed=seed)


def _batch(task):
    from gnomix1000g.tracts import haplotype_tracts,window_lengths_cm
    start,end,path=task
    st=_STATE;m=st['model'];win=st['win']
    iids=st['samples'].IID.values[start:end].astype(str)
    X=model_input(st['pgen'],np.arange(start,end),len(st['samples']),int(st['match']['pgen_variants']),m.C,st['match'])
    labels,support,swaps,agreement=calls(m,X,iids,st['chrom'],st['seed'])
    tracts=haplotype_tracts(labels,support,win,iids,st['chrom'])
    counts=tracts.groupby(['sample','haplotype']).n_windows.sum()
    if len(counts)!=2*len(iids) or not (counts==m.W).all():
        raise ValueError('Tracts do not cover all windows on both haplotypes')
    weights=window_lengths_cm(win)
    totals=np.zeros((len(iids),len(m.population_order)),np.float64)
    for a in range(len(m.population_order)):
        totals[:,a]=((labels.reshape(-1,2,m.W)==a)*weights).sum(axis=(1,2))
    path=Path(path)
    tracts.to_csv(path.with_suffix('.tracts.tsv.gz'),sep='\t',index=False,compression='gzip')
    np.savez_compressed(path.with_suffix('.npz'),samples=iids,totals=totals,swaps=swaps,raw_fix_diploid_agreement=agreement,
                        called_posterior=support.astype(np.float32),labels=labels)
    return {'start':start,'end':end,'tract_rows':len(tracts),'gnofix_completed':len(iids)}


def run_infer(args,config,manifest):
    started=time.time()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    parts=out/'parts';parts.mkdir(exist_ok=True)
    samples=read_psam(args.psam)
    cohort=pd.read_csv(args.samples,sep='\t',dtype={'IID':str,'FID':str})
    left=args.shard*config['shard_size']
    expected=cohort.iloc[left:left+config['shard_size']]
    if list(samples.IID)!=list(expected.IID):
        raise ValueError('Inference shard is not the expected cohort slice')
    tasks=[(s,min(s+config['batch_size'],len(samples)),str(parts/f'part_{s:05d}')) for s in range(0,len(samples),config['batch_size'])]
    with ProcessPoolExecutor(config['infer_cores'],mp_context=get_context('spawn'),initializer=_init,
                             initargs=(args.runtime,args.model,args.pgen,args.match,args.psam,args.windows,args.chrom,config['seeds']['inference'])) as ex:
        stats=list(ex.map(_batch,tasks))
    totals=np.zeros((len(samples),len(config['ancestries'])),np.float64)
    win=pd.read_csv(args.windows,sep='\t')
    W=len(win)
    switch_rows=[]
    completion=[]
    with gzip.open(out/'tracts_gnofix.tsv.gz','wt') as merged:
        for j,(start,end,prefix) in enumerate(tasks):
            prefix=Path(prefix)
            with gzip.open(prefix.with_suffix('.tracts.tsv.gz'),'rt') as f:
                header=f.readline()
                if j==0: merged.write(header)
                for line in f: merged.write(line)
            with np.load(prefix.with_suffix('.npz')) as d:
                if list(d['samples'])!=list(samples.IID.iloc[start:end]):
                    raise ValueError('Inference batch changed sample order')
                totals[start:end]=d['totals']
                s,w=np.nonzero(np.diff(d['swaps'].astype(np.int8),axis=1,prepend=0)!=0)
                switch_rows.append(pd.DataFrame({'sample':d['samples'][s],'chrom':args.chrom,'window':w,
                                                  'pos_hg19':win.spos_hg19.values[w],'pos_hg38':win.spos_hg38.values[w],
                                                  'cM':win.sgpos.values[w]}))
                completion.extend({'sample':iid,'chrom':args.chrom,'gnofix_completed':True,'n_windows_per_haplotype':W,
                                   'n_switches':int(np.count_nonzero(s==k)), 'raw_fix_diploid_agreement':float(d['raw_fix_diploid_agreement'][k])}
                                  for k,iid in enumerate(d['samples']))
    pd.concat(switch_rows,ignore_index=True).to_csv(out/'gnofix_switches.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(completion).to_csv(out/'sample_status.tsv.gz',sep='\t',index=False,compression='gzip')
    np.savez_compressed(out/'ancestry_totals.npz',samples=samples.IID.values.astype(str),totals=totals,populations=np.asarray(config['ancestries']))
    if sum(x['gnofix_completed'] for x in stats)!=len(samples):
        raise ValueError('Incomplete Gnofix processing')
    summary={'stage':'infer','chrom':args.chrom,'shard':args.shard,'samples':len(samples),'windows':W,
             'gnofix_completed':len(samples),'tract_rows':sum(x['tract_rows'] for x in stats),
             'wall_seconds':time.time()-started,'seeds':config['seeds']}
    # Keep compact window calls/posteriors as audit checkpoints; duplicate tract parts are temporary.
    for p in parts.glob('*.tracts.tsv.gz'):p.unlink()
    complete(out,manifest['manifest_id'],summary)
