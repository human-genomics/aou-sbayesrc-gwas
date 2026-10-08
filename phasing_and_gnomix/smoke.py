"""Exercise real Gnomix, shard I/O, and output generation on public demo data."""
from __future__ import annotations
import argparse
from pathlib import Path
import tempfile
from types import SimpleNamespace
import numpy as np
import pandas as pd
from .common import configure_runtime,read_json,write_json,verify_complete,run,read_psam
from .models import load_model,window_table
from .infer import calls
from .worker import run_union_infer


def beagle_smoke(runtime, config):
    from .tests.real_tools import validate
    with tempfile.TemporaryDirectory(prefix='union-beagle-smoke-') as temp:
        validate(runtime, Path(temp))
    return {'real_beagle_phase_passed':True, 'observed_genotypes_preserved':True,
            'shard_phase_verified':True, 'real_union_selection_passed':True}


def smoke(runtime,model_path,meta_path,config,output):
    configure_runtime(runtime)
    import pgenlib
    from src.utils import read_vcf,vcf_to_npy
    m=load_model(model_path)
    vcf=read_vcf(str(Path(runtime)/'gnomix/demo/data/small_query_chr22.vcf.gz'),chm='22',fields='*')
    X=vcf_to_npy(vcf,m.snp_pos,m.snp_ref,miss_fill=0,verbose=False)
    samples=list(vcf['samples']);n=len(samples)
    expected=calls(m,X,samples,22,config['seeds']['inference'])
    with tempfile.TemporaryDirectory(prefix='lai-public-smoke-') as temp:
        p=Path(temp)
        with pgenlib.PgenWriter(str(p/'input.pgen').encode(),n,m.C+2,False,2,True) as writer:
            for start in range(0,m.C,4096):
                writer.append_alleles_batch(np.ascontiguousarray(X[:,start:start+4096].T,dtype=np.int32),all_phased=True)
            extra=np.tile(np.array([0,1],dtype=np.int32),n)
            writer.append_alleles_batch(np.array([extra,1-extra],dtype=np.int32),all_phased=True)
        pd.DataFrame({'#FID':['0']*n,'IID':samples}).to_csv(p/'input.psam',sep='\t',index=False)
        pd.DataFrame({'FID':['0']*n,'IID':samples}).to_csv(p/'samples.tsv',sep='\t',index=False)
        np.savez(p/'match.npz',model_idx=np.arange(m.C),pgen_idx=np.arange(m.C),flip=np.zeros(m.C,bool),pgen_variants=np.asarray(m.C+2))
        meta=dict(np.load(meta_path))
        win=window_table(meta,np.arange(m.C),22)
        win.to_csv(p/'windows.tsv',sep='\t',index=False)
        cfg=dict(config,batch_size=4,infer_cores=2)
        args=SimpleNamespace(runtime=str(Path(runtime).resolve()),model=str(Path(model_path).resolve()),
                             pgen=str(p/'input.pgen'),psam=str(p/'input.psam'),samples=str(p/'samples.tsv'),
                             match=str(p/'match.npz'),windows=str(p/'windows.tsv'),chrom=22,shard=0,output=str(p/'out'))
        good=np.flatnonzero(meta['pos38']>0)
        idx=[int(good[0]),int(good[-1]),m.C,m.C+1]
        windows=np.minimum(np.asarray([idx[0],idx[1],idx[1],idx[1]])//m.M,m.W-1)
        positions=[int(meta['pos38'][idx[0]]),int(meta['pos38'][idx[1]]),int(meta['pos38'][idx[1]])+1,int(meta['pos38'][idx[1]])+2]
        mapping=pd.DataFrame({'chrom':22,'pos':positions,'id':['demo'+str(i) for i in idx],
                              'ref_source':['A']*4,'alt_source':['C']*4,'source_variant_index':idx,'window':windows})
        mapping.to_csv(p/'scoring.tsv.gz',sep='\t',index=False)
        import os
        previous=os.environ.get('SCORING_MAP')
        os.environ['SCORING_MAP']=str(p/'scoring.tsv.gz')
        try:
            run_union_infer(args,cfg,{'manifest_id':'public-smoke'},{'policy_id':'public-smoke'})
        finally:
            if previous is None:os.environ.pop('SCORING_MAP',None)
            else:os.environ['SCORING_MAP']=previous
        with pgenlib.PgenReader(str(p/'out/scoring_genotypes_gnofix.pgen').encode(),raw_sample_ct=n) as reader:
            for k,row in mapping.iterrows():
                before=X[:,row.source_variant_index].astype(np.int32) if row.source_variant_index<m.C else (extra if row.source_variant_index==m.C else 1-extra)
                corrected=before.reshape(n,2).copy()
                mask=expected[2][:,row.window]
                corrected[mask]=corrected[mask,::-1]
                actual_alleles=np.empty(2*n,np.int32)
                reader.read_alleles(k,actual_alleles)
                np.testing.assert_array_equal(actual_alleles,corrected.ravel())
        assert read_json(p/'out/COMPLETE.json')['independent_tract_audit']
        assert verify_complete(p/'out','public-smoke')
        actual=np.concatenate([np.load(part)['labels'] for part in sorted((p/'out/parts').glob('*.npz'))])
        assert np.array_equal(actual,expected[0])
        tracts=pd.read_csv(p/'out/tracts_gnofix.tsv.gz',sep='\t')
        assert list(tracts.columns)==['sample','haplotype','chrom','start_hg38','end_hg38','start_hg19','end_hg19','start_cM','end_cM','ancestry','n_windows','mean_posterior']
        assert (tracts.groupby(['sample','haplotype']).n_windows.sum()==m.W).all()
        with np.load(p/'out/ancestry_totals.npz') as totals:
            assert (totals['totals']>=0).all() and (totals['totals'].sum(axis=1)>0).all()
            assert totals['samples'].dtype.kind=='U'
        status=pd.read_csv(p/'out/sample_status.tsv.gz',sep='\t')
        assert status.gnofix_completed.all() and len(status)==n
    phasing=beagle_smoke(runtime,config)
    write_json(output,{'public_end_to_end_inference_passed':True,'all_samples_gnofix':True,'tract_schema_valid':True,'batch_labels_identical':True,'extra_union_sites_preserve_model_predictions':True,'scoring_phase_verified':True,'independent_tract_audit_passed':True,**phasing})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('runtime','model','meta','config','output'):p.add_argument('--'+name,required=True)
    a=p.parse_args()
    smoke(a.runtime,a.model,a.meta,read_json(a.config),a.output)
