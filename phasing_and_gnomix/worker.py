"""Frozen union extraction, Beagle, Gnomix/Gnofix and scoring-genotype worker."""
from __future__ import annotations
import argparse,hashlib,os
from pathlib import Path
import numpy as np
import pandas as pd
from phasing_and_gnomix.common import complete,configure_runtime,read_json,read_psam,sha256,write_json


def audit_outputs(a):
    """Check tract coverage and posteriors against independent saved checkpoints."""
    from .audit import audit_unit_data
    out = Path(a.output)
    samples = read_psam(a.psam).IID.to_numpy().astype(str)
    def checkpoints():
        for path in sorted((out/'parts').glob('*.npz')):
            start = int(path.stem.split('_')[1])
            with np.load(path, allow_pickle=False) as z:
                part = {name:z[name] for name in z.files}
            yield start, start+len(part['samples']), part
    with np.load(out/'ancestry_totals.npz', allow_pickle=False) as z:
        totals = {name:z[name] for name in z.files}
    _, report = audit_unit_data(
        pd.read_csv(out/'tracts_gnofix.tsv.gz', sep='\t', dtype={'sample':str}),
        pd.read_csv(out/'sample_status.tsv.gz', sep='\t', dtype={'sample':str}),
        pd.read_csv(out/'gnofix_switches.tsv.gz', sep='\t', dtype={'sample':str}),
        totals, checkpoints(), pd.read_csv(a.windows, sep='\t'), samples, a.chrom)
    write_json(out/'tract_audit.json', report)


def run_union_phase(a,cfg,manifest,policy):
    from phasing_and_gnomix.phase import run_phase
    row=policy['chromosomes'][str(a.chrom)]
    if row['threads']!=cfg['beagle_threads']:raise ValueError('Beagle threads differ from frozen chromosome settings')
    a.heap_gb=row['heap_gb']
    run_phase(a,cfg,manifest)
    out=Path(a.output);stats=read_json(out/'phasing.json');selection=read_json(out/'extraction.json')
    stats.update(union_variants=stats['variants'],variants=selection['retained_sites'],
                 extra_variants=selection['extra_sites'],union_policy_id=policy['policy_id'],
                 single_original_vcf_export=True,direct_union_variant_selection=True,
                 scoring_variants=len(pd.read_csv(out/'scoring_variants.tsv.gz',sep='\t')),
                 variant_selection='Gnomix union SBayesRC', beagle_threads=cfg['beagle_threads'])
    write_json(out/'phasing.json',stats)
    complete(out,manifest['manifest_id'],{'stage':'phase','chrom':a.chrom,**stats})


def correct_scoring_genotypes(pgen,psam,mapping,swaps,out):
    import pgenlib
    samples=read_psam(psam);n,v=len(samples),len(mapping);dest=Path(out)/'scoring_genotypes_gnofix';expected=hashlib.sha256()
    with pgenlib.PgenReader(str(pgen).encode(),raw_sample_ct=n) as reader, \
         pgenlib.PgenWriter(str(dest.with_suffix('.pgen')).encode(),n,v,nonref_flags=False,allele_ct_limit=2,hardcall_phase_present=True) as writer:
        for start in range(0,v,1024):
            chunk=mapping.iloc[start:start+1024];a=np.empty((len(chunk),2*n),np.int32);phase=np.empty((len(chunk),n),np.uint8)
            reader.read_alleles_and_phasepresent_list(chunk.source_variant_index.to_numpy(np.uint32),a,phase)
            if (a<0).any() or (a>1).any() or ((a[:,::2]!=a[:,1::2]) & (phase==0)).any():raise ValueError('Scoring input is not fully phased')
            pairs=a.reshape(len(chunk),n,2);mask=swaps[:,chunk.window.to_numpy()].T
            pairs[mask]=pairs[mask,::-1];expected.update(a.tobytes());writer.append_alleles_batch(a,all_phased=True)
    actual=hashlib.sha256()
    with pgenlib.PgenReader(str(dest.with_suffix('.pgen')).encode(),raw_sample_ct=n) as reader:
        for start in range(0,v,1024):
            stop=min(v,start+1024);a=np.empty((stop-start,2*n),np.int32);phase=np.empty((stop-start,n),np.uint8)
            reader.read_alleles_and_phasepresent_range(start,stop,a,phase)
            if ((a[:,::2]!=a[:,1::2]) & (phase==0)).any():raise ValueError('Scoring output lost phase')
            actual.update(a.tobytes())
    if actual.digest()!=expected.digest():raise ValueError('Gnofix scoring readback differs')
    pd.DataFrame({'#CHROM':'chr'+mapping.chrom.astype(str),'POS':mapping.pos,'ID':mapping.id,'REF':mapping.ref_source,'ALT':mapping.alt_source}).to_csv(dest.with_suffix('.pvar'),sep='\t',index=False)
    samples[['FID','IID']].rename(columns={'FID':'#FID'}).to_csv(dest.with_suffix('.psam'),sep='\t',index=False)
    return {'samples':n,'variants':v,'all_dosages_preserved':True,'all_required_swaps_readback_verified':True,'phased_alleles_sha256':actual.hexdigest()}


def run_union_infer(a, cfg, manifest, policy):
    from phasing_and_gnomix.infer import run_infer
    run_infer(a, cfg, manifest)
    out = Path(a.output)
    previous = read_json(out/'COMPLETE.json')
    (out/'COMPLETE.json').unlink()
    samples = read_psam(a.psam)
    mapping = pd.read_csv(os.environ['SCORING_MAP'],sep='\t')
    win = pd.read_csv(a.windows,sep='\t')
    swaps = np.zeros((len(samples),len(win)),dtype=bool)
    seen = np.zeros(len(samples),np.uint8)
    for p in sorted((out/'parts').glob('*.npz')):
        start = int(p.stem.split('_')[1])
        with np.load(p) as part:
            end = start+len(part['samples'])
            if list(part['samples']) != list(samples.IID.iloc[start:end]):
                raise ValueError('Gnofix checkpoint sample order differs')
            swaps[start:end] = part['swaps']
            seen[start:end] += 1
    if not (seen == 1).all():
        raise ValueError('Missing or duplicate sample Gnofix corrections')
    report = correct_scoring_genotypes(a.pgen, a.psam, mapping, swaps, out)
    report.update(union_policy_id=policy['policy_id'], scoring_map_sha256=sha256(os.environ['SCORING_MAP']),
                  ancestry_labels=cfg['ancestries'], window_assignment='exact model window for matched sites; nearest reliably lifted model SNP for additional sites',
                  model_window_calls='parts/*.npz labels and called_posterior; index with scoring_variants.window',
                  sbayesrc_allele_orientation='PVAR retains source REF/ALT; scoring map records allele reversal relative to SBayesRC')
    write_json(out/'scoring_manifest.json',report)
    audit_outputs(a)
    detail = {k:v for k,v in previous.items() if k not in ['files','manifest_id']}
    detail.update(union_policy_id=policy['policy_id'], scoring_gnofix=report, independent_tract_audit=True)
    complete(out,manifest['manifest_id'],detail)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage',choices=['phase','infer'])
    for name in ['runtime','config','manifest','output','samples']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--chrom',required=True,type=int)
    for name in ['source','psam','keep','meta','scratch','map','model','pgen','match','windows','resource-plan','resource-plan-sha256']:
        p.add_argument('--'+name)
    p.add_argument('--heap-gb',type=int)
    p.add_argument('--shard',type=int)
    a = p.parse_args()
    configure_runtime(a.runtime)
    cfg, manifest = read_json(a.config), read_json(a.manifest)
    policy = manifest['union_policy']
    if cfg != manifest['config'] or str(a.chrom) not in policy['chromosomes']:
        raise ValueError('Mixed union policy, frozen cohort or scientific config')
    from .common import fingerprint
    if fingerprint({k:v for k,v in manifest.items() if k!='manifest_id'}) != manifest['manifest_id']:
        raise ValueError('Manifest fingerprint differs')
    if sha256(a.samples) != manifest['cohort_sha256']:
        raise ValueError('Frozen cohort changed')
    if sha256(a.model) != manifest['input_files'][f'models/chr{a.chrom}.pkl']['sha256']:
        raise ValueError('Unverified pretrained model')
    if a.stage == 'phase':
        cfg = {**cfg, 'beagle_threads': policy['chromosomes'][str(a.chrom)]['threads']}
        if sha256(os.environ['SBAYES']) != policy['chromosomes'][str(a.chrom)]['sbayesrc_csv_sha256']:
            raise ValueError('SBayesRC input checksum differs')
        run_union_phase(a,cfg,manifest,policy)
    else:
        run_union_infer(a,cfg,manifest,policy)


if __name__ == '__main__':
    main()
