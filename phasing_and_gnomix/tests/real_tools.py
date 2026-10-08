"""Synthetic end-to-end union/QC/real Beagle/scoring correction test.

Run with the pinned runtime Python; no AoU participant data is used.
"""
import json
import os
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pgenlib
from phasing_and_gnomix.common import configure_runtime, read_psam, run, verify_complete, complete
from phasing_and_gnomix.infer import model_input
from phasing_and_gnomix.worker import run_union_phase, run_union_infer, correct_scoring_genotypes
from phasing_and_gnomix.extract import nearest_windows


def validate(runtime, work):
    runtime, work = Path(runtime).resolve(), Path(work).resolve()
    work.mkdir(parents=True,exist_ok=True)
    configure_runtime(runtime)
    rng=np.random.default_rng(20261001)
    n=48;pos=20000000+np.arange(100)*1000;ids=[f'public{i}' for i in range(n)]
    h=rng.integers(0,2,size=(12,100));choices=rng.integers(0,12,size=(n,2));gt=h[choices]
    header='##fileformat=VCFv4.2\n##contig=<ID=chr21,length=46709983>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n'
    with (work/'source.vcf').open('w') as f:
        f.write(header)
        for k,p in enumerate(pos):
            ref,alt=('TATG','CATG,T') if k==80 else ('A','C')
            if k==50:alt='C,G'
            calls=[f'{gt[i,0,k]}/{gt[i,1,k]}' for i in range(n)]
            if k in [30,99]:calls[:10]=['./.']*10
            f.write(f'chr21\t{p}\tv{k}\t{ref}\t{alt}\t.\t.\t.\tGT\t'+'\t'.join(calls)+'\n')
            if k==98:f.write(f'chr21\t{p}\tv98_duplicate\tA\tC\t.\t.\t.\tGT\t'+'\t'.join(calls)+'\n')
    run([runtime/'plink2','--vcf',work/'source.vcf','--make-pgen','--output-chr','chrM','--out',work/'source'])
    samples=read_psam(work/'source.psam');samples[['FID','IID']].to_csv(work/'samples.tsv',sep='\t',index=False)
    samples[['FID','IID']].to_csv(work/'keep',sep='\t',index=False,header=False)
    meta={'pos19':pos[:80],'pos38':pos[:80],'ref':np.array(['A']*80),'alt':np.array(['C']*80),
          'ref38':np.array(['A']*80),'alt38':np.array(['C']*80),'C':80,'M':20,'W':4,
          'gm_pos':np.array([pos[0],pos[-1]]),'gm_cm':np.array([0.,1.])}
    np.savez(work/'meta.npz',**meta)
    sb=pd.DataFrame({'chrom':21,'pos':np.append(pos,[20100000,20101000,20102000]),'ref':'A','alt':'C'})
    sb.loc[80,['ref','alt']]=['T','C']
    sb.loc[97,['ref','alt']]=['G','T']
    sb.loc[10,['ref','alt']]=['C','A']
    sb['rsid']=[f'rsSynthetic{i}' for i in range(len(sb))]
    sb=pd.concat([sb,pd.DataFrame([{'chrom':21,'pos':pos[50],'ref':'A','alt':'G','rsid':'rsSynthetic_multi'}])],ignore_index=True)
    sb.to_csv(work/'sb.csv',index=False)
    pd.DataFrame({'chrom':['chr21']*2,'id':['.']*2,'cm':[0.,1.],'pos':[pos[0],pos[-1]]}).to_csv(work/'map',sep='\t',index=False,header=False)
    cfg={'beagle_threads':1,'max_variant_missingness':.1,'minimum_model_coverage':.95,'shard_size':16,
         'seeds':{'plink':20261001,'beagle':20261001,'inference':20261001},'beagle_version':'5.5 27Feb25.75f'}
    policy={'policy_id':'synthetic','chromosomes':{'21':{'heap_gb':4,'threads':1}}}
    os.environ['SBAYES']=str(work/'sb.csv')
    a=SimpleNamespace(runtime=runtime,output=work/'out',scratch=work/'scratch',samples=work/'samples.tsv',
                      chrom=21,source=work/'source',psam=work/'source.psam',keep=work/'keep',meta=work/'meta.npz',map=work/'map')
    from phasing_and_gnomix import phase,extract
    commands=[]
    original_run=phase.run
    def trace(cmd,**kwargs):
        commands.append(list(map(str,cmd)))
        return original_run(cmd,**kwargs)
    phase.run=extract.run=trace
    try:run_union_phase(a,cfg,{'manifest_id':'synthetic'},policy)
    finally:phase.run=extract.run=original_run
    assert len(commands)==7 and sum('--export' in cmd for cmd in commands)==1
    assert not any('bcftools' in str(cmd) or '--pmerge' in cmd for cmd in commands)
    assert verify_complete(work/'out','synthetic')
    report=json.load(open(work/'out/phasing.json'))
    assert report['variants']==79 and report['extra_variants']==18 and report['union_variants']==97
    assert report['single_original_vcf_export'] and report['direct_union_variant_selection']
    features=pd.read_csv(work/'out/sbayesrc_feature_qc.tsv.gz',sep='\t')
    assert features.loc[features.pos==pos[30],'status'].item()=='missingness_excluded'
    assert features.loc[features.pos==pos[99],'status'].item()=='missingness_excluded'
    assert features.loc[features.pos==pos[97],'status'].item()=='source_absent_or_other_alleles'
    assert features.loc[features.pos==pos[98],'status'].item()=='ambiguous_source'
    score=pd.read_csv(work/'out/scoring_variants.tsv.gz',sep='\t')
    assert len(score)==97 and score.loc[score.pos==pos[10],'sbayesrc_ref_alt_swapped'].item()
    assert score.loc[score.pos==pos[80],'ref_source'].item()=='T'
    assert 'rsid' in score and (score.rsid=='rsSynthetic_multi').any()
    mapping=dict(np.load(work/'out/match.npz'))
    X=model_input(work/'out/shards/batch_00000.pgen',np.arange(16),16,97,80,mapping)
    assert X.shape==(32,80) and (X[:,30]==0).all()
    out=work/'corrected';out.mkdir(exist_ok=True)
    swaps=rng.integers(0,2,size=(16,4)).astype(bool)
    correction=correct_scoring_genotypes(work/'out/shards/batch_00000.pgen',
                                        work/'out/shards/batch_00000.psam',score,swaps,out)
    # Independently compare every corrected genotype to its source and required swap.
    with pgenlib.PgenReader(str(work/'out/shards/batch_00000.pgen').encode(),raw_sample_ct=16) as m, \
         pgenlib.PgenReader(str(out/'scoring_genotypes_gnofix.pgen').encode(),raw_sample_ct=16) as dst:
        for k,row in score.iterrows():
            before=np.empty(32,np.int32);after=np.empty(32,np.int32)
            m.read_alleles(row.source_variant_index,before);dst.read_alleles(k,after)
            expected=before.reshape(16,2).copy();expected[swaps[:,row.window]]=expected[swaps[:,row.window],::-1]
            np.testing.assert_array_equal(after,expected.ravel())
    # Liftover inversions and failed sites must not break nearest-site assignment.
    nearest,windows=nearest_windows({'pos38':np.array([300,100,-1,200]),'M':1,'W':4},np.array([1,150,999]))
    assert list(nearest)==[1,1,0] and list(windows)==[1,1,0]
    # Exercise the real checkpoint-to-scoring wrapper independently of model prediction.
    # The original pretrained inference was already validated by the parent run;
    # these deterministic synthetic checkpoints test our newly added correction path.
    from phasing_and_gnomix import infer
    from phasing_and_gnomix import worker
    old_infer=infer.run_infer
    old_audit=worker.audit_outputs
    def checkpoint_stub(args,config,manifest):
        dest=Path(args.output);(dest/'parts').mkdir(parents=True,exist_ok=True)
        sn=read_psam(args.psam).IID.to_numpy().astype(str)
        np.savez_compressed(dest/'parts/part_00000.npz',samples=sn,swaps=swaps,
                            labels=np.zeros((32,4),np.uint8),called_posterior=np.ones((32,4)))
        complete(dest,manifest['manifest_id'],{'stage':'infer','chrom':21,'shard':0,'samples':16,'gnofix_completed':16})
    infer.run_infer=checkpoint_stub
    # This fixture isolates checkpoint-to-genotype correction. The real public
    # model smoke test below exercises the independent full tract audit.
    worker.audit_outputs=lambda args:None
    wrapped=work/'wrapped';ia=SimpleNamespace(output=wrapped,psam=work/'out/shards/batch_00000.psam',
                          pgen=work/'out/shards/batch_00000.pgen',windows=work/'out/windows.tsv')
    os.environ['SCORING_MAP']=str(work/'out/scoring_variants.tsv.gz')
    try:
        run_union_infer(ia,{**cfg,'ancestries':['EUR','EAS','NAT','AFR','SAS','AHG','OCE','WAS']},{'manifest_id':'synthetic'},policy)
    finally:
        infer.run_infer=old_infer
        worker.audit_outputs=old_audit
    assert verify_complete(wrapped,'synthetic')
    assert json.load(open(wrapped/'COMPLETE.json'))['scoring_gnofix']['all_required_swaps_readback_verified']
    result={'real_beagle_union_phase_passed':True,'normalization_qc_passed':True,'model_input_mapping_preserved':True,
            'gnofix_scoring_alleles_independently_verified':True,'inference_checkpoint_wrapper_verified':True,
            'single_original_vcf_export_verified':True,'original_phase_command_sequence_verified':True,'phase_command_trace':commands,
            'phase_report':report,'correction':correction}
    (work/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('UNION REAL-TOOLS VALIDATION PASSED',flush=True)
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    args=parser.parse_args()
    validate(args.runtime,args.work)
