"""Public chr22 regression checks and AoU missing-feature sensitivity."""
from __future__ import annotations
import argparse
from pathlib import Path
import tempfile
from types import SimpleNamespace
import numpy as np
import pandas as pd
from .common import configure_runtime,read_json,write_json,derived_seed
from .models import load_model
from .infer import calls,model_input


def demo(runtime,model_path,config,output,match_path=None):
    configure_runtime(runtime)
    from src.utils import read_vcf,vcf_to_npy
    from gnomix1000g.infer import gnomix_calls,window_of_snp
    model=load_model(model_path)
    vcf=read_vcf(str(Path(runtime)/'gnomix/demo/data/small_query_chr22.vcf.gz'),chm='22',fields='*')
    X=vcf_to_npy(vcf,model.snp_pos,model.snp_ref,verbose=False)
    # Upstream's position-only intersection leaves a few duplicated-position
    # features at code 2 even in its "complete" demo. Report this explicitly.
    upstream_missing=int(np.count_nonzero(~np.isin(X,[0,1])))
    official=model.predict_proba(X).argmax(axis=-1)
    X[~np.isin(X,[0,1])]=config['model_fill']
    baseline=model.predict_proba(X).argmax(axis=-1)
    rows=[]
    rng=np.random.Generator(np.random.PCG64(config['seeds']['validation']))
    masks={'random_5pct':rng.random(model.C)<.05,'random_20pct':rng.random(model.C)<.20}
    if match_path:
        absent=np.ones(model.C,bool)
        absent[np.load(match_path)['model_idx']]=False
        masks['actual_aou_mask']=absent
    for name,mask in masks.items():
        for fill in (0,2):
            x=X.copy();x[:,mask]=fill
            labels=model.predict_proba(x).argmax(axis=-1)
            rows.append({'mask':name,'fill':fill,'fraction_masked':float(mask.mean()),
                         'window_agreement':float((labels==baseline).mean())})
    labels,support,swaps,_=calls(model,X,list(vcf['samples']),22,config['seeds']['inference'])
    reference,_=gnomix_calls(model,X,window_of_snp(model.C,model.M,model.W))
    if not np.array_equal(labels,reference['lab_fix']) or not np.array_equal(swaps,reference['swap']):
        raise ValueError('Our Gnofix path differs from the pinned 1000G adapter')
    expected_support=np.take_along_axis(reference['prob_fix'],labels[...,None].astype(np.int64),axis=2)[...,0]/255
    if np.max(np.abs(expected_support-support))>1/255:
        raise ValueError('Posterior difference exceeds reference quantization precision')
    # Confirm batch boundaries do not alter corrected labels or phase.
    midpoint=len(vcf['samples'])//2
    split=[calls(model,X[:2*midpoint],list(vcf['samples'][:midpoint]),22,config['seeds']['inference']),
           calls(model,X[2*midpoint:],list(vcf['samples'][midpoint:]),22,config['seeds']['inference'])]
    if not np.array_equal(labels,np.concatenate([x[0] for x in split])) or not np.array_equal(swaps,np.concatenate([x[2] for x in split])):
        raise ValueError('Inference changes across batch boundaries')
    report={'gnofix_reference_labels_equal':True,'gnofix_reference_swaps_equal':True,'batch_invariant':True,
            'upstream_demo_missing_alleles':upstream_missing,
            'reference_fill_vs_official_demo_window_agreement':float((baseline==official).mean()),
            'seed':config['seeds']['validation'],'missing_feature_sensitivity':rows}
    write_json(output,report)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('runtime','model','config','output'):p.add_argument('--'+n,required=True)
    p.add_argument('--match')
    a=p.parse_args()
    demo(a.runtime,a.model,read_json(a.config),a.output,a.match)
