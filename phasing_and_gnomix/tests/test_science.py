from pathlib import Path
from types import SimpleNamespace
import gzip
import numpy as np
import pandas as pd
import pytest

from phasing_and_gnomix.cohort import select
from phasing_and_gnomix.common import derived_seed,complete,verify_complete,read_json,read_psam,run
from phasing_and_gnomix.extract import match,candidate_ids,passing_missingness,prepare,pvar_records
from phasing_and_gnomix.models import assert_snv_alleles
from phasing_and_gnomix.infer import reconstruct_haplotypes,model_input
from phasing_and_gnomix.phase import check_and_shard
from phasing_and_gnomix.normalize_pvar_alleles import right_trim


def test_cohort_union_and_determinism():
    p=pd.DataFrame({'FID':['0']*30,'IID':[f'i{i:02}' for i in range(30)]})
    eur=set(p.IID[:20]);fit=set(p.IID[:15])
    a=select(p,eur,fit,5,20261001);b=select(p,eur,set(reversed(sorted(fit))),5,20261001)
    pd.testing.assert_frame_equal(a,b)
    assert set(p.IID[20:])<=set(a.IID)
    assert len(a)==15 and (a.reason=='random_fit_pca').sum()==5
    assert a.source_index.is_monotonic_increasing
    with pytest.raises(ValueError):select(p,eur,fit|{'absent'},5,1)


def test_seed_independent_of_python_hash_and_batch():
    assert derived_seed(20261001,22,'sample-A')==derived_seed(20261001,22,'sample-A')
    assert len({derived_seed(20261001,c,s) for c in [1,22] for s in ['a','b']})==4


def test_original_allele_audit_before_truncation():
    assert_snv_alleles(['A','G'],['C','T'])
    for r,a in [(['AT'],['A']),(['N'],['G']),(['A'],['A'])]:
        with pytest.raises(ValueError):assert_snv_alleles(r,a)


def test_normalize_before_snv_filter(tmp_path):
    p=tmp_path/'x.pvar'
    p.write_text('#CHROM\tPOS\tID\tREF\tALT\nchr22\t100\ta\tTATG\tCATG,T\nchr22\t200\tb\tA\tAT\n')
    ids=tmp_path/'ids'
    assert candidate_ids(p,{'pos38':np.array([100,200])},ids)==1
    assert ids.read_text()=='a\n'
    assert right_trim('TATG','CATG')==('T','C')


def test_position_is_not_the_allele_matching_key():
    meta={'pos38':np.array([100,100,200,200,300,-1]),'ref38':np.array(['A','A','C','T','G','A']),
          'alt38':np.array(['G','T','T','C','A','C'])}
    panel=pd.DataFrame({'pidx':[0,1,2,3],'pos':[100,100,200,300],'id':['a','b','c','d'],
                        'pref':['A','A','C','A'],'palt':['G','T','T','G']})
    m,s=match(meta,panel)
    assert list(m.i)==[0,1,4]
    assert list(m.flip)==[False,False,True]
    assert list(s)==['matched','matched','ambiguous','ambiguous','matched','failed_liftover']


def test_missingness_exact_boundary():
    assert list(passing_missingness([0,1,2,0],[10,10,10,0],.1))==[True,True,False,False]
    assert list(passing_missingness([12345,12346],[123459,123459],.1))==[True,False]


def test_remainder_uses_true_model_window():
    X=np.array([[0]*11,[1]*11],np.int8)
    s=np.array([[False,True,False]])
    Y=reconstruct_haplotypes(X,s,11,3,3)
    assert list(Y[0])==[0,0,0,1,1,1,0,0,0,0,0]
    np.testing.assert_array_equal(Y[0]+Y[1],X[0]+X[1])


def test_completion_rejects_truncation_and_mixed_run(tmp_path):
    (tmp_path/'data').write_text('abcdef')
    complete(tmp_path,'id',{'samples':42})
    assert verify_complete(tmp_path,'id')
    with pytest.raises(ValueError):verify_complete(tmp_path,'other')
    (tmp_path/'data').write_text('abc')
    assert not verify_complete(tmp_path,'id')
