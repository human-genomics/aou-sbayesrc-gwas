import copy
import io
import json
from pathlib import Path
import tarfile
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest

from phasing_and_gnomix.common import read_json, write_json
from phasing_and_gnomix.prepare import freeze_input, source_archive
from phasing_and_gnomix.scheduler import command, eligible_chromosomes, pending_decision, require_marker
from phasing_and_gnomix.settings import PACKAGE, REPO, mounted_workspace, resolve, source_files


def manifest_fixture():
    cfg=read_json(PACKAGE/'config.json')
    cfg.update(seeds=read_json(PACKAGE/'seeds.json'),project='synthetic-project',region='us-central1',
               source_uri='gs://synthetic-source/acaf_threshold/pgen',service_account='synthetic@example.invalid',
               network='projects/synthetic-project/global/networks/network',subnetwork='synthetic-subnet',shard_size=2)
    resources=read_json(PACKAGE/'resources.json')
    return {'manifest_id':'synthetic','run_id':'synthetic-union','config':cfg,'resources':resources,
            'output_uri':'gs://synthetic-workspace/union','cohort':{'cohort':3,'shared_psam':'acaf_threshold.chr2.psam'},
            'union_policy':{'policy_id':'union'}, 'input_files':{'code.tar.gz':{'sha256':'code'},'runtime.tar.gz':{'sha256':'runtime'}}}


def marker_fixture(m,stage,c,b=0):
    value={'manifest_id':m['manifest_id'],'stage':stage,'chrom':c,'shard':b,'union_policy_id':'union',
           'samples':3 if stage=='phase' else min(2,3-2*b)}
    if stage=='phase':
        value.update(direct_union_variant_selection=True,single_original_vcf_export=True,
                     beagle_threads=m['resources']['chromosomes'][str(c)]['threads'])
    elif stage=='infer':
        value.update(gnofix_completed=value['samples'],independent_tract_audit=True,
                     scoring_gnofix={'samples':value['samples'],'all_dosages_preserved':True,
                                     'all_required_swaps_readback_verified':True})
    else:value['public_end_to_end_inference_passed']=True
    return value


def test_worker_commands_use_union_and_correct_phase(tmp_path):
    m=manifest_fixture()
    for c in range(1,23):
        args,out,unit=command(m,'phase',c)
        assert '--use-private-address' in args and '--unique-job-id' in args
        assert args[args.index('--machine-type')+1]==m['resources']['chromosomes'][str(c)]['machine_type']
        assert f'SBAYES={m["output_uri"]}/inputs/sbayesrc/chr{c}.csv' in args
        assert f'PSAM={m["config"]["source_uri"]}/acaf_threshold.chr2.psam' in args
        assert out.endswith(f'/phased/chr{c}')
        infer,iout,_=command(m,'infer',c,1)
        assert f'PGEN={out}/shards/batch_00001.pgen' in infer
        assert f'SCORING_MAP={out}/scoring_variants.tsv.gz' in infer
        assert iout.endswith(f'/inference/chr{c}/batch_00001')


def test_pilot_requires_full_inference_before_fanout():
    assert eligible_chromosomes({})==[]
    assert eligible_chromosomes({'smoke_done':True})==[22]
    assert eligible_chromosomes({'smoke_done':True,'pilot_complete':True})==list(range(22,0,-1))


def test_no_duplicate_or_uncertain_submissions():
    job={'unit':'unit'}
    row={'labels':{'lai-unit':'unit','job-id':'attempted'},'status':{'state':'FAILED'}}
    assert pending_decision(job,{},row,300,True,'attempted')==('adopt','attempted')
    assert pending_decision(job,{'attempted':row},None,300,True,'attempted')==('adopt','attempted')
    assert pending_decision(job,{},None,10,True,'attempted')==('wait',None)
    assert pending_decision(job,{},None,300,True,'attempted')==('absent',None)
    for cloud,direct,transient,jid in [({'attempted':row,'duplicate':row},None,True,'attempted'),
                                      ({},None,False,'attempted'),({},None,True,None)]:
        with pytest.raises(ValueError):pending_decision(job,cloud,direct,300,transient,jid)


def test_completion_requires_every_person_and_corrected_scoring_phase():
    m=manifest_fixture()
    marker=marker_fixture(m,'infer',22,0)
    require_marker(marker,m,'infer',22,0)
    for path,value in [('samples',2+1),('gnofix_completed',0),('union_policy_id','old-model-only'),
                       ('independent_tract_audit',False)]:
        bad=copy.deepcopy(marker);bad[path]=value
        with pytest.raises(ValueError):require_marker(bad,m,'infer',22,0)
    bad=copy.deepcopy(marker);bad['scoring_gnofix']['all_required_swaps_readback_verified']=False
    with pytest.raises(ValueError):require_marker(bad,m,'infer',22,0)


def test_source_archive_is_deterministic_and_excludes_working_data(tmp_path):
    files=source_files()
    assert all(p.suffix in {'.py','.sh','.json','.txt','.md'} or p.name=='.gitignore' for p in files)
    assert not any('data' in p.relative_to(PACKAGE).parts for p in files)
    first,second=tmp_path/'first.tar.gz',tmp_path/'second.tar.gz'
    source_archive(first,files);source_archive(second,list(reversed(files)))
    assert first.read_bytes()==second.read_bytes()
    with tarfile.open(first) as archive:
        assert set(archive.getnames())=={str(p.relative_to(REPO)) for p in files}


def test_frozen_cache_survives_later_preparation(tmp_path):
    source=tmp_path/'mutable';source.write_text('first input')
    first=freeze_input(source,tmp_path/'work')
    source.write_text('second input')
    second=freeze_input(source,tmp_path/'work')
    assert first['local_path']!=second['local_path']
    assert Path(first['local_path']).read_text()=='first input'
    assert freeze_input(source,tmp_path/'work')==second


def test_preparation_lock_prevents_cache_mutation(monkeypatch,tmp_path):
    import fcntl
    import phasing_and_gnomix.prepare as prepare
    monkeypatch.setattr(prepare,'resolve',lambda args: ({},{},'gs://synthetic/run'))
    entered=[]
    monkeypatch.setattr(prepare,'prepare_locked',lambda *args:entered.append(True))
    args=SimpleNamespace(work=tmp_path)
    with (tmp_path/'prepare.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):prepare.prepare(args)
        assert not entered
    prepare.prepare(args)
    assert entered==[True]


def test_reference_cache_requires_success_and_refreshes_unverified_urls(monkeypatch,tmp_path):
    from phasing_and_gnomix.prepare import require_public_checks
    from phasing_and_gnomix import bootstrap
    from phasing_and_gnomix.common import sha256
    (tmp_path/'references').mkdir()
    write_json(tmp_path/'references/demo_validation.json',{})
    with pytest.raises(ValueError):require_public_checks(tmp_path)
    target=tmp_path/'chr22.gmap.gz';target.write_bytes(b'old map commit')
    downloads=[]
    def fetch(args):
        downloads.append(args)
        Path(args[-1]).write_bytes(b'new map commit')
    monkeypatch.setattr(bootstrap,'run',fetch)
    bootstrap.fetch('https://example.invalid/pinned-map',target)
    assert target.read_bytes()==b'new map commit' and len(downloads)==1
    bootstrap.fetch('https://example.invalid/pinned-map',target,sha256(target))
    assert len(downloads)==1


def test_workspace_discovery_and_custom_settings(monkeypatch,tmp_path):
    import phasing_and_gnomix.settings as settings
    monkeypatch.setenv('WORKSPACE_BUCKET','gs://stale-cloned-bucket')
    monkeypatch.setattr(settings,'mounted_workspace',lambda:'gs://mounted-workspace')
    monkeypatch.setattr(settings,'gcloud_value',lambda key:'synthetic@example.invalid' if key=='account' else 'synthetic-project')
    args=SimpleNamespace(run_id='test',config=PACKAGE/'config.json',seeds=PACKAGE/'seeds.json',
                         resources=PACKAGE/'resources.json',work=tmp_path,region='us-central1',
                         project='synthetic-project',workspace_bucket=None,pipeline_uri=None,data_version='v9',
                         source_uri='gs://test/source',output_uri=None,service_account=None,network=None,subnetwork=None)
    cfg,resources,output=resolve(args)
    assert cfg['pipeline_uri']=='gs://mounted-workspace/sbayesrc_genotypes'
    assert output.endswith('/phasing_and_gnomix/test')
    args.workspace_bucket='gs://explicit-workspace'
    assert resolve(args)[0]['pipeline_uri'].startswith('gs://explicit-workspace/')
    args.run_id='../invalid'
    with pytest.raises(ValueError):resolve(args)
    assert mounted_workspace('bucket /home/jupyter/workspace/workspace-bucket fuse.gcsfuse rw 0 0')=='gs://bucket'


def test_all22_aggregation_rejects_missing_and_mixed_outputs(tmp_path):
    from phasing_and_gnomix.results import aggregate
    m=manifest_fixture();pops=m['config']['ancestries']
    (tmp_path/'cohort').mkdir()
    pd.DataFrame({'FID':['0']*3,'IID':['synthetic0','synthetic1','synthetic2'],
                  'reason':['synthetic']*3,'sample_index':range(3)}).to_csv(tmp_path/'cohort/samples.tsv',sep='\t',index=False)
    class Store:
        def __init__(self):self.data={};self.missing=True;self.mixed=False
        def completed(self,uri,mid):
            if uri.endswith('/results'):
                return json.loads(self.data[uri+'/COMPLETE.json']) if uri+'/COMPLETE.json' in self.data else None
            c=int(uri.split('/chr')[1].split('/')[0]);b=int(uri.split('/batch_')[1]) if '/batch_' in uri else 0
            if self.missing and c==22:return None
            value=marker_fixture(m,'infer' if '/inference/' in uri else 'phase',c,b)
            if self.mixed and c==22:value['union_policy_id']='old-model-only'
            return value
        def blob(self,uri):
            b=int(uri.split('/batch_')[1].split('/')[0]);c=int(uri.split('/chr')[1].split('/')[0])
            samples=np.array(['synthetic0','synthetic1'] if b==0 else ['synthetic2'])
            totals=np.zeros((len(samples),8));totals[:,0]=c;totals[:,7]=2*c
            buf=io.BytesIO();np.savez(buf,samples=samples,populations=np.array(pops),totals=totals)
            return SimpleNamespace(download_as_bytes=lambda:buf.getvalue())
        def put(self,path,uri):self.data[uri]=Path(path).read_bytes()
    store=Store()
    with pytest.raises(ValueError):aggregate(m,tmp_path,store)
    assert not store.data
    store.missing=False;store.mixed=True
    with pytest.raises(ValueError):aggregate(m,tmp_path,store)
    assert not store.data
    store.mixed=False;aggregate(m,tmp_path,store)
    result=pd.read_csv(tmp_path/'results/global_ancestry_gnofix.tsv.gz',sep='\t')
    np.testing.assert_allclose(result.EUR,1/3);np.testing.assert_allclose(result.WAS,2/3)
    np.testing.assert_allclose(result[pops].sum(axis=1),1)
    assert (result.completed_chromosomes==22).all()
    assert len(pd.read_csv(tmp_path/'results/scoring_file_index.tsv.gz',sep='\t'))==44
    before=dict(store.data);aggregate(m,tmp_path,store);assert store.data==before
    # Reuse the verified byte-identical local spool after a partial cloud upload.
    store.data.pop(m['output_uri']+'/results/COMPLETE.json');aggregate(m,tmp_path,store)
    assert store.data==before
