"""Exercise complete coordinator transitions with fake cloud I/O; no paid jobs.

Run with the notebook Python (Google client dependencies installed):
python -m phasing_and_gnomix.tests.exercise_scheduler
"""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch, Mock
import tempfile
import time

from phasing_and_gnomix.common import read_json, write_json
from phasing_and_gnomix import scheduler
from phasing_and_gnomix.settings import PACKAGE


def exercise(directory):
    cfg=read_json(PACKAGE/'config.json')
    cfg.update(project='synthetic-project',region='us-central1',service_account='fake@example.invalid',
               source_uri='gs://synthetic-source/pgen',network='fake-network',subnetwork='fake-subnet',shard_size=2)
    manifest={'manifest_id':'synthetic','run_id':'synthetic','config':cfg,
              'resources':read_json(PACKAGE/'resources.json'),'output_uri':'gs://synthetic-workspace/union',
              'cohort':{'cohort':3,'shared_psam':'acaf_threshold.chr2.psam'},'union_policy':{'policy_id':'union'},
              'input_files':{'code.tar.gz':{'sha256':'code'},'runtime.tar.gz':{'sha256':'runtime'}},
              'source_objects':{}}
    for c in range(1,23):
        for ext in ('pgen','pvar','psam'):
            manifest['source_objects'][cfg['source_uri']+f'/acaf_threshold.chr{c}.{ext}']={'generation':1}
    rd=directory/'runs/synthetic';rd.mkdir(parents=True)
    write_json(rd/'manifest.json',manifest)
    cloud={};submissions=[]
    class Store:
        def json(self,uri):return manifest
        def metadata(self,uri):return {'generation':1}
        def completed(self,uri,mid):
            state=read_json(rd/'state.json')
            job=next(j for j in state['jobs'].values() if j['output']==uri)
            stage,c,b=job['stage'],job['chrom'],job['shard']
            value={'manifest_id':'synthetic','stage':stage,'chrom':c,'shard':b,'union_policy_id':'union'}
            if stage=='smoke':value['public_end_to_end_inference_passed']=True
            elif stage=='phase':value.update(samples=3,direct_union_variant_selection=True,single_original_vcf_export=True,
                        beagle_threads=manifest['resources']['chromosomes'][str(c)]['threads'])
            else:
                n=min(2,3-2*b)
                value.update(samples=n,gnofix_completed=n,independent_tract_audit=True,
                             scoring_gnofix={'samples':n,'all_dosages_preserved':True,'all_required_swaps_readback_verified':True})
            return value
    class Writer:
        def __init__(self,gcs,path,uri):self.path=path
        def save(self,state):write_json(self.path,state)
    def submit(args,stdout,**kwargs):
        state=read_json(rd/'state.json')
        pending=[j for j in state['jobs'].values() if j['job_id'] is None]
        assert len(pending)==1 and pending[0]['state']=='SUBMITTING'
        unit=pending[0]['unit'];jid='fake'+str(len(submissions))
        submissions.append((unit,list(args)))
        stdout.write('Job properties:\n  job-id: '+jid+'\nLaunched job-id: '+jid+'\n')
        cloud[jid]={'labels':{'lai-unit':unit,'lai-run':'phasing-gnomix-synthetic','job-id':jid},
                    'status':{'state':'RUNNING'}}
        return SimpleNamespace(returncode=0)
    def succeed(stage,c,b=0):
        state=read_json(rd/'state.json');job=state['jobs'][f'{stage}:{c}:{b}']
        cloud[job['job_id']]['status']['state']='SUCCEEDED'
    args=SimpleNamespace(work=directory,run_id='synthetic',once=True)
    with patch.object(scheduler,'GCS',lambda p:Store()),patch.object(scheduler,'StateWriter',Writer), \
         patch.object(scheduler,'statuses',lambda cfg,label:cloud), \
         patch.object(scheduler,'verify_frozen',lambda *a:None),patch.object(scheduler.subprocess,'run',submit):
        scheduler.loop(args)
        assert [u for u,_ in submissions]==['smoke-c22-b0-a0']
        succeed('smoke',22);scheduler.loop(args)
        assert [u for u,_ in submissions][-1]=='phase-c22-b0-a0' and len(submissions)==2
        succeed('phase',22);scheduler.loop(args)
        assert len(submissions)==4 and not read_json(rd/'state.json')['pilot_complete']
        succeed('infer',22,0);scheduler.loop(args)
        assert len(submissions)==4 and not read_json(rd/'state.json')['pilot_complete']
        succeed('infer',22,1);scheduler.loop(args)
        assert len(submissions)==25 and read_json(rd/'state.json')['pilot_complete']
        state=read_json(rd/'state.json')
        assert set(j['chrom'] for j in state['jobs'].values() if j['stage']=='phase')==set(range(1,23))
        failed=state['jobs']['phase:1:0'];cloud[failed['job_id']]['status']['state']='FAILED'
        succeed('phase',2);scheduler.loop(args)
        state=read_json(rd/'state.json')
        assert 'phase:1:0' in state['failures'] and len(submissions)==27
        assert sum(unit.startswith('phase-c1-') for unit,_ in submissions)==1
        assert 'infer:2:0' in state['jobs'] and 'infer:2:1' in state['jobs']
    print('COORDINATOR TRANSITIONS PASSED: reservation, smoke, full chr22 pilot, all22 fan-out, downstream overlap, no paid failure retry')


def exercise_recovery(directory):
    from google.api_core.exceptions import NotFound, ServiceUnavailable
    manifest={'run_id':'synthetic','config':{'project':'synthetic-project','region':'us-central1'}}
    unit='infer-c22-b0-a0'
    (directory/(unit+'.submission.log')).write_text(
        'Job properties:\n  job-id: attempted\ngoogle.api_core.exceptions.ServiceUnavailable: 503 Bad Gateway\n')
    job={'unit':unit,'job_id':None,'submitted_at':time.time()-300,'state':'SUBMITTING'}
    client=Mock()
    with patch('google.cloud.batch_v1.BatchServiceClient',return_value=client):
        client.get_job.side_effect=ServiceUnavailable('temporary read error')
        try:scheduler.reconcile_pending(manifest,directory,job,{})
        except ServiceUnavailable:pass
        else:raise AssertionError('Read error was interpreted as absence')
        assert job['job_id'] is None and 'first_confirmed_absent' not in job
        client.get_job.side_effect=NotFound('authoritative404')
        assert not scheduler.reconcile_pending(manifest,directory,job,{})
        job['first_confirmed_absent']-=31
        assert scheduler.reconcile_pending(manifest,directory,job,{})
        # An accepted FAILED task is adopted even when directGET lags the list.
        cloud={'attempted':{'labels':{'job-id':'attempted','lai-unit':unit},'status':{'state':'FAILED'}}}
        assert not scheduler.reconcile_pending(manifest,directory,job,cloud)
        assert job['job_id']=='attempted' and 'first_confirmed_absent' not in job
    print('SUBMISSION RECOVERY PASSED: API read errors held, repeated404 required, accepted failed task adopted')


if __name__=='__main__':
    with tempfile.TemporaryDirectory(prefix='public-union-scheduler-') as temp:
        exercise(Path(temp))
        exercise_recovery(Path(temp))
