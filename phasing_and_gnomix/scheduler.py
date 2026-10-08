"""Detached, resumable Batch coordination with an end-to-end chr22 pilot."""
import fcntl
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
import sys
import time

from .cloud import GCS, StatusLookupError, statuses
from .common import log, read_json, write_json
from .prepare import verify_frozen
from .settings import PACKAGE, REPO
from .state import StateWriter

TRANSIENT = re.compile(r'google\.api_core\.exceptions\.(?:ServiceUnavailable|DeadlineExceeded|InternalServerError|TooManyRequests|BadGateway|GatewayTimeout):[^\n]*\s*$')


def run_label(manifest):
    return 'phasing-gnomix-' + manifest['run_id']


def command(manifest, stage, chrom=22, shard=0, attempt=0):
    cfg, resources = manifest['config'], manifest['resources']
    root = manifest['output_uri']
    inputs = root + '/inputs/'
    unit = f'{stage}-c{chrom}-b{shard}-a{attempt}'
    files = {name: inputs + path for name, path in {
        'RUNTIME': 'runtime.tar.gz', 'CODE': 'code.tar.gz', 'CONFIG': 'config.json',
        'MANIFEST': 'manifest.json', 'SAMPLES': 'cohort/samples.tsv',
        'MODEL': f'models/chr{chrom}.pkl'}.items()}
    env = {'STAGE': stage, 'CHROM': str(chrom),
           'CODE_SHA256': manifest['input_files']['code.tar.gz']['sha256'],
           'RUNTIME_SHA256': manifest['input_files']['runtime.tar.gz']['sha256']}
    allocation = resources['chromosomes'][str(chrom)] if stage == 'phase' else resources['inference']
    phase = root + f'/phased/chr{chrom}'
    if stage == 'phase':
        files.update(PGEN=cfg['source_uri']+f'/acaf_threshold.chr{chrom}.pgen',
                     PVAR=cfg['source_uri']+f'/acaf_threshold.chr{chrom}.pvar',
                     PSAM=cfg['source_uri']+'/'+manifest['cohort']['shared_psam'],
                     KEEP=inputs+'cohort/cohort.keep', META=inputs+f'models/chr{chrom}.meta.npz',
                     MAP=inputs+f'maps/chr{chrom}.map', SBAYES=inputs+f'sbayesrc/chr{chrom}.csv')
        env['HEAP'] = str(allocation['heap_gb'])
        output, disk, hours = phase, resources['phase_disk_gb'], resources['phase_timeout_hours']
    elif stage == 'infer':
        files.update(PGEN=phase+f'/shards/batch_{shard:05d}.pgen', PSAM=phase+f'/shards/batch_{shard:05d}.psam',
                     MATCH=phase+'/match.npz', WINDOWS=phase+'/windows.tsv', SCORING_MAP=phase+'/scoring_variants.tsv.gz')
        env['SHARD'] = str(shard)
        output, disk, hours = root+f'/inference/chr{chrom}/batch_{shard:05d}', allocation['disk_gb'], allocation['timeout_hours']
    elif stage == 'smoke':
        files['META'] = inputs+'models/chr22.meta.npz'
        output, disk, hours = root+'/preflight/cloud_smoke', allocation['disk_gb'], allocation['timeout_hours']
    else:
        raise ValueError('Unknown cloud stage')
    args = [sys.executable, '-u', str(PACKAGE/'dsub_launcher.py'), '--provider', 'google-batch',
            '--project', cfg['project'], '--location', cfg['region'], '--regions', cfg['region'],
            '--service-account', cfg['service_account'], '--use-private-address',
            '--network', cfg['network'], '--subnetwork', cfg['subnetwork'], '--user-project', cfg['project'],
            '--logging', root+'/logs/dsub', '--name', 'lai-'+unit, '--unique-job-id',
            '--image', cfg['image'], '--script', str(PACKAGE/('cloud_smoke.sh' if stage == 'smoke' else 'worker.sh')),
            '--machine-type', allocation['machine_type'], '--disk-size', str(disk), '--disk-type', 'pd-ssd',
            '--boot-disk-size', '50', '--timeout', str(hours)+'h',
            '--label', 'lai-run='+run_label(manifest), 'lai-unit='+unit,
            '--env', *[f'{k}={v}' for k,v in env.items()],
            '--input', *[f'{k}={v}' for k,v in files.items()], '--output-recursive', 'OUTDIR='+output]
    return args, output, unit


def require_marker(marker, manifest, stage, chrom, shard=0):
    if marker is None or marker.get('manifest_id') != manifest['manifest_id'] or marker.get('stage') != stage:
        raise ValueError('Missing or incorrect completion marker')
    if stage == 'smoke':
        if not marker.get('public_end_to_end_inference_passed'):
            raise ValueError('Cloud public-data smoke test failed')
        return
    cfg = manifest['config']
    n = manifest['cohort']['cohort']
    expected = n if stage == 'phase' else min(cfg['shard_size'], n-shard*cfg['shard_size'])
    if (marker.get('chrom') != chrom or marker.get('samples') != expected
            or marker.get('union_policy_id') != manifest['union_policy']['policy_id']):
        raise ValueError('Mixed cohort, chromosome or variant-set output')
    if stage == 'phase':
        if (not marker.get('direct_union_variant_selection') or not marker.get('single_original_vcf_export')
                or marker.get('beagle_threads') != manifest['resources']['chromosomes'][str(chrom)]['threads']):
            raise ValueError('Phase output differs from the frozen union workflow')
    else:
        score = marker.get('scoring_gnofix', {})
        if (marker.get('shard') != shard or marker.get('gnofix_completed') != expected
                or score.get('samples') != expected or not score.get('all_required_swaps_readback_verified')
                or not score.get('all_dosages_preserved') or not marker.get('independent_tract_audit')):
            raise ValueError('Incomplete Gnofix, scoring genotype or tract validation')


def pending_decision(job, cloud, direct, age, transient, attempted_id):
    matches = {jid: row for jid, row in cloud.items() if row['labels'].get('lai-unit') == job['unit']}
    if direct is not None:
        labels = direct['labels']
        if labels.get('lai-unit') != job['unit'] or labels.get('job-id') != attempted_id:
            raise ValueError('Direct Batch identity mismatch')
        matches[attempted_id] = direct
    if len(matches) > 1:
        raise ValueError('Duplicate accepted cloud requests for one unit')
    if matches:
        jid = next(iter(matches))
        if attempted_id and jid != attempted_id:
            raise ValueError('Unexpected accepted ID for a reserved unit')
        return 'adopt', jid
    if age < 180:
        return 'wait', None
    if not attempted_id or not transient:
        raise ValueError('Unaccepted submission needs review; see its submission log')
    return 'absent', None


def reconcile_pending(manifest, directory, job, cloud):
    from google.cloud import batch_v1
    from google.api_core.exceptions import NotFound
    path = directory / (job['unit']+'.submission.log')
    content = path.read_text() if path.exists() else ''
    found = re.findall(r'^\s*job-id:\s*([a-z0-9-]+)\s*$', content, re.M)
    if len(found) > 1:
        raise ValueError('Ambiguous attempted Batch IDs')
    attempted = found[0] if found else None
    direct = None
    if attempted:
        cfg = manifest['config']
        client = batch_v1.BatchServiceClient(transport='rest')
        try:
            try:
                value = client.get_job(name=f"projects/{cfg['project']}/locations/{cfg['region']}/jobs/{attempted}-0-0",
                                       timeout=30, retry=None)
                if value.labels.get('lai-run') != run_label(manifest):
                    raise ValueError('Attempted ID belongs to another run')
                direct = {'labels': dict(value.labels), 'status': {'state': value.status.state.name}}
            except NotFound:
                pass  # Only an authoritative404 can establish absence.
        finally:
            client.transport.close()
    action, jid = pending_decision(job, cloud, direct, time.time()-job['submitted_at'],
                                   bool(TRANSIENT.search(content)), attempted)
    if action == 'adopt':
        job.update(job_id=jid, state='SUBMITTED')
        job.pop('first_confirmed_absent', None)
    elif action == 'absent':
        first = job.setdefault('first_confirmed_absent', time.time())
        if time.time()-first >= 30:
            return True
    else:
        job.pop('first_confirmed_absent', None)
    return False


def eligible_chromosomes(state):
    if not state.get('smoke_done'):
        return []
    return list(range(22, 0, -1)) if state.get('pilot_complete') else [22]


def loop(args):
    directory = args.work.resolve()/'runs'/args.run_id
    manifest = read_json(directory/'manifest.json')
    cfg = manifest['config']
    gcs = GCS(cfg['project'])
    lock = (directory/'coordinator.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    verify_frozen(manifest, directory, gcs)
    if gcs.json(manifest['output_uri']+'/inputs/manifest.json') != manifest:
        raise ValueError('Run prepare to finish staging all immutable inputs')
    path = directory/'state.json'
    state = read_json(path) if path.exists() else {'manifest_id': manifest['manifest_id'], 'jobs': {},
              'phase_done': {}, 'infer_done': {}, 'pilot_complete': False, 'smoke_done': False}
    if state['manifest_id'] != manifest['manifest_id']:
        raise ValueError('State belongs to another run')
    writer = StateWriter(gcs, path, manifest['output_uri']+'/state.json')
    n = manifest['cohort']['cohort']
    nshards = (n+cfg['shard_size']-1)//cfg['shard_size']

    def launch(stage, c=22, b=0):
        key = f'{stage}:{c}:{b}'
        count = state.setdefault('submission_retries', {}).get(key, 0)
        if count >= 8:
            raise ValueError('Repeated submission errors require review: '+key)
        # Keep the logical unit label stable across API-only retries so a late
        # accepted request can never become an unrelated, invisible unit.
        command_args, output, unit = command(manifest, stage, c, b)
        if stage == 'phase':
            uris = [cfg['source_uri']+f'/acaf_threshold.chr{c}.{ext}' for ext in ('pgen', 'pvar')]
            uris += [cfg['source_uri']+'/'+manifest['cohort']['shared_psam']]
            for uri in uris:
                if gcs.metadata(uri) != manifest['source_objects'][uri]:
                    raise ValueError('Source object changed since cohort preparation')
        job = {'job_id': None, 'unit': unit, 'stage': stage, 'chrom': c, 'shard': b,
               'output': output, 'state': 'SUBMITTING', 'submitted_at': time.time()}
        state['jobs'][key] = job
        writer.save(state)  # Reserve before the API call, including on notebook interruption.
        logpath = directory/(unit+'.submission.log')
        with logpath.open('w') as output_log:
            result = subprocess.run(command_args, stdout=output_log, stderr=subprocess.STDOUT, cwd=REPO)
        text = logpath.read_text()
        match = re.search(r'Launched job-id:\s*(\S+)', text)
        if result.returncode or not match:
            log('Submission interrupted; reconciling the recorded request before retry: '+unit)
            return False
        job.update(job_id=match[1], state='SUBMITTED')
        writer.save(state)
        log('Submitted '+unit+': '+match[1])
        return True

    from google.api_core import exceptions
    from requests.exceptions import ConnectionError, Timeout
    transient = (StatusLookupError, exceptions.ServiceUnavailable, exceptions.DeadlineExceeded,
                 exceptions.TooManyRequests, exceptions.InternalServerError, ConnectionError, Timeout)
    while True:
        try:
            cloud = statuses(cfg, run_label(manifest))
            units = Counter(row['labels'].get('lai-unit') for row in cloud.values())
            if any(unit and count > 1 for unit,count in units.items()):
                raise ValueError('Multiple cloud jobs share one logical unit; reconcile before continuing')
            for key, job in list(state['jobs'].items()):
                if job['state'] == 'DONE':
                    continue
                if job['job_id'] is None:
                    if reconcile_pending(manifest, directory, job, cloud):
                        state.setdefault('submission_history', []).append(dict(job))
                        count = state.get('submission_retries', {}).get(key, 0)
                        current = directory/(job['unit']+'.submission.log')
                        (directory/(job['unit']+f'.unaccepted_{count}.log')).write_bytes(current.read_bytes())
                        state.setdefault('submission_retries', {})[key] = state.get('submission_retries', {}).get(key, 0)+1
                        del state['jobs'][key]
                    writer.save(state)
                    continue
                observed = cloud.get(job['job_id'])
                if observed is None:
                    continue  # Missing API records never mean successful or unaccepted.
                job['state'] = observed['status']['state']
                if job['state'] == 'SUCCEEDED':
                    marker = gcs.completed(job['output'], manifest['manifest_id'])
                    require_marker(marker, manifest, job['stage'], job['chrom'], job['shard'])
                    compact = {k:v for k,v in marker.items() if k != 'files'}
                    if job['stage'] == 'smoke':
                        state['smoke_done'] = True
                    elif job['stage'] == 'phase':
                        state['phase_done'][str(job['chrom'])] = compact
                    else:
                        state['infer_done'][f"{job['chrom']}:{job['shard']}"] = compact
                    job['state'] = 'DONE'
                    writer.save(state)
                elif job['state'] in ('FAILED', 'DELETION_IN_PROGRESS'):
                    state.setdefault('failures', {})[key] = observed
                    writer.save(state)
                    # Continue unrelated chromosomes, but never automatically
                    # rerun paid failed compute or claim the dataset is complete.
            if '22' in state['phase_done'] and all(f'22:{b}' in state['infer_done'] for b in range(nshards)):
                state['pilot_complete'] = True
            if not any(job['job_id'] is None for job in state['jobs'].values()):
                requests = []
                if not state['smoke_done'] and 'smoke:22:0' not in state['jobs']:
                    requests.append(('smoke', 22, 0))
                for c in eligible_chromosomes(state):
                    if str(c) not in state['phase_done'] and f'phase:{c}:0' not in state['jobs']:
                        requests.append(('phase', c, 0))
                for c in sorted(map(int, state['phase_done'])):
                    for b in range(nshards):
                        if f'{c}:{b}' not in state['infer_done'] and f'infer:{c}:{b}' not in state['jobs']:
                            requests.append(('infer', c, b))
                for request in requests:
                    if not launch(*request):
                        break
            writer.save(state)
            log(f"{len(state['phase_done'])}/22 chromosomes phased; {len(state['infer_done'])}/{22*nshards} Gnofix shards verified; {len(state.get('failures', {}))} failed tasks")
            if len(state['phase_done']) == 22 and len(state['infer_done']) == 22*nshards:
                from .results import aggregate
                aggregate(manifest, directory, gcs)
                state['complete'] = True
                writer.save(state)
                log('All22 union results verified and published in the workspace bucket')
                return
        except transient as error:
            log('Temporary cloud error; preserving reservations and retrying: '+type(error).__name__)
        if args.once:
            return
        time.sleep(30)


def status(args):
    directory = args.work.resolve()/'runs'/args.run_id
    manifest = read_json(directory/'manifest.json')
    path = directory/'state.json'
    state = read_json(path) if path.exists() else {'jobs': {}}
    cloud = statuses(manifest['config'], run_label(manifest))
    rows = []
    for c in range(1, 23):
        job = state['jobs'].get(f'phase:{c}:0')
        record = cloud.get(job['job_id']) if job and job.get('job_id') else None
        phase = record['status']['state'] if record else job['state'] if job else 'NOT_SUBMITTED'
        inference = [j for j in state['jobs'].values() if j['stage']=='infer' and j['chrom']==c]
        done = sum(cloud.get(j['job_id'], {}).get('status', {}).get('state')=='SUCCEEDED' or j['state']=='DONE' for j in inference)
        rows.append({'chrom':c, 'phase':phase, 'ram_gib':manifest['resources']['chromosomes'][str(c)]['ram_gib'],
                     'inference_submitted':len(inference), 'inference_succeeded':done})
    print(json.dumps({'output_uri':manifest['output_uri'], 'pilot_complete':state.get('pilot_complete', False),
                      'complete':state.get('complete', False), 'chromosomes':rows}, indent=2))
