"""Resolve notebook settings without embedding a researcher's workspace in Git."""
import os
from pathlib import Path
import re
import subprocess

from .common import read_json

PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parent
DEFAULT_WORK = REPO / 'data/phasing_and_gnomix'
SOURCE_NAMES = (
    '.gitignore', '__init__.py', 'README.md', 'run.sh', 'run.py', 'settings.py',
    'prepare.py', 'scheduler.py', 'results.py', 'common.py', 'cloud.py', 'state.py',
    'bootstrap.py', 'models.py', 'cohort.py', 'extract.py', 'phase.py', 'infer.py',
    'worker.py', 'worker.sh', 'dsub_launcher.py', 'normalize_pvar_alleles.py',
    'validation.py', 'smoke.py', 'cloud_smoke.sh', 'audit.py',
    'config.json', 'resources.json', 'seeds.json', 'requirements.txt',
    'requirements-orchestrator.txt',
)


def source_files():
    # Explicit source-only allowlist: runtime, logs, cohorts and results cannot
    # enter a worker bundle even if they are mistakenly put beside the code.
    files = [PACKAGE / name for name in SOURCE_NAMES]
    files.extend(sorted((PACKAGE / 'tests').glob('*.py')))
    if not all(path.is_file() for path in files):
        raise ValueError('Incomplete phasing_and_gnomix source checkout')
    return files


def gcloud_value(key):
    value = subprocess.check_output(['gcloud', 'config', 'get-value', key],
                                    text=True, stderr=subprocess.DEVNULL).strip()
    if not value or value == '(unset)':
        raise ValueError('gcloud ' + key + ' is unset; use an AoU Jupyter terminal')
    return value


def mounted_workspace(mount_text=None):
    target = os.environ.get('WORKSPACE_BUCKET_MOUNT', '/home/jupyter/workspace/workspace-bucket')
    text = Path('/proc/mounts').read_text() if mount_text is None else mount_text
    for line in text.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[1].replace('\\040', ' ') == target and 'gcsfuse' in fields[2]:
            return 'gs://' + fields[0].removeprefix('gs://')
    return None


def require_gcs(uri):
    if not re.fullmatch(r'gs://[a-z0-9][a-z0-9._-]+(?:/[A-Za-z0-9_./-]+)?', uri) or '/..' in uri:
        raise ValueError('Expected a gs:// URI with an ordinary bucket/object path')
    return uri.rstrip('/')


def resolve(args):
    if not re.fullmatch(r'[a-z][a-z0-9-]{0,47}', args.run_id):
        raise ValueError('Run ID must start with a letter and use at most48 lowercase letters/digits/hyphens')
    cfg = read_json(args.config)
    cfg['seeds'] = read_json(args.seeds)
    resources = read_json(args.resources)
    if cfg['chromosomes'] != list(range(1, 23)):
        raise ValueError('This workflow requires all22 autosomes')
    if cfg['model_fill'] != 0 or cfg['pilot_chromosome'] != 22:
        raise ValueError('Reference filling and the chr22 pilot are required')
    if cfg['ancestries'] != ['EUR', 'EAS', 'NAT', 'AFR', 'SAS', 'AHG', 'OCE', 'WAS']:
        raise ValueError('The pinned model requires its original eight ancestry labels')
    if not 0 <= cfg['max_variant_missingness'] <= 1 or not 0 < cfg['minimum_model_coverage'] <= 1:
        raise ValueError('Invalid missingness/coverage threshold')
    if set(cfg['seeds']) != {'cohort', 'plink', 'beagle', 'inference', 'validation'}:
        raise ValueError('Missing reproducibility seed')
    if any(type(seed) is not int or not 0 <= seed < 2**31 for seed in cfg['seeds'].values()):
        raise ValueError('Seeds must be integers in [0, 2^31)')
    if any(type(cfg[k]) is not int or cfg[k] <= 0 for k in ('european_sample_size', 'shard_size', 'batch_size', 'infer_cores')):
        raise ValueError('Sample/shard/process settings must be positive integers')
    thresholds = cfg['european_classifier']
    if set(thresholds) != {'European_min', 'African_max', 'American_max', 'East_Asian_max', 'Oceanian_max'} or any(
            not isinstance(value, (int, float)) or not 0 <= value <= 1 for value in thresholds.values()):
        raise ValueError('European classifier thresholds must match the documented rule and lie in [0, 1]')
    if set(resources['chromosomes']) != set(map(str, range(1, 23))):
        raise ValueError('Specify resources for every chromosome')
    for row in resources['chromosomes'].values():
        if any(type(row[key]) is not int or row[key] <= 0 for key in ('ram_gib', 'heap_gb', 'threads', 'vcpus')):
            raise ValueError('Chromosome resources must be positive integers')
        if not 0 < row['heap_gb'] < row['ram_gib'] or not 0 < row['threads'] <= row['vcpus']:
            raise ValueError('Java heap must fit VM RAM and threads must fit CPUs')
    if any(type(resources[key]) is not int or resources[key] <= 0 for key in ('phase_disk_gb', 'phase_timeout_hours')):
        raise ValueError('Phasing disk size and timeout must be positive integers')
    if any(type(resources['inference'][key]) is not int or resources['inference'][key] <= 0
           for key in ('ram_gib', 'vcpus', 'disk_gb', 'timeout_hours')):
        raise ValueError('Inference resources must be positive integers')
    if cfg['infer_cores'] > resources['inference']['vcpus']:
        raise ValueError('Inference process count exceeds its VM CPU count')
    project = args.project or os.environ.get('GOOGLE_PROJECT') or gcloud_value('project')
    bucket = args.workspace_bucket or mounted_workspace() or os.environ.get('WORKSPACE_BUCKET_URI') or os.environ.get('WORKSPACE_BUCKET')
    if args.pipeline_uri:
        pipeline = require_gcs(args.pipeline_uri)
    elif bucket:
        pipeline = require_gcs(bucket) + '/' + os.environ.get('SBAYESRC_OUTPUT_PREFIX', 'sbayesrc_genotypes').strip('/')
    else:
        raise ValueError('Cannot discover the workspace bucket; pass --workspace-bucket gs://YOUR_BUCKET')
    version = args.data_version or os.environ.get('AOU_DATA_VERSION', 'v9')
    if version not in ('v8', 'v9'):
        raise ValueError('Use v8 or v9 matching get_genotypes.sh')
    source = args.source_uri or os.environ.get('AOU_PGEN_GS_DIR') or f'gs://vwb-aou-datasets-controlled/{version}/wgs/short_read/snpindel/acaf_threshold/pgen'
    cfg.update(project=project, pipeline_uri=pipeline, source_uri=require_gcs(source), data_version=version,
               region=args.region, service_account=args.service_account or gcloud_value('account'),
               network=args.network or f'projects/{project}/global/networks/network',
               subnetwork=args.subnetwork or f'projects/{project}/regions/{args.region}/subnetworks/subnetwork')
    output = require_gcs(args.output_uri or (pipeline + '/phasing_and_gnomix/' + args.run_id))
    if output == pipeline or output.startswith(cfg['source_uri']):
        raise ValueError('Output must be a new run directory in the workspace bucket')
    if args.work.resolve().is_relative_to(PACKAGE):
        raise ValueError('Place private working data outside the public source directory; the default is data/phasing_and_gnomix')
    return cfg, resources, output
