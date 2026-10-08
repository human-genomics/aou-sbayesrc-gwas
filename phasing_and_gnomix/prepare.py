"""Prepare public tools, freeze the cohort and stage immutable run inputs."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import gzip
import os
from pathlib import Path
import re
import subprocess
import shutil
import tarfile

from .cloud import GCS
from .common import fingerprint, log, read_json, run, sha256, write_json
from .settings import PACKAGE, REPO, resolve, source_files


def runtime_environment(work):
    rt = Path(work) / 'runtime'
    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = ':'.join(str(rt / x / 'lib') for x in ('libgomp', 'libstdcxx', 'libgcc')) + ':' + env.get('LD_LIBRARY_PATH', '')
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    return env


def references(work, config_path):
    """Reuse only a cache whose source/config and every public input still match."""
    work = Path(work)
    cfg = read_json(config_path)
    key = fingerprint({'config': cfg, 'source': {str(p.relative_to(REPO)): sha256(p) for p in source_files()}})
    ready = work / 'references/READY.json'
    if ready.exists():
        record = read_json(ready)
        if record['key'] == key and all((work / name).is_file() and sha256(work / name) == digest
                                      for name, digest in record['files'].items()):
            require_public_checks(work)
            return
    # Bootstrap runs in the notebook interpreter; models use their own pinned
    # Python3.9 runtime, independently of packages in the notebook environment.
    import sys
    run([sys.executable, '-m', 'phasing_and_gnomix.bootstrap', '--work', work, '--config', config_path])
    python = work / 'runtime/python/bin/python3'
    env = runtime_environment(work)
    run([python, '-m', 'phasing_and_gnomix.models', '--work', work, '--config', config_path], env=env)
    run([python, '-m', 'phasing_and_gnomix.validation', '--runtime', work/'runtime',
         '--model', work/'references/models/chr22.pkl', '--config', config_path,
         '--output', work/'references/demo_validation.json'], env=env)
    run([python, '-m', 'phasing_and_gnomix.smoke', '--runtime', work/'runtime',
         '--model', work/'references/models/chr22.pkl', '--meta', work/'references/models/chr22.meta.npz',
         '--config', config_path, '--output', work/'references/smoke_validation.json'], env=env)
    require_public_checks(work)
    inputs = read_json(work / 'references/public_inputs.json')
    inputs.update({str(p.relative_to(work)): sha256(p) for p in (work/'references/models').iterdir() if p.is_file()})
    for name in ('demo_validation.json', 'smoke_validation.json'):
        inputs['references/' + name] = sha256(work/'references'/name)
    write_json(ready, {'key': key, 'files': inputs})


def require_public_checks(work):
    """Never stage paid-work inputs from an incomplete validation report."""
    required = {
        'demo_validation.json': (
            'gnofix_reference_labels_equal', 'gnofix_reference_swaps_equal', 'batch_invariant'),
        'smoke_validation.json': (
            'all_samples_gnofix', 'batch_labels_identical', 'extra_union_sites_preserve_model_predictions',
            'independent_tract_audit_passed', 'observed_genotypes_preserved',
            'public_end_to_end_inference_passed', 'real_beagle_phase_passed',
            'real_union_selection_passed', 'scoring_phase_verified', 'shard_phase_verified',
            'tract_schema_valid'),
    }
    for name, keys in required.items():
        report = read_json(Path(work)/'references'/name)
        if any(report.get(key) is not True for key in keys):
            raise ValueError('Incomplete or failed public validation: ' + name)


def download(gcs, uri, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = gcs.blob(uri)
    blob.reload()
    record = {'generation': blob.generation, 'size': blob.size, 'crc32c': blob.crc32c}
    tmp = path.with_name(path.name + '.part')
    blob.download_to_filename(str(tmp), if_generation_match=blob.generation)
    tmp.replace(path)
    return record


def freeze_metadata(gcs, cfg, directory):
    source = directory / 'source_metadata'
    pipeline = directory / 'pipeline_metadata'
    published = [b for b in gcs.objects(cfg['source_uri'] + '/acaf_threshold.chr')
                 if re.search(r'/acaf_threshold\.chr(?:[1-9]|1[0-9]|2[0-2])\.psam$', b.name)]
    if not published or len({(b.size, b.crc32c) for b in published}) != 1 or not published[0].size:
        raise ValueError('Published autosome PSAMs must be nonempty and byte-identical')
    published.sort(key=lambda b: int(re.search(r'chr(\d+)\.psam$', b.name)[1]))
    first = published[0]
    objects = {}
    bucket = cfg['source_uri'].removeprefix('gs://').split('/', 1)[0]
    for b in published:
        objects['gs://' + bucket + '/' + b.name] = {'generation': b.generation, 'size': b.size, 'crc32c': b.crc32c}
    uri = 'gs://' + bucket + '/' + first.name
    if download(gcs, uri, source / Path(first.name).name) != objects[uri]:
        raise ValueError('PSAM changed during preparation')
    for name in ('europeans/classified_european_iids.txt', 'pca_eur/fit_pca_iids.txt', 'statgen/aou_admixture_k6.tsv'):
        uri = cfg['pipeline_uri'] + '/' + name
        objects[uri] = download(gcs, uri, pipeline / name)
    uris = [cfg['source_uri'] + f'/acaf_threshold.chr{c}.{ext}'
            for c in cfg['chromosomes'] for ext in ('pgen', 'pvar')]
    with ThreadPoolExecutor(8) as pool:
        objects.update(zip(uris, pool.map(gcs.metadata, uris)))
    return source, pipeline, objects


def source_archive(path, files):
    with Path(path).open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode='w') as archive:
            for file in sorted(files):
                info = archive.gettarinfo(str(file), arcname=str(file.relative_to(REPO)))
                info.uid = info.gid = 0
                info.uname = info.gname = ''
                info.mtime = 0
                with file.open('rb') as stream:
                    archive.addfile(info, stream)


def freeze_input(path, work):
    """A later bootstrap/preparation cannot mutate another run's local inputs."""
    path = Path(path)
    digest = sha256(path)
    target = Path(work)/'artifacts'/digest
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if sha256(target) != digest:
            raise ValueError('Corrupted content-addressed input cache')
    else:
        temporary = target.with_suffix('.part')
        shutil.copyfile(path, temporary)
        if sha256(temporary) != digest:
            raise ValueError('Input changed while being frozen')
        temporary.replace(target)
    return {'local_path':str(target), 'sha256':digest, 'bytes':target.stat().st_size}


def verify_frozen(manifest, directory, gcs):
    if fingerprint({k: v for k, v in manifest.items() if k != 'manifest_id'}) != manifest['manifest_id']:
        raise ValueError('Run manifest fingerprint differs')
    for name, digest in manifest['source_code'].items():
        if sha256(REPO / name) != digest:
            raise ValueError('Code changed since preparation; restore it or use a new run ID: ' + name)
    if sha256(directory / 'cohort/samples.tsv') != manifest['cohort_sha256']:
        raise ValueError('Frozen sample list changed')
    for name, row in manifest['input_files'].items():
        if sha256(Path(row['local_path'])) != row['sha256']:
            raise ValueError('Frozen staged input changed: ' + name)
    existing = gcs.json(manifest['output_uri'] + '/inputs/manifest.json')
    if existing is not None and existing != manifest:
        raise ValueError('GCS run prefix already contains a different manifest')


def stage(manifest, directory, gcs):
    verify_frozen(manifest, directory, gcs)
    root = manifest['output_uri'] + '/inputs/'
    for name, row in manifest['input_files'].items():
        log('Staging ' + name)
        gcs.put(Path(row['local_path']), root + name)
    gcs.put(directory / 'manifest.json', root + 'manifest.json')
    log('Prepared immutable inputs: ' + manifest['output_uri'])


def prepare(args):
    cfg, resources, output = resolve(args)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    # Preparation shares a public-reference cache across runs. Serialize it so
    # two terminals cannot replace a cache or a partial upload simultaneously.
    with (work/'prepare.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return prepare_locked(args, cfg, resources, output)


def prepare_locked(args, cfg, resources, output):
    work = args.work.resolve()
    directory = work / 'runs' / args.run_id
    directory.mkdir(parents=True, exist_ok=True)
    gcs = GCS(cfg['project'])
    if (directory / 'manifest.json').exists():
        manifest = read_json(directory / 'manifest.json')
        if cfg != manifest['config'] or resources != manifest['resources'] or output != manifest['output_uri']:
            raise ValueError('Prepared run settings differ; use a new run ID')
        stage(manifest, directory, gcs)
        return
    # Refuse a collision before downloading references or touching source lists.
    if gcs.json(output + '/inputs/manifest.json') is not None:
        raise ValueError('This GCS run already exists; resume with its original local run directory')
    write_json(directory / 'config.json', cfg)
    references(work, directory / 'config.json')
    source, pipeline, source_objects = freeze_metadata(gcs, cfg, directory)
    run([work/'runtime/python/bin/python3', '-m', 'phasing_and_gnomix.cohort',
         '--config', directory/'config.json', '--source', source, '--pipeline', pipeline,
         '--dest', directory/'cohort'], env=runtime_environment(work))
    cohort = read_json(directory / 'cohort/cohort.json')
    import pandas as pd
    sb = pd.read_csv(work / 'downloads/sbayesrc_hg38.csv')
    if list(sb.columns) != ['chrom', 'pos', 'ref', 'alt', 'rsid']:
        raise ValueError('Unexpected SBayesRC hg38 release schema')
    scoring_dir = directory / 'sbayesrc'
    scoring_dir.mkdir(exist_ok=True)
    policy = {'variant_set': 'gnomix_union_sbayesrc_hg38', 'chromosomes': {}}
    for c in cfg['chromosomes']:
        path = scoring_dir / f'chr{c}.csv'
        part = sb[sb.chrom == c]
        if part.empty:
            raise ValueError('Empty SBayesRC chromosome')
        part.to_csv(path, index=False)
        policy['chromosomes'][str(c)] = {**resources['chromosomes'][str(c)], 'sbayesrc_csv_sha256': sha256(path)}
    policy['policy_id'] = fingerprint(policy)
    files = source_files()
    source_archive(directory / 'code.tar.gz', files)
    # Runtime archive is content-addressed in this run's immutable manifest.
    archive = work / 'runtime.tar.gz'
    if not archive.exists():
        run(['tar', '-czf', archive, '-C', work/'runtime', '.'])
    inputs = {'runtime.tar.gz': archive, 'code.tar.gz': directory/'code.tar.gz', 'config.json': directory/'config.json'}
    for folder, prefix in ((directory/'cohort', 'cohort'), (scoring_dir, 'sbayesrc'),
                           (work/'references/models', 'models')):
        inputs.update({prefix + '/' + p.name: p for p in folder.iterdir() if p.is_file()})
    inputs.update({'maps/' + p.name: p for p in (work/'references/maps').glob('*.map')})
    manifest = {'run_id': args.run_id, 'output_uri': output, 'config': cfg, 'resources': resources,
                'cohort': cohort, 'cohort_sha256': sha256(directory/'cohort/samples.tsv'),
                'source_objects': source_objects, 'union_policy': policy,
                'source_code': {str(p.relative_to(REPO)): sha256(p) for p in files},
                'public_inputs': read_json(work/'references/public_inputs.json'),
                'model_audit': read_json(work/'references/models/audit.json'),
                'public_demo_validation': read_json(work/'references/demo_validation.json'),
                'public_smoke': read_json(work/'references/smoke_validation.json'),
                'input_files': {name: freeze_input(path, work) for name, path in inputs.items()}}
    manifest['manifest_id'] = fingerprint(manifest)
    write_json(directory / 'manifest.json', manifest)
    stage(manifest, directory, gcs)
