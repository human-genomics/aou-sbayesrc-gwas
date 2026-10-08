"""Small deterministic I/O helpers; participant identifiers never enter public logs."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.part')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def derived_seed(seed, *keys):
    # Stable across Python versions, processes, shards, and retry order.
    return int.from_bytes(hashlib.sha256(json.dumps([int(seed), *map(str, keys)], separators=(',', ':')).encode()).digest()[:4], 'big')


def log(message):
    print(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), message, flush=True)


def run(args, **kwargs):
    log('exec: ' + ' '.join(map(str, args)))
    return subprocess.run(list(map(str, args)), check=True, **kwargs)


def read_psam(path):
    import pandas as pd
    df = pd.read_csv(path, sep=r'\s+', dtype=str, keep_default_na=False)
    df.columns = df.columns.str.lstrip('#')
    if 'FID' not in df:
        df['FID'] = '0'
    if 'IID' not in df or df.IID.duplicated().any() or (df.IID == '').any():
        raise ValueError('PSAM must have unique nonempty IIDs')
    return df


def complete(directory, manifest_id, details):
    directory = Path(directory)
    files = {str(p.relative_to(directory)): {'sha256': sha256(p), 'bytes': p.stat().st_size}
             for p in sorted(directory.rglob('*')) if p.is_file() and p.name != 'COMPLETE.json'}
    write_json(directory / 'COMPLETE.json', {'manifest_id': manifest_id, 'files': files, **details})


def verify_complete(directory, manifest_id, checksums=True):
    directory = Path(directory)
    path = directory / 'COMPLETE.json'
    if not path.exists():
        return False
    d = read_json(path)
    if d['manifest_id'] != manifest_id:
        raise ValueError(f'Run fingerprint mismatch: {directory}')
    for name, record in d['files'].items():
        p = directory / name
        if not p.is_file() or p.stat().st_size != record['bytes']:
            return False
        if checksums and sha256(p) != record['sha256']:
            return False
    return True


def configure_runtime(runtime):
    import sys
    runtime = Path(runtime)
    for name in ('gnomix', 'adapter'):
        sys.path.insert(0, str(runtime / name))
    os.environ['GNOMIX1000G_TOOLS'] = str(runtime)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = '1'
    os.environ['MPLBACKEND'] = 'Agg'
    os.environ['MPLCONFIGDIR'] = '/tmp/gnomix-matplotlib'
