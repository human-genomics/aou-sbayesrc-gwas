"""Freeze the exact complement-plus-random-Europeans cohort."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from .common import read_psam, sha256, write_json


def keep_ids(path):
    d = pd.read_csv(path, sep=r'\s+', header=None, dtype=str, keep_default_na=False)
    if d.shape[1] != 2 or d[1].duplicated().any():
        raise ValueError('Expected a unique two-column FID/IID keep list')
    return set(d[1])


def select(psam, european_ids, fit_ids, size, seed):
    wgs = set(psam.IID)
    if not fit_ids <= european_ids <= wgs:
        raise ValueError('Required set relationship: fit_pca <= classified_european <= WGS')
    if len(fit_ids) < size:
        raise ValueError('Insufficient PCA-fit Europeans')
    candidates = np.array(sorted(fit_ids), dtype=str)
    rng = np.random.Generator(np.random.PCG64(seed))
    selected = set(rng.choice(candidates, size=size, replace=False))
    non_eur = wgs - european_ids
    result = psam.loc[psam.IID.isin(non_eur | selected), ['FID', 'IID']].copy()
    result['source_index'] = np.flatnonzero(psam.IID.isin(non_eur | selected))
    result['reason'] = np.where(result.IID.isin(selected), 'random_fit_pca', 'not_classified_european')
    result['sample_index'] = np.arange(len(result))
    assert len(result) == len(non_eur) + size
    return result


def build(config, source_mount, pipeline_mount, dest):
    source_mount, pipeline_mount, dest = map(Path, (source_mount, pipeline_mount, dest))
    dest.mkdir(parents=True, exist_ok=True)
    paths = {'europeans': pipeline_mount / 'europeans/classified_european_iids.txt',
             'fit_pca': pipeline_mount / 'pca_eur/fit_pca_iids.txt',
             'admixture': pipeline_mount / 'statgen/aou_admixture_k6.tsv'}
    psams = [source_mount / f'acaf_threshold.chr{c}.psam' for c in config['chromosomes']]
    published = [p for p in psams if p.is_file()]
    if not published or len({sha256(p) for p in published}) != 1:
        raise ValueError('Published autosome PSAM files must be nonempty and identical')
    paths['psam'] = published[0]
    eur, fit = keep_ids(paths['europeans']), keep_ids(paths['fit_pca'])
    a = pd.read_csv(paths['admixture'], sep='\t', dtype={'IID': str, 'FID': str})
    psam = read_psam(paths['psam'])
    if a.IID.duplicated().any() or set(a.IID) != set(psam.IID):
        raise ValueError('ADMIXTURE sample universe differs from WGS')
    thresholds = config['european_classifier']
    classified = set(a.loc[(a.European >= thresholds['European_min'])
                           & (a.African <= thresholds['African_max'])
                           & (a.American <= thresholds['American_max'])
                           & (a.East_Asian <= thresholds['East_Asian_max'])
                           & (a.Oceanian <= thresholds['Oceanian_max']), 'IID'])
    if classified != eur:
        raise ValueError('Saved European list does not match the K=6 classifier')
    result = select(psam, eur, fit, config['european_sample_size'], config['seeds']['cohort'])
    result.to_csv(dest / 'samples.tsv', sep='\t', index=False)
    result[['FID', 'IID']].to_csv(dest / 'cohort.keep', sep='\t', header=False, index=False)
    result.loc[result.reason == 'random_fit_pca', ['FID', 'IID']].to_csv(dest / 'random_fit_pca.keep', sep='\t', header=False, index=False)
    summary = {'wgs': len(psam), 'classified_european': len(eur), 'fit_pca': len(fit),
               'not_classified_european': len(psam)-len(eur), 'selected_european': config['european_sample_size'],
               'cohort': len(result), 'seed': config['seeds']['cohort'],
               'inputs': {k: {'path': str(p), 'sha256': sha256(p)} for k,p in paths.items()},
               'shared_psam': published[0].name}
    write_json(dest / 'cohort.json', summary)
    return summary


if __name__ == '__main__':
    import argparse
    from .common import read_json
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('config','source','pipeline','dest'):p.add_argument('--'+name,required=True)
    a=p.parse_args()
    build(read_json(a.config),a.source,a.pipeline,a.dest)
