"""Audit original checkpoints and prepare their fixed SNV/window mapping."""
from __future__ import annotations

import argparse
import gc
import pickle
from pathlib import Path
import shutil
import tarfile
import numpy as np
import pandas as pd

from .common import configure_runtime, read_json, write_json, log, sha256


def load_model(path):
    # Imported after configure_runtime; never use an arbitrary unpinned pickle.
    from gnomix1000g.prepare import _quiet_xgboost_teardown
    _quiet_xgboost_teardown()
    with open(path, 'rb') as f:
        m = pickle.load(f)
    m.base.base_multithread = False
    m.base.log_inference = False
    m.base.n_jobs = 1
    m.smooth.n_jobs = 1
    m.smooth.model.n_jobs = 1
    m.smooth.model.get_booster().set_param('nthread', 1)
    return m


def assert_snv_alleles(ref, alt):
    # Check original strings, never a truncating U1 conversion.
    ref, alt = np.asarray(ref).astype(str), np.asarray(alt).astype(str)
    if ref.ndim != 1 or alt.shape != ref.shape:
        raise ValueError('Unexpected model allele dimensions')
    if not (np.isin(ref, list('ACGT')) & np.isin(alt, list('ACGT')) & (ref != alt)).all():
        raise ValueError('Original checkpoint includes a non-SNV or invalid allele')
    return ref, alt


def prepare(work, config):
    work = Path(work)
    configure_runtime(work / 'runtime')
    from gnomix1000g.liftover import Chain, complement
    chain = Chain(work / 'downloads/hg19ToHg38.over.chain.gz')
    dest = work / 'references/models'
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(work / 'downloads/pretrained_gnomix_models.tar.gz', 'r:gz') as tar:
        wanted = {f'pretrained_gnomix_models/chr{c}/model_chm_{c}.pkl': c for c in config['chromosomes']}
        for member in tar:
            if member.name in wanted:
                c = wanted[member.name]
                path = dest / f'chr{c}.pkl'
                # This path can be reached after a cached checkpoint failed
                # its checksum. Restore it from the verified public archive.
                with tar.extractfile(member) as src, open(path.with_suffix('.part'), 'wb') as out:
                    shutil.copyfileobj(src, out, 8 << 20)
                path.with_suffix('.part').replace(path)
    rows = []
    for c in config['chromosomes']:
        log(f'Auditing original chr{c} checkpoint')
        m = load_model(dest / f'chr{c}.pkl')
        ref, alt = assert_snv_alleles(m.snp_ref, m.snp_alt)
        if list(m.population_order) != config['ancestries'] or len(ref) != m.C or m.W != m.C // m.M:
            raise ValueError('Unexpected checkpoint dimensions or ancestry order')
        pos19 = np.asarray(m.snp_pos, dtype=np.int64)
        ch, pos38, minus = chain.lift(f'chr{c}', pos19)
        valid = ch == f'chr{c}'
        pos38[~valid] = -1
        np.savez_compressed(dest / f'chr{c}.meta.npz', pos19=pos19, pos38=pos38,
                            ref=ref, alt=alt, ref38=np.where(minus, complement(ref), ref),
                            alt38=np.where(minus, complement(alt), alt), minus=minus,
                            C=m.C, M=m.M, W=m.W, populations=np.asarray(m.population_order),
                            gm_pos=np.asarray(m.gen_map_df.pos, dtype=float), gm_cm=np.asarray(m.gen_map_df.pos_cm, dtype=float))
        rows.append({'chrom': c, 'sites': int(m.C), 'non_snv': 0, 'lifted_same_chrom': int(valid.sum()),
                     'shared_position_sites': int(pd.Series(pos19).duplicated(keep=False).sum()),
                     'checkpoint_sha256': sha256(dest / f'chr{c}.pkl')})
        del m
        gc.collect()
    if sum(r['sites'] for r in rows) != config['expected_model_sites']:
        raise ValueError('Model archive has an unexpected total feature count')
    write_json(dest / 'audit.json', rows)


def window_table(meta, model_idx, chrom):
    from gnomix1000g.tracts import hg38_discordant
    C, M, W = (int(meta[k]) for k in ('C', 'M', 'W'))
    starts = np.arange(W) * M
    ends = np.append(starts[1:], C) - 1
    w = np.minimum(np.arange(C) // M, W-1)
    df = pd.DataFrame({'chrom': chrom, 'window': np.arange(W), 'spos_hg19': meta['pos19'][starts], 'epos_hg19': meta['pos19'][ends]})
    df['sgpos'] = np.round(np.interp(df.spos_hg19, meta['gm_pos'], meta['gm_cm']), 5)
    df['egpos'] = np.round(np.interp(df.epos_hg19, meta['gm_pos'], meta['gm_cm']), 5)
    df['n_model_snps'] = np.bincount(w, minlength=W)
    df['n_query_snps'] = np.bincount(w[model_idx], minlength=W)
    lifted = np.flatnonzero(meta['pos38'] > 0)
    first = np.full(W, -1, dtype=np.int64)
    last = first.copy()
    if len(lifted):
        lw = w[lifted]
        left = np.searchsorted(lw, np.arange(W), side='left')
        right = np.searchsorted(lw, np.arange(W), side='right') - 1
        has = right >= left
        first[has], last[has] = meta['pos38'][lifted[left[has]]], meta['pos38'][lifted[right[has]]]
    df['spos_hg38'], df['epos_hg38'] = np.minimum(first,last), np.maximum(first,last)
    df['hg38_discordant'] = hg38_discordant(first,last,(df.epos_hg19-df.spos_hg19).values)
    return df


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    a=p.parse_args()
    prepare(a.work.resolve(), read_json(a.config))
