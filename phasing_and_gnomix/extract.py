"""Select normalized model SNVs and apply only cohort-wide missingness QC."""
from __future__ import annotations

import gzip
import os
from pathlib import Path
import numpy as np
import pandas as pd
from phasing_and_gnomix.normalize_pvar_alleles import right_trim, normalize
from phasing_and_gnomix.common import read_psam, run, write_json
from phasing_and_gnomix.models import window_table


def pvar_records(path):
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t', 5)
            yield fields[:5]


def pairkey(pos, ref, alt):
    return f'{int(pos)}:{min(ref, alt)}:{max(ref, alt)}'


def nearest_windows(meta, positions):
    """Nearest reliably lifted model site; upstream site wins a distance tie."""
    valid = np.flatnonzero(meta['pos38'] > 0)
    order = np.argsort(meta['pos38'][valid], kind='stable')
    idx = valid[order]
    pos = meta['pos38'][idx]
    right = np.minimum(np.searchsorted(pos, positions), len(pos)-1)
    left = np.maximum(right-1, 0)
    chosen = np.where(np.abs(positions-pos[left]) <= np.abs(pos[right]-positions), left, right)
    model_idx = idx[chosen]
    return model_idx, np.minimum(model_idx // int(meta['M']), int(meta['W'])-1)



def candidate_ids(pvar, meta, output, sb=None):
    positions = set(meta['pos38'][meta['pos38'] > 0].tolist())
    if sb is not None:positions.update(sb.pos.astype(int).tolist())
    count = 0
    seen = set()
    with open(output, 'w') as f:
        for chrom, pos, vid, ref, alt in pvar_records(pvar):
            if int(pos) in positions and any(len(r) == len(a) == 1 and r in 'ACGT' and a in 'ACGT'
                                             for r,a in (right_trim(ref, x) for x in alt.split(','))):
                if vid == '.' or vid in seen:
                    raise ValueError('Candidate source record IDs must be unique and nonmissing')
                seen.add(vid)
                f.write(vid + '\n')
                count += 1
    if not count:
        raise ValueError('No candidate source records')
    return count


def panel_frame(path):
    rows = [(i, int(pos), vid, ref, alt) for i,(_,pos,vid,ref,alt) in enumerate(pvar_records(path))
            if ref in ('A','C','G','T') and alt in ('A','C','G','T') and ref != alt]
    return pd.DataFrame(rows, columns=['pidx','pos','id','pref','palt'])


def match(meta, panel):
    model = pd.DataFrame({'i': np.arange(len(meta['pos38'])), 'pos': meta['pos38'],
                          'mref': meta['ref38'], 'malt': meta['alt38']})
    joined = model[model.pos > 0].merge(panel, on='pos')
    same = (joined.mref == joined.pref) & (joined.malt == joined.palt)
    swapped = (joined.mref == joined.palt) & (joined.malt == joined.pref)
    exact = joined.loc[same | swapped].copy()
    exact['flip'] = swapped.loc[exact.index].values
    ambiguous = exact.i.duplicated(keep=False) | exact.pidx.duplicated(keep=False)
    status = np.full(len(model), 'absent', dtype='U24')
    status[model.pos <= 0] = 'failed_liftover'
    status[np.unique(joined.i)] = 'allele_mismatch'
    status[np.unique(exact.loc[ambiguous, 'i'])] = 'ambiguous'
    result = exact.loc[~ambiguous, ['i','pidx','id','flip']].sort_values('i')
    status[result.i.values] = 'matched'
    return result, status


def passing_missingness(missing, observed, ceiling):
    missing, observed = np.asarray(missing), np.asarray(observed)
    return (observed > 0) & (missing <= ceiling * observed)


def prepare(chrom, source, psam, keep, meta_path, output, scratch, plink, config):
    source, output, scratch = map(Path, (source, output, scratch))
    plink_args = [plink, '--seed', str(config['seeds']['plink'])]
    output.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    meta = dict(np.load(meta_path, allow_pickle=False))
    sb=pd.read_csv(os.environ['SBAYES'])
    if not (sb.chrom==chrom).all():raise ValueError('Wrong SBayesRC chromosome')
    sb['key']=[pairkey(x.pos,x.ref,x.alt) for x in sb.itertuples(index=False)]
    if sb.key.duplicated().any() or not sb.ref.isin(list('ACGT')).all() or not sb.alt.isin(list('ACGT')).all():raise ValueError('SBayesRC must contain unique SNVs')
    extract = scratch / 'candidate.ids'
    n_candidates = candidate_ids(source.with_suffix('.pvar'), meta, extract, sb)
    candidates = scratch / 'candidates'
    # PLINK applies --set-all-var-ids before --extract. Subset original IDs in
    # a separate pass before assigning split-variant IDs.
    run([*plink_args, '--pgen', source.with_suffix('.pgen'), '--pvar', source.with_suffix('.pvar'), '--psam', psam,
         '--keep', keep, '--extract', extract, '--make-pgen', 'erase-phase',
         '--no-pheno', '--output-chr', 'chrM', '--threads', config['beagle_threads'], '--out', candidates])
    split = scratch / 'split'
    run([*plink_args, '--pfile', candidates, '--make-pgen', 'multiallelics=-', 'erase-phase',
         '--set-all-var-ids', '@:#:$r:$a', '--new-id-max-allele-len', '10000',
         '--no-pheno', '--output-chr', 'chrM', '--threads', config['beagle_threads'], '--out', split])
    # All transformations happen on worker-local files, never mounted source data.
    normalized = scratch / 'normalized.pvar'
    trim_stats = normalize(split.with_suffix('.pvar'), normalized)
    normalized.replace(split.with_suffix('.pvar'))
    panel = panel_frame(split.with_suffix('.pvar'))
    matched, status = match(meta, panel)
    if matched.empty:
        raise ValueError('No unambiguous model matches')
    ids = scratch / 'model.ids'
    panel['key']=[pairkey(x.pos,x.pref,x.palt) for x in panel.itertuples(index=False)]
    duplicate_keys=set(panel.loc[panel.key.duplicated(keep=False),'key'])
    extra=panel[panel.key.isin(set(sb.key)) & ~panel.key.duplicated(keep=False) & ~panel.id.isin(matched.id)].copy()
    sb['status']='source_absent_or_other_alleles'
    sb.loc[sb.key.isin(duplicate_keys),'status']='ambiguous_source'
    sb.loc[sb.key.isin(panel.loc[panel.id.isin(matched.id),'key']) | sb.key.isin(extra.key),'status']='missingness_excluded'
    pd.concat([matched.id,extra.id]).to_csv(ids, header=False, index=False)
    metrics = scratch / 'missingness'
    run([*plink_args, '--pfile', split, '--extract', ids, '--missing', 'variant-only', '--out', metrics])
    miss = pd.read_csv(metrics.with_suffix('.vmiss'), sep=r'\s+')
    if miss.ID.duplicated().any():
        raise ValueError('Ambiguous normalized IDs survived matching')
    miss['pass'] = passing_missingness(miss.MISSING_CT, miss.OBS_CT, config['max_variant_missingness'])
    passing = set(miss.loc[miss['pass'], 'ID'])
    status[matched.loc[~matched.id.isin(passing), 'i'].values] = 'missingness_excluded'
    matched = matched.loc[matched.id.isin(passing)].copy()
    extra=extra[extra.id.isin(passing)].copy()
    sb.loc[sb.key.isin(extra.key),'status']='extra_passed'
    sb.loc[sb.key.isin(panel.loc[panel.id.isin(matched.id),'key']),'status']='model_passed'
    sb.to_csv(output/'sbayesrc_feature_qc.tsv.gz',sep='\t',index=False)
    coverage = len(matched) / int(meta['C'])
    # Write the exclusion audit even when coverage fails.
    feature_qc = pd.DataFrame({'model_index': np.arange(len(status)), 'pos_hg19': meta['pos19'],
                              'pos_hg38': meta['pos38'], 'ref_model': meta['ref'], 'alt_model': meta['alt'], 'status': status})
    feature_qc.to_csv(output / 'feature_qc.tsv.gz', sep='\t', index=False)
    miss.to_csv(output / 'variant_missingness.tsv.gz', sep='\t', index=False)
    if coverage < config['minimum_model_coverage']:
        raise ValueError(f'chr{chrom}: retained model coverage {coverage:.4f} below required floor')
    pd.concat([matched.id,extra.id]).to_csv(ids,header=False,index=False)
    qc = scratch / 'qc'
    run([*plink_args, '--pfile', split, '--extract', ids, '--make-pgen', 'erase-phase', '--sort-vars',
         '--no-pheno', '--output-chr', 'chrM', '--out', qc])
    out_panel = panel_frame(qc.with_suffix('.pvar'))
    final, _ = match(meta, out_panel)
    if len(final) != len(matched) or len(out_panel) != len(final)+len(extra):
        raise ValueError('Final extraction changed the model match set')
    np.savez_compressed(output / 'match.npz', model_idx=final.i.values, pgen_idx=final.pidx.values,
                        flip=final.flip.values.astype(bool), pos38=meta['pos38'],pgen_variants=np.asarray(len(out_panel)))
    out_panel['key']=[pairkey(x.pos,x.pref,x.palt) for x in out_panel.itertuples(index=False)]
    out_panel['model_index'],out_panel['window']=nearest_windows(meta,out_panel.pos.to_numpy())
    out_panel.loc[final.pidx.values,'model_index']=final.i.values
    out_panel.loc[final.pidx.values,'window']=np.minimum(final.i.values//int(meta['M']),int(meta['W'])-1)
    scoring=out_panel.merge(sb[['key','ref','alt','rsid']],on='key',validate='one_to_one')
    scoring['chrom']=chrom;scoring['sbayesrc_ref_alt_swapped']=scoring.pref!=scoring.ref
    scoring=scoring.rename(columns={'pidx':'source_variant_index','pref':'ref_source','palt':'alt_source'}).sort_values(['pos','id'],kind='stable').reset_index(drop=True)
    scoring['scoring_variant_index']=np.arange(len(scoring));scoring['source']='union'
    scoring.to_csv(output/'scoring_variants.tsv.gz',sep='\t',index=False)
    win = window_table(meta, final.i.values, chrom)
    win.to_csv(output / 'windows.tsv', sep='\t', index=False)
    summary = {'chrom': chrom, 'candidate_records': n_candidates, 'normalization': trim_stats,
               'model_sites': int(meta['C']), 'retained_sites': len(final), 'coverage': coverage,
               'samples': len(read_psam(qc.with_suffix('.psam'))),
               'feature_status': pd.Series(status).value_counts().astype(int).to_dict(),
               'missingness_threshold': config['max_variant_missingness']}
    summary.update(union_sites=len(out_panel),extra_sites=len(extra),variant_selection='Gnomix union SBayesRC',sbayesrc_status_counts=sb.status.value_counts().astype(int).to_dict())
    write_json(output / 'extraction.json', summary)
    return qc
