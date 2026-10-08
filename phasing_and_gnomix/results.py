"""Aggregate only complete union results with validated Gnofix scoring phase."""
import io
from pathlib import Path
import numpy as np
import pandas as pd

from .common import complete, read_json, verify_complete, write_json
from .scheduler import require_marker


def aggregate(manifest, directory, gcs):
    root = manifest['output_uri']
    uri = root+'/results'
    existing = gcs.completed(uri, manifest['manifest_id'])
    if existing is not None:
        if not existing.get('scoring_complete') or existing.get('chromosomes') != 22:
            raise ValueError('Invalid final completion marker')
        return
    out = Path(directory)/'results'
    if (out/'COMPLETE.json').exists() and verify_complete(out, manifest['manifest_id']):
        # Retry identical compressed bytes after an interrupted GCS upload.
        for path in sorted(out.iterdir(), key=lambda p:p.name=='COMPLETE.json'):
            gcs.put(path, uri+'/'+path.name)
        return
    cfg = manifest['config']
    n, size = manifest['cohort']['cohort'], cfg['shard_size']
    samples = pd.read_csv(Path(directory)/'cohort/samples.tsv', sep='\t', dtype={'FID':str, 'IID':str})
    if len(samples) != n:
        raise ValueError('Frozen cohort size changed')
    totals = np.zeros((n,len(cfg['ancestries'])), np.float64)
    seen = np.zeros(n, np.uint8)
    scores, tracts, chromosomes = [], [], []
    for c in range(1,23):
        phase = root+f'/phased/chr{c}'
        require_marker(gcs.completed(phase, manifest['manifest_id']), manifest, 'phase', c)
        chromosomes.append({'chrom':c, 'variant_set':'gnomix_union_sbayesrc_hg38',
            'beagle_pgen_uri':phase+f'/chr{c}.pgen', 'beagle_pvar_uri':phase+f'/chr{c}.pvar',
            'beagle_psam_uri':phase+f'/chr{c}.psam', 'windows_uri':phase+'/windows.tsv',
            'scoring_variant_map_uri':phase+'/scoring_variants.tsv.gz'})
        for b,left in enumerate(range(0,n,size)):
            right = min(left+size,n)
            unit = root+f'/inference/chr{c}/batch_{b:05d}'
            require_marker(gcs.completed(unit, manifest['manifest_id']), manifest, 'infer', c, b)
            with np.load(io.BytesIO(gcs.blob(unit+'/ancestry_totals.npz').download_as_bytes()), allow_pickle=False) as z:
                if list(z['samples']) != list(samples.IID.iloc[left:right]) or list(z['populations']) != cfg['ancestries']:
                    raise ValueError('Ancestry identities/order differ')
                if z['totals'].shape != (right-left,len(cfg['ancestries'])) or not np.isfinite(z['totals']).all() or (z['totals']<0).any():
                    raise ValueError('Invalid ancestry totals')
                totals[left:right] += z['totals']
                seen[left:right] += 1
            shared = {'chrom':c, 'shard':b, 'sample_start_index':left, 'sample_end_index_exclusive':right}
            tracts.append({**shared, 'tracts_uri':unit+'/tracts_gnofix.tsv.gz'})
            scores.append({**shared, 'pgen_uri':unit+'/scoring_genotypes_gnofix.pgen',
                'pvar_uri':unit+'/scoring_genotypes_gnofix.pvar', 'psam_uri':unit+'/scoring_genotypes_gnofix.psam',
                'scoring_variant_map_uri':phase+'/scoring_variants.tsv.gz', 'calls_uri':unit+'/parts/',
                'scoring_manifest_uri':unit+'/scoring_manifest.json'})
    if not (seen==22).all() or not (totals.sum(axis=1)>0).all():
        raise ValueError('Every sample must have complete calls on all22 chromosomes')
    out.mkdir(parents=True, exist_ok=True)
    fractions = totals/totals.sum(axis=1)[:,None]
    result = samples[['FID','IID','reason']].copy()
    for k,pop in enumerate(cfg['ancestries']):
        result[pop], result[pop+'_pct'] = fractions[:,k], 100*fractions[:,k]
    result['total_cM'], result['completed_chromosomes'] = totals.sum(axis=1)/2, seen
    tables = {'global_ancestry_gnofix':result, 'scoring_file_index':pd.DataFrame(scores),
              'tract_file_index':pd.DataFrame(tracts), 'chromosome_index':pd.DataFrame(chromosomes),
              'sample_index':samples.assign(shard=samples.sample_index//size)}
    for name,frame in tables.items():
        frame.to_csv(out/(name+'.tsv.gz'), sep='\t', index=False, float_format='%.10g',
                     compression={'method':'gzip', 'mtime':0})
    write_json(out/'DATASET_VERSION.json', {'manifest_id':manifest['manifest_id'],
        'variant_set':'gnomix_union_sbayesrc_hg38', 'union_policy_id':manifest['union_policy']['policy_id'],
        'ancestries':cfg['ancestries'], 'scoring_phase':'Beagle plus Gnofix swaps',
        'posterior':'Probability of the called ancestry per model window',
        'additional_sbayesrc_sites':'Nearest reliably lifted model site selects the same window for swap, ancestry and posterior'})
    complete(out, manifest['manifest_id'], {'chromosomes':22, 'samples':n, 'scoring_complete':True,
             'gnofix_completed_sample_chromosomes':int(seen.sum()), 'union_policy_id':manifest['union_policy']['policy_id']})
    for path in sorted(out.iterdir(), key=lambda p:p.name=='COMPLETE.json'):
        gcs.put(path, uri+'/'+path.name)
