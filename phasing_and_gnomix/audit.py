"""Independent consistency checks of tracts, posteriors and Gnofix swaps."""
import io
import numpy as np
import pandas as pd

ANCESTRIES = ['EUR', 'EAS', 'NAT', 'AFR', 'SAS', 'AHG', 'OCE', 'WAS']
TRACT_COLUMNS = ['sample', 'haplotype', 'chrom', 'start_hg38', 'end_hg38',
                 'start_hg19', 'end_hg19', 'start_cM', 'end_cM', 'ancestry',
                 'n_windows', 'mean_posterior']


def require(condition, message):
    if not condition:
        raise ValueError('Artifact audit: ' + message)


def table(raw, compressed=True):
    return pd.read_csv(io.BytesIO(raw), sep='\t', compression='gzip' if compressed else None,
                       dtype={'sample': str, 'IID': str, 'FID': str}, keep_default_na=False)


def arrays(raw):
    with np.load(io.BytesIO(raw), allow_pickle=False) as data:
        return {name: data[name] for name in data.files}


def audit_unit_data(tracts, status, switches, saved_totals, parts, windows, samples, chrom):
    """Reconstruct labels/totals from tracts and compare independent saved files.

    parts yields (sample_start, sample_end, checkpoint_dict) in sample order.
    No Gnomix/1000G inference or tract-writing function is used here.
    """
    n, W = len(samples), len(windows)
    require(n > 0 and W > 0, 'empty sample or window set')
    require(list(tracts.columns) == TRACT_COLUMNS, 'tract schema differs')
    require(not tracts.empty, 'empty tract table')
    require(len(set(samples)) == n, 'duplicate expected samples')
    require(np.array_equal(windows.window, np.arange(W)), 'window index differs')
    require((windows.chrom == chrom).all(), 'window chromosome differs')
    require(not windows.spos_hg19.duplicated().any() and not windows.epos_hg19.duplicated().any(),
            'ambiguous model-window endpoints')
    require((tracts.chrom == chrom).all(), 'tract chromosome differs')
    ids = {iid: i for i, iid in enumerate(samples)}
    si = tracts['sample'].map(ids)
    hi = tracts.haplotype.map({'A': 0, 'B': 1})
    ai = tracts.ancestry.map({a: i for i, a in enumerate(ANCESTRIES)})
    require(si.notna().all() and hi.notna().all() and ai.notna().all(),
            'unknown sample, haplotype, or ancestry')
    si, hi, ai = (x.to_numpy(dtype=np.int64) for x in (si, hi, ai))
    hap = 2*si + hi
    require((np.diff(hap) >= 0).all(), 'tract sample/haplotype order differs')
    require(np.array_equal(np.unique(hap), np.arange(2*n)), 'missing sample haplotype')
    start = tracts.start_hg19.map(dict(zip(windows.spos_hg19, range(W))))
    end = tracts.end_hg19.map(dict(zip(windows.epos_hg19, range(W))))
    require(start.notna().all() and end.notna().all(), 'tract endpoint is not a model-window endpoint')
    start, end = start.to_numpy(np.int64), end.to_numpy(np.int64)
    counts = tracts.n_windows.to_numpy()
    require(np.isfinite(counts).all() and np.array_equal(counts, end-start+1) and (counts > 0).all(),
            'tract window count differs')
    counts = counts.astype(np.int64)
    first = np.r_[True, np.diff(hap) != 0]
    last = np.r_[np.diff(hap) != 0, True]
    require((start[first] == 0).all() and (end[last] == W-1).all(), 'incomplete chromosome endpoints')
    continuation = np.flatnonzero(~first)
    require((start[continuation] == end[continuation-1]+1).all(), 'overlapping or missing tract windows')
    require((ai[continuation] != ai[continuation-1]).all(), 'adjacent equal-ancestry tracts were not merged')
    require(np.allclose(tracts.start_cM, windows.sgpos.to_numpy()[start], rtol=0, atol=1e-9)
            and np.allclose(tracts.end_cM, windows.egpos.to_numpy()[end], rtol=0, atol=1e-9),
            'tract genetic coordinates differ')
    posterior = tracts.mean_posterior.to_numpy(dtype=float)
    require(np.isfinite(posterior).all() and ((posterior >= 0) & (posterior <= 1)).all(),
            'invalid tract posterior')

    # Use interval searches rather than the writer's prefix/suffix span algorithm.
    require(windows.hg38_discordant.isin([True, False]).all(), 'invalid liftover flags')
    good = np.flatnonzero(~windows.hg38_discordant.to_numpy(dtype=bool))
    left = np.searchsorted(good, start)
    right = np.searchsorted(good, end, side='right')-1
    usable = (left <= right) & (left < len(good)) & (right >= 0)
    s38, e38 = np.full(len(tracts), -1, np.int64), np.full(len(tracts), -1, np.int64)
    s38[usable] = windows.spos_hg38.to_numpy()[good[left[usable]]]
    e38[usable] = windows.epos_hg38.to_numpy()[good[right[usable]]]
    inverted = e38 < s38
    s38[inverted] = -1
    e38[inverted] = -1
    require(np.array_equal(tracts.start_hg38, s38) and np.array_equal(tracts.end_hg38, e38),
            'tract GRCh38 coordinates differ')

    weights = np.diff(np.r_[windows.sgpos.to_numpy(), windows.egpos.iloc[-1]])
    require(np.isfinite(weights).all() and (weights >= 0).all() and weights.sum() > 0,
            'invalid genetic window lengths')
    cumulative = np.r_[0., np.cumsum(weights)]
    totals = np.zeros((n, len(ANCESTRIES)), np.float64)
    np.add.at(totals, (si, ai), cumulative[end+1]-cumulative[start])
    require(list(saved_totals['samples']) == list(samples)
            and list(saved_totals['populations']) == ANCESTRIES, 'saved totals identities differ')
    require(saved_totals['totals'].shape == totals.shape
            and np.allclose(saved_totals['totals'], totals, rtol=0, atol=1e-7),
            'saved ancestry totals disagree with tracts')

    require(list(status['sample']) == list(samples) and (status.chrom == chrom).all(),
            'sample completion identities differ')
    require(status.gnofix_completed.eq(True).all() and (status.n_windows_per_haplotype == W).all(),
            'a sample lacks complete Gnofix processing')
    labels = np.repeat(ai.astype(np.uint8), counts).reshape(2*n, W)
    next_sample = 0
    expected_sw_samples, expected_sw_windows, agreements = [], [], []
    for lo, high, part in parts:
        require(lo == next_sample and lo < high <= n, 'checkpoint sample coverage differs')
        next_sample = high
        require(list(part['samples']) == list(samples[lo:high]), 'checkpoint sample identities differ')
        require(np.array_equal(part['labels'], labels[2*lo:2*high]), 'tract labels disagree with Gnofix checkpoint')
        require(part['totals'].shape == totals[lo:high].shape
                and np.allclose(part['totals'], totals[lo:high], rtol=0, atol=1e-7),
                'checkpoint ancestry totals disagree with tracts')
        support = part['called_posterior']
        require(support.shape == (2*(high-lo), W) and np.isfinite(support).all()
                and ((support >= 0) & (support <= 1)).all(), 'invalid checkpoint posteriors')
        sums = np.pad(np.cumsum(support, axis=1, dtype=np.float64), ((0, 0), (1, 0)))
        selected = (si >= lo) & (si < high)
        rows = hap[selected]-2*lo
        means = (sums[rows, end[selected]+1]-sums[rows, start[selected]])/counts[selected]
        # Tracts round the original probabilities to 4 decimals; checkpoints
        # store float32. The pinned upstream writer also accumulates float32
        # model probabilities in float32 before casting its prefixes to float64.
        # Cancellation in short, late tracts can exceed decimal-rounding error.
        # Check that exact numeric path as well, rather than widening tolerance
        # around the independent float64 mean for unrelated discrepancies.
        sums32 = np.pad(np.cumsum(support, axis=1, dtype=np.float32).astype(np.float64),
                        ((0, 0), (1, 0)))
        means32 = (sums32[rows, end[selected]+1]-sums32[rows, start[selected]])/counts[selected]
        matches64 = np.abs(means-posterior[selected]) <= 0.0000501
        matches32 = np.abs(np.round(means32, 4)-posterior[selected]) <= 1e-12
        require(np.all(matches64 | matches32), 'tract posterior differs from checkpoint')
        swaps = part['swaps']
        require(swaps.shape == (high-lo, W) and np.isin(swaps, [0, 1]).all(), 'invalid Gnofix swap matrix')
        ss, ww = np.nonzero(np.diff(swaps.astype(np.int8), axis=1, prepend=0) != 0)
        expected_sw_samples.extend(np.asarray(samples)[lo+ss].tolist())
        expected_sw_windows.extend(ww.tolist())
        agree = part['raw_fix_diploid_agreement']
        require(agree.shape == (high-lo,) and np.isfinite(agree).all()
                and ((agree >= 0) & (agree <= 1)).all(), 'invalid Gnofix agreement measure')
        agreements.extend(agree.tolist())
    require(next_sample == n, 'missing Gnofix checkpoints')
    require(list(switches['sample']) == expected_sw_samples
            and list(switches.window) == expected_sw_windows
            and (switches.chrom == chrom).all(), 'stored switches disagree with Gnofix tracker')
    ww = np.asarray(expected_sw_windows, dtype=np.int64)
    require(np.array_equal(switches.pos_hg19, windows.spos_hg19.to_numpy()[ww])
            and np.array_equal(switches.pos_hg38, windows.spos_hg38.to_numpy()[ww])
            and np.allclose(switches.cM.to_numpy(dtype=float), windows.sgpos.to_numpy()[ww], rtol=0, atol=1e-9),
            'switch coordinates differ')
    switch_counts = pd.Series(expected_sw_samples, dtype=str).value_counts()
    require(np.array_equal(status.n_switches, switch_counts.reindex(samples, fill_value=0)), 'sample switch counts differ')
    require(np.allclose(status.raw_fix_diploid_agreement, agreements, rtol=0, atol=1e-12),
            'sample agreement differs from checkpoint')
    return totals, {'tract_rows': len(tracts), 'window_calls': 2*n*W, 'switches': len(ww)}
