#!/usr/bin/env python3
"""wilson_sweep.py: Wilson lower-bound rules against plain occupancy cutoffs, on no-filter pileups.

    python wilson_sweep.py --mod m6A --bed 08_01=<nofilter>/C_d212_v52/08_01_25_RNA004_HEK293.bed --bed 09_09=... \\
        --bed 10_06=... --ortho-dir /projects/Genometechlab/andrew/orthogonal --fasta GRCh38.fa --mask-bed clusters.bed \\
        --out-dir wilson/v5.2_m6A

Every rule sees positions with >= 20 valid reads (masked loci and alt/decoy contigs removed). Families:
  pct_<c>            percent modified >= c, for c from 10 to 50 (pct_20 is the current rule)
  pct_20_min<d>      percent modified >= 20 and >= d valid reads (d = 30, 50): how much of a Wilson gain is
                     just dropping low-depth sites
  wilson<conf>_t<t>  the <conf>% (90, 95, 99) Wilson lower bound of the modified fraction >= t%
  wilsonbg<conf>_d<d>  (--background empirical) the same lower bound >= the site's background + d points
                     (d = 2.5 to 15). Background: cmp_rules.background on the run itself, i.e. the 5-mer's pooled
                     modified fraction after dropping the top (1 - --trim) of that context's sites by percent modified,
                     falling back to the 3-mer, then the modification's pooled rate, below --min-context-sites sites.
                     Modified sites that survive the trim raise their context's background, so true sites in
                     heavily modified contexts (DRACH for m6A) pay a penalty; background_<run>.tsv records the table.
  wilsonbg<conf>_t10_d<d>  lower bound >= max(10%, background + d points), d = 2.5, 5, 10, conf 95 and 90: wilson<conf>_t10
                     wherever the background is below 10% - d
  wilsonbg_givt07_<conf>_d<d>, wilsonbg_givt07_<conf>_t10_d<d>  (--background givt07:<table>) the same two families with
                     the background from a genomic-IVT 9-mer table (columns mer, mod, false_positive_0.7_threshold, ...,
                     support; 9-mers on the RNA strand, centred on the site): false_positive_0.7_threshold for the site's
                     9-mer and the modification's code, or the empirical 5-mer background where support < 200 or the
                     9-mer is missing. It is a lower bound on no-filter background: the table counts calls with
                     probability >= 0.70 as modified and everything else as canonical, over reference-matched 9-mers.
                     Valid only for the model the table was built from (the GM12878 table: v5.2, Dorado 1.0.0).
  fpcorr_pct20       (with givt07) the pipeline's cov20_fpcorr_pct20: percent modified - 100 * the gIVT 0.70 rate >= 20
                     (bin/fp_correct.py + mod_stages.sh fpcorr): any 9-mer in the table is corrected, whatever its
                     support; sites the table cannot correct are never called, as in the pipeline (no fallback)
--background can be given more than once. background_summary.tsv: per run and background, median background at DRACH
and other sites (5-mer on the RNA strand; DRACH is m6A's motif, the split is reported for every modification).
--fail-as-canonical (for pileups at a fixed threshold, e.g. 0.70): depth N_valid + N_fail and occupancy
N_mod / (N_valid + N_fail), so calls below the threshold count as canonical, as in the gIVT tables; every rule
(including --min-valid) then works on that depth and occupancy.
For each rule and run: calls, precision, recall against the modification's main truth set, or --truth-set (any set
cmp_ortho builds for the modification, e.g. m6A-SAC-seq, miCLIP-CITS) (precision is a lower bound: the maps are incomplete); recall split by the truth set's own level (does a rule gain real
low-stoichiometry sites, and how much of the map is low-level at all?); precision by depth band; agreement
between runs, separately at positions where both runs have < 100 or >= 100 valid reads (depth repeats between
runs, so pooled agreement favours rules that call by depth). The comparison that matters is at matched recall:
each rule's precision minus the plain-cutoff curve's precision at the same recall (interpolated per run).
Outputs in --out-dir: sweep.tsv, recall_by_level.tsv, by_depth.tsv, between_runs.tsv, matched.tsv, sweep.pdf/png
(and background_<run>.tsv, background_summary.tsv with --background).
Needs cmp_rules.py and cmp_ortho.py in the same directory.
"""
import argparse
import itertools
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmp_ortho as co  # noqa: E402
import cmp_rules as cr  # noqa: E402

PCT = [10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 35, 40, 50]
CONF = {90: 1.645, 95: 1.960, 99: 2.576}
TARGETS = [5, 7.5, 10, 12.5, 15, 17.5, 20]
DELTAS = [2.5, 5, 7.5, 10, 12.5, 15]
T10_CONF, T10_DELTAS = [95, 90], [2.5, 5, 10]
GIVT_MIN_SUPPORT = 200
LEVEL_BINS = [0, 10, 20, 40, 101]
LEVEL_LABELS = ['<10%', '10-20%', '20-40%', '>=40%']
DEPTH_BINS = [20, 50, 100, 300, np.inf]
DEPTH_LABELS = ['20-49', '50-99', '100-299', '300+']


def rules(d, p0s=None, givt_fp=None):
    """p0s: {rule prefix: per-site background rate (fraction)}, e.g. {'wilsonbg': empirical, 'wilsonbg_givt07_': gIVT};
    None leaves out the background families. givt_fp: per-site gIVT rate, NaN where the table cannot correct
    (adds fpcorr_pct20)."""
    pct, k, n = d['pct'].to_numpy(), d['nmod'].to_numpy(), d['nvalid'].to_numpy()
    out = {}
    for c in PCT:
        out[f'pct_{c:g}'] = pct >= c
    for dmin in (30, 50):
        out[f'pct_20_min{dmin}'] = (pct >= 20) & (n >= dmin)
    lbs = {conf: cr.wilson_lb(k, n, z) for conf, z in CONF.items()}
    for conf, lb in lbs.items():
        for t in TARGETS:
            out[f'wilson{conf}_t{t:g}'] = lb >= t / 100
    for pre, p0 in (p0s or {}).items():
        for conf, lb in lbs.items():
            for dl in DELTAS:
                out[f'{pre}{conf}_d{dl:g}'] = lb >= p0 + dl / 100
        for conf in T10_CONF:
            for dl in T10_DELTAS:
                out[f'{pre}{conf}_t10_d{dl:g}'] = lbs[conf] >= np.maximum(0.10, p0 + dl / 100)
    if givt_fp is not None:
        with np.errstate(invalid='ignore'):
            out['fpcorr_pct20'] = pct - 100 * givt_fp >= 20          # NaN rate: False, as mod_stages.sh drops it
    return out


def family(r):
    if r.startswith('wilsonbg'):
        return r.rsplit('_d', 1)[0]
    return 'pct' if r.startswith('pct_') and '_min' not in r else ('depth floor' if '_min' in r else r.split('_')[0])


def load_givt(path, code):
    """gIVT 9-mer table -> (false_positive_0.7_threshold, support), indexed by 9-mer, for one modification code."""
    t = pd.read_csv(path, sep='\t', usecols=['mer', 'mod', 'false_positive_0.7_threshold', 'support'], dtype={'mod': str})
    t = t[t['mod'] == code].set_index('mer')
    if not len(t):
        sys.exit(f'{path}: no rows for modification code {code}')
    return t['false_positive_0.7_threshold'], t['support']


def givt_p0(table, km9, fallback):
    """per-site gIVT background for (n,) 9-mer strings; the fallback where support < GIVT_MIN_SUPPORT or the 9-mer is
    missing. Returns the rates and the mask of sites that used the table."""
    fp, sup = table
    k = pd.Series(km9)
    use = (k.map(sup).fillna(0) >= GIVT_MIN_SUPPORT).to_numpy()
    return np.where(use, k.map(fp).to_numpy(dtype=float), fallback), use


def drach(k5):
    """(n, 5) uint8 5-mers on the RNA strand -> DRACH at the centre (D = A/G/U, R = A/G, H = A/C/U)."""
    return (np.isin(k5[:, 0], list(b'AGT')) & np.isin(k5[:, 1], list(b'AG')) & (k5[:, 2] == ord('A'))
            & (k5[:, 3] == ord('C')) & np.isin(k5[:, 4], list(b'ACT')))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--mod', default='m6A', choices=sorted(cr.TRUTH))
    ap.add_argument('--bed', action='append', required=True, metavar='NAME=PATH')
    ap.add_argument('--ortho-dir')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--truth-tsv', help='testing: chrom, pos, strand[, level]')
    ap.add_argument('--truth-set', help="cmp_ortho set to score against (default: the modification's main set)")
    ap.add_argument('--fail-as-canonical', action='store_true',
                    help='depth N_valid + N_fail, occupancy N_mod / (N_valid + N_fail): for fixed-threshold pileups')
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-valid', type=int, default=20)
    ap.add_argument('--background', action='append', default=[], metavar='empirical|givt07:TABLE',
                    help='add the wilsonbg rules for this background (repeat for both); see the description')
    ap.add_argument('--trim', type=float, default=0.9, help='for --background empirical, as in cmp_rules.py')
    ap.add_argument('--min-context-sites', type=int, default=30, help='for --background empirical, as in cmp_rules.py')
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    use_emp, givt = False, None
    for b in a.background:
        if b == 'empirical':
            use_emp = True
        elif b.startswith('givt07:'):
            givt = load_givt(b.partition(':')[2], co.CODE[a.mod])
        else:
            ap.error(f'--background {b}: expected empirical or givt07:<table>')
    os.makedirs(a.out_dir, exist_ok=True)
    runs = cr.pairs(a.bed)
    st = a.truth_set or cr.TRUTH[a.mod]
    if a.truth_tsv:
        truth = pd.read_csv(a.truth_tsv, sep='\t')
    else:
        listed, _ = co.load_sets(a.ortho_dir, '10ng_mRNA_input')
        listed = listed[listed['mod'] == a.mod]
        truth, _ = co.check_against_reference(listed, co.Fasta(a.fasta))
        truth = co.add_consensus(truth)
        sets = sorted(truth['set'].unique())
        truth = truth[truth['set'] == st]
        if not len(truth):
            sys.exit(f'no {a.mod} sites in truth set {st}; sets: {", ".join(sets)}')
    if 'level' not in truth:
        truth['level'] = np.nan
    mask = pd.read_csv(a.mask_bed, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                       comment='#') if a.mask_bed else None
    ids = {}
    for line in open(a.fasta + '.fai'):
        ids[line.split('\t')[0]] = len(ids)
    tk = co.keys(truth['chrom'], truth['pos'].to_numpy(), truth['strand'], ids)
    lvl = pd.Series(truth['level'].to_numpy(), index=tk)
    lvl = lvl[~lvl.index.duplicated()]
    fa = cr.FastaMM(a.fasta) if a.background else None

    sweep, by_level, by_depth, store, bg_summary = [], [], [], {}, []
    for name, path in runs:
        d = cr.read_pileup(path, co.CODE[a.mod], a.min_valid, a.contig_regex, mask, a.fail_as_canonical)
        keys = co.keys(d['chrom'], d['pos'].to_numpy(), d['strand'], ids)
        in_t = np.isin(keys, tk)
        tl = lvl.reindex(keys).to_numpy()
        lbin = pd.cut(tl, LEVEL_BINS, right=False, labels=LEVEL_LABELS)
        dbin = pd.cut(d['nvalid'].to_numpy(), DEPTH_BINS, right=False, labels=DEPTH_LABELS)
        ncov = int(in_t.sum())
        p0s, givt_fp = {}, None
        if fa is not None:
            chrom, pos, strand = d['chrom'].to_numpy(), d['pos'].to_numpy(), d['strand'].to_numpy()
            km = fa.kmers(chrom, pos, strand)
            p0, bgt = cr.background(km.view('S5').ravel().astype(str), d['nvalid'].to_numpy(), d['nmod'].to_numpy(),
                                    d['pct'].to_numpy(), a.trim, a.min_context_sites)
            bgt.to_csv(f'{a.out_dir}/background_{name}.tsv', sep='\t', index=False, float_format='%.5g')
            is_drach = drach(km)
            bgs = [('empirical', p0, np.ones(len(d), dtype=bool))]
            if use_emp:
                p0s['wilsonbg'] = p0
            if givt is not None:
                km9 = fa.kmers(chrom, pos, strand, k=9).view('S9').ravel().astype(str)
                pg, used = givt_p0(givt, km9, p0)
                p0s['wilsonbg_givt07_'] = pg
                givt_fp = pd.Series(km9).map(givt[0]).to_numpy(dtype=float)   # as fp_correct.py: any support
                bgs.append(('givt07', pg, used))
            for bname, p, used in bgs:
                bg_summary.append({'run': name, 'background': bname, 'sites': len(d), 'table_used': used.mean(),
                                   'median_drach': np.median(p[is_drach]) if is_drach.any() else np.nan,
                                   'median_other': np.median(p[~is_drach]), 'median_all': np.median(p),
                                   'max': p.max(), 'drach_sites': int(is_drach.sum())})
                print(f'{name}: {bname} background median {100 * np.median(p):.2f}%, max {100 * p.max():.2f}%; DRACH '
                      f'{100 * bg_summary[-1]["median_drach"]:.2f}%, other {100 * bg_summary[-1]["median_other"]:.2f}%'
                      + (f'; table used at {used.mean():.2%} of sites' if bname != 'empirical' else ''), file=sys.stderr)
        calls = rules(d, p0s, givt_fp)
        for r, c in calls.items():
            sweep.append({'run': name, 'rule': r, 'family': family(r), 'called': int(c.sum()),
                          'precision': in_t[c].mean() if c.any() else np.nan, 'recall': in_t[c].sum() / ncov if ncov else np.nan,
                          'truth_covered': ncov})
            for lb_ in LEVEL_LABELS:
                sel = in_t & (lbin == lb_)
                by_level.append({'run': name, 'rule': r, 'truth_level': lb_, 'truth_sites': int(sel.sum()),
                                 'recall': c[sel].mean() if sel.any() else np.nan})
            for db in DEPTH_LABELS:
                sel = c & (dbin == db)
                by_depth.append({'run': name, 'rule': r, 'depth_band': db, 'called': int(sel.sum()),
                                 'precision': in_t[sel].mean() if sel.any() else np.nan})
        store[name] = (keys, d['nvalid'].to_numpy(), calls)
        print(f'{name}: {len(d):,} positions, {ncov:,} {st} sites covered; truth level of covered sites: '
              + ', '.join(f'{lb_} {int(((lbin == lb_) & in_t).sum()):,}' for lb_ in LEVEL_LABELS), file=sys.stderr)

    between = []
    for x, y in itertools.combinations(store, 2):
        kx, nx, cx = store[x]; ky, ny, cy = store[y]
        _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
        deep = np.minimum(nx[ix], ny[iy]) >= 100
        for r in cx:
            A, B_ = cx[r][ix], cy[r][iy]
            for stratum, sel in (('both < 100 reads', ~deep), ('both >= 100 reads', deep)):
                a_, b_ = A & sel, B_ & sel
                between.append({'rule': r, 'run_a': x, 'run_b': y, 'stratum': stratum, 'called_a': int(a_.sum()),
                                'called_b': int(b_.sum()),
                                'repeat': np.nanmean([(a_ & b_).sum() / a_.sum() if a_.any() else np.nan,
                                                      (a_ & b_).sum() / b_.sum() if b_.any() else np.nan])})
    S, L, D, Bt = pd.DataFrame(sweep), pd.DataFrame(by_level), pd.DataFrame(by_depth), pd.DataFrame(between)

    # matched recall: each rule's precision minus the plain-cutoff curve's precision at the same recall, per run
    m = []
    for run, g in S.groupby('run'):
        curve = g[g['family'] == 'pct'].sort_values('recall')
        for _, row in g.iterrows():
            ref = np.interp(row['recall'], curve['recall'], curve['precision'], left=np.nan, right=np.nan)
            m.append({'run': run, 'rule': row['rule'], 'family': row['family'], 'called': row['called'],
                      'recall': row['recall'], 'precision': row['precision'], 'pct_curve_precision_at_recall': ref,
                      'precision_gain_at_matched_recall': row['precision'] - ref})
    M = pd.DataFrame(m)
    for df, nm in ((S, 'sweep'), (L, 'recall_by_level'), (D, 'by_depth'), (Bt, 'between_runs'), (M, 'matched')):
        df.to_csv(f'{a.out_dir}/{nm}.tsv', sep='\t', index=False, float_format='%.4g')
    if bg_summary:
        pd.DataFrame(bg_summary).to_csv(f'{a.out_dir}/background_summary.tsv', sep='\t', index=False, float_format='%.5g')

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), layout='constrained')
    ax = axes[0]
    mean = S.groupby(['rule', 'family'], as_index=False)[['recall', 'precision', 'called']].mean()
    style = {'pct': ('k', 'o', 'plain cutoff (10-50%)'), 'wilson90': ('#9ecae1', 's', 'Wilson 90%'),
             'wilson95': ('#3182bd', 's', 'Wilson 95%'), 'wilson99': ('#08519c', 's', 'Wilson 99%'),
             'depth floor': ('#e6550d', 'D', 'pct 20, >= 30 or 50 reads'),
             'wilsonbg90': ('#a1d99b', '^', 'Wilson 90% over 5-mer background'),
             'wilsonbg95': ('#31a354', '^', 'Wilson 95% over 5-mer background'),
             'wilsonbg99': ('#006d2c', '^', 'Wilson 99% over 5-mer background'),
             'wilsonbg95_t10': ('#31a354', 'v', 'Wilson 95% >= max(10%, 5-mer background + d)'),
             'wilsonbg_givt07_95': ('#756bb1', '^', 'Wilson 95% over gIVT 0.7 background'),
             'wilsonbg_givt07_95_t10': ('#756bb1', 'v', 'Wilson 95% >= max(10%, gIVT 0.7 background + d)'),
             'fpcorr': ('#e6550d', 'X', 'pipeline gIVT correction (pct - fp >= 20)')}
    single = ('depth floor', 'fpcorr')                                # one point or two, not curves
    for fam, (col, mk, lab) in style.items():
        g = mean[mean['family'] == fam].sort_values('recall')
        if g.empty:
            continue
        ax.plot(g['recall'], g['precision'], linestyle='none' if fam in single else '-', marker=mk, ms=3.5, color=col, lw=1, label=lab)
    for r, lab in (('pct_20', 'pct 20'), ('wilson95_t10', 'wilson10')):
        row = mean[mean['rule'] == r]
        if len(row):
            ax.annotate(lab, (row['recall'].iloc[0], row['precision'].iloc[0]), fontsize=7, xytext=(4, 4), textcoords='offset points')
    ax.set_xlabel('recall'); ax.set_ylabel('precision (lower bound)'); ax.set_title(f'{a.mod} against {st}, mean of runs', fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    ax = axes[1]
    mm = M.groupby(['rule', 'family'], as_index=False)[['recall', 'precision_gain_at_matched_recall']].mean()
    for fam, (col, mk, lab) in style.items():
        if fam == 'pct':
            continue
        g = mm[mm['family'] == fam].sort_values('recall')
        if g.empty:
            continue
        ax.plot(g['recall'], 100 * g['precision_gain_at_matched_recall'], linestyle='none' if fam in single else '-', marker=mk,
                ms=3.5, color=col, lw=1, label=lab)
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlabel('recall'); ax.set_ylabel('precision gain over plain cutoffs\nat the same recall (points)')
    ax.set_title('Better than a plain cutoff? (above 0 = yes)', fontsize=8); ax.legend(frameon=False, fontsize=6.5)
    fig.savefig(f'{a.out_dir}/sweep.pdf'); fig.savefig(f'{a.out_dir}/sweep.png', dpi=200); plt.close(fig)

    with pd.option_context('display.width', 200, 'display.max_rows', 100):
        keep = ['pct_15', 'pct_17.5', 'pct_20', 'pct_22.5', 'pct_25', 'pct_20_min30', 'pct_20_min50',
                'wilson90_t10', 'wilson95_t7.5', 'wilson95_t10', 'wilson95_t12.5', 'wilson95_t15', 'wilson99_t10',
                'wilsonbg95_d2.5', 'wilsonbg95_d5', 'wilsonbg95_d7.5', 'wilsonbg95_d10', 'wilsonbg95_t10_d5',
                'wilsonbg95_t10_d10', 'wilsonbg_givt07_95_d5', 'wilsonbg_givt07_95_d7.5', 'wilsonbg_givt07_95_t10_d5',
                'wilsonbg_givt07_90_d10', 'wilsonbg_givt07_95_d10', 'fpcorr_pct20']
        t = M.groupby('rule')[['called', 'recall', 'precision', 'precision_gain_at_matched_recall']].agg(['mean', 'min', 'max'])
        print(t.reindex([k for k in keep if k in t.index]).round(4).to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()