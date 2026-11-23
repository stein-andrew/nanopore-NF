#!/usr/bin/env python3
"""cmp_plots.py: compare modification calls between basecalling setups run on the SAME reads.

    python cmp_plots.py --sample GM12878 \\
        --run A_d212_v6=pileup064/A_d212_v6/GM12878.bed \\
        --run B_d130_v52=pileup064/B_d130_v52/GM12878.bed \\
        --run C_d212_v52=pileup064/C_d212_v52/GM12878.bed \\
        --out plots/GM12878 [--mask-bed clusters.bed] [--bootstrap 200] [--style my.mplstyle]

Input: modkit pileup bedMethyl (18 columns), one per run; a <prefix>.thresholds.tsv next to each
(written by cmp_pileup.sh) is picked up for the threshold panel. Occupancy is modkit's raw percent
(no gIVT correction: the 9-mer table is v5.2.0's). Rows on unplaced/alt contigs (--contig-regex) and
inside --mask-bed intervals (e.g. the rRNA-like clusters) are dropped for every run alike.

Definitions (defaults --min-cov 20 --min-pct 20):
  covered     N_valid_cov >= min-cov in that run
  valid site  covered and percent >= min-pct
  comparable  covered in every run being compared, so a site can only differ by its call, not by
              its coverage. Venns are drawn on the comparable set (three runs) and unrestricted.
  Jaccard     |valid X & valid Y| / |valid X | valid Y| within the comparable set
  ICC(2,1)    two-way random, absolute agreement, single measure (Shrout & Fleiss 1979), on
              comparable sites that are valid in at least one run (all-zero sites would dominate)
  r, RMSE     over all comparable sites, every occupancy (what the heatmaps show)
  95% CIs     --bootstrap N: resample 100 kb genomic blocks with replacement (sites in one gene are
              not independent; blocks stand in for genes without needing an annotation)

Outputs in --out: figures/{venn_comparable,venn_unrestricted,heatmaps,summary}.pdf, one heatmap PDF
per modification, tables/*.tsv, legends.md.
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle

MODS = [('a', 'm6A'), ('17596', 'Inosine'), ('69426', '2OmeA'), ('m', 'm5C'),
        ('19228', '2OmeC'), ('17802', 'Psi'), ('19227', '2OmeU'), ('19229', '2OmeG')]
MOD_COLORS = {'m6A': '#1f77b4', 'm5C': '#2ca02c', 'Inosine': '#ff7f0e', 'Psi': '#d62728',
              '2OmeA': '#9467bd', '2OmeC': '#8c564b', '2OmeG': '#e377c2', '2OmeU': '#7f7f7f'}
VENN_PATCH_COLORS = ['skyblue', 'lightgreen', 'salmon', 'purple', 'gold', 'orange', 'gray']
RUN_COLORS = ['skyblue', 'lightgreen', 'salmon', '#e377c2', '#bcbd22']
PAIR_COLORS = ['#0072B2', '#009E73', '#D55E00']
KNOWN_LABELS = {'A_d212_v6': 'v6.0.0 + v2 mods\n(dorado 2.1.2)', 'B_d130_v52': 'v5.2.0\n(dorado 1.3.0)',
                'C_d212_v52': 'v5.2.0\n(dorado 2.1.2)'}
SHORT = {'A_d212_v6': 'v6 / 2.1.2', 'B_d130_v52': 'v5.2 / 1.3.0', 'C_d212_v52': 'v5.2 / 2.1.2'}
# the three comparisons that each change one thing, when runs A, B, C are present
PREFERRED_PAIRS = [('A_d212_v6', 'C_d212_v52', 'models (same dorado)'),
                   ('C_d212_v52', 'B_d130_v52', 'dorado version (same models)'),
                   ('A_d212_v6', 'B_d130_v52', 'new vs production')]
BLOCK = 100_000


# ---------------------------------------------------------------- style
def set_style(style_path):
    if style_path:
        plt.style.use(style_path)
        return f'style sheet {style_path}'
    try:
        from genometechlab_plotting import setup_style
        setup_style('inline')
        return 'genometechlab_plotting'
    except Exception:
        plt.rcParams.update({'axes.grid': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                             'axes.spines.top': True, 'axes.spines.right': True})
        return 'matplotlib defaults (no grid)'


# ---------------------------------------------------------------- loading
def read_bedmethyl(path, min_cov, contig_re, mask):
    cols = [0, 1, 3, 5, 9, 10, 16]
    names = ['chrom', 'start', 'code', 'strand', 'cov', 'pct', 'ndiff']
    parts = []
    reader = pd.read_csv(path, sep='\t', header=None, usecols=cols, names=names, chunksize=5_000_000,
                         dtype={'chrom': str, 'start': np.int64, 'code': str, 'strand': str,
                                'cov': np.int64, 'pct': np.float64, 'ndiff': np.int64})
    n_all = 0
    for ch in reader:
        n_all += len(ch)
        ch = ch[ch['cov'] >= min_cov]
        if contig_re is not None:               # test each contig name once, not every row
            bad = [c for c in ch['chrom'].unique() if contig_re.search(c)]
            ch = ch[~ch['chrom'].isin(bad)]
        parts.append(ch)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=names)
    n_masked = 0
    for c, s, e in mask:
        hit = (df['chrom'] == c) & (df['start'] >= s) & (df['start'] < e)
        n_masked += int(hit.sum())
        df = df[~hit]
    return df, n_all, n_masked


def read_mask(path):
    if not path:
        return []
    out = []
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(('#', 'track', 'browser')):
                f = line.split('\t')
                out.append((f[0], int(f[1]), int(f[2])))
    return out


def read_thresholds(bed):
    p = re.sub(r'\.bed$', '', bed) + '.thresholds.tsv'
    if not os.path.exists(p):
        return None
    t = pd.read_csv(p, sep='\t', header=None, names=['base', 'threshold'])
    return dict(zip(t['base'], t['threshold']))


# ---------------------------------------------------------------- statistics
def icc_stats(x):
    """Sufficient statistics of an n x k matrix for ICC(2,1): n, sum x^2, sum of squared row sums,
    column sums (k)."""
    return np.concatenate([[x.shape[0], (x ** 2).sum(), (x.sum(axis=1) ** 2).sum()], x.sum(axis=0)])


def icc_from_stats(s, k):
    n, sxx, srr = s[..., 0], s[..., 1], s[..., 2]
    csum = s[..., 3:3 + k]
    T = csum.sum(axis=-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        cf = T ** 2 / (n * k)
        ssr = srr / k - cf
        ssc = (csum ** 2).sum(axis=-1) / n - cf
        sse = sxx - cf - ssr - ssc
        msr, msc, mse = ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))
        icc = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    return np.where(n > 2, icc, np.nan)


def block_bootstrap(block_ids, per_site, reps, rng):
    """per_site: n x m array of additive per-site quantities. Returns (reps x m) totals over
    blocks resampled with replacement."""
    codes, inv = np.unique(block_ids, return_inverse=True)
    sums = np.zeros((len(codes), per_site.shape[1]))
    np.add.at(sums, inv, per_site)
    w = rng.multinomial(len(codes), np.full(len(codes), 1 / len(codes)), size=reps)
    return w @ sums


def ci(arr):
    arr = arr[np.isfinite(arr)]
    return (np.nan, np.nan) if len(arr) < 10 else tuple(np.percentile(arr, [2.5, 97.5]))


def pair_metrics(sub, x, y, min_pct, reps, rng):
    both = sub[f'cov_{x}'].notna() & sub[f'cov_{y}'].notna()
    d = sub[both]
    px, py = d[f'pct_{x}'].to_numpy(), d[f'pct_{y}'].to_numpy()
    vx, vy = px >= min_pct, py >= min_pct
    inter, union = (vx & vy), (vx | vy)
    out = {'n_comparable': len(d), 'n_valid_x': int(vx.sum()), 'n_valid_y': int(vy.sum()),
           'n_both_valid': int(inter.sum()), 'n_either_valid': int(union.sum())}
    out['jaccard'] = inter.sum() / union.sum() if union.sum() else np.nan
    m = union
    xm = np.column_stack([px[m], py[m]])
    out['n_icc'] = int(m.sum())
    out['icc'] = float(icc_from_stats(icc_stats(xm), 2)) if m.sum() > 2 else np.nan
    if len(d) > 2:
        out['r_all'] = float(np.corrcoef(px, py)[0, 1]) if px.std() > 0 and py.std() > 0 else np.nan
        out['rmse_all'] = float(np.sqrt(np.mean((px - py) ** 2)))
        out['mean_diff_all'] = float(np.mean(py - px))
    else:
        out['r_all'] = out['rmse_all'] = out['mean_diff_all'] = np.nan
    for k in ('jaccard', 'icc'):
        out[f'{k}_lo'] = out[f'{k}_hi'] = np.nan
    if reps and len(d) > 2:
        blocks = (d['chrom'].astype(str) + ':' + (d['start'] // BLOCK).astype(str)).to_numpy()
        jac = block_bootstrap(blocks, np.column_stack([inter, union]).astype(float), reps, rng)
        with np.errstate(divide='ignore', invalid='ignore'):
            out['jaccard_lo'], out['jaccard_hi'] = ci(jac[:, 0] / jac[:, 1])
        if m.sum() > 2:
            per = np.column_stack([np.ones(m.sum()), (xm ** 2).sum(1), xm.sum(1) ** 2, xm])
            out['icc_lo'], out['icc_hi'] = ci(icc_from_stats(block_bootstrap(blocks[m], per, reps, rng), 2))
    return out, d, px, py


def discordant(sub, x, y, min_delta, top):
    both = sub[f'cov_{x}'].notna() & sub[f'cov_{y}'].notna()
    d = sub[both].copy()
    d['delta'] = d[f'pct_{y}'] - d[f'pct_{x}']
    d = d[d['delta'].abs() >= min_delta]
    if d.empty:
        return d
    d = d.reindex(d['delta'].abs().sort_values(ascending=False).index).head(top)
    out = d[['chrom', 'start', 'strand']].copy()
    out.insert(1, 'pos_1based', d['start'] + 1)
    for r in (x, y):
        out[f'cov_{r}'] = d[f'cov_{r}'].astype(int)
        out[f'pct_{r}'] = d[f'pct_{r}']
        out[f'mismatch_frac_{r}'] = (d[f'ndiff_{r}'] / (d[f'cov_{r}'] + d[f'ndiff_{r}'])).round(4)
    out['delta_y_minus_x'] = d['delta']
    return out.drop(columns='start')


def three_way(sub, runs, min_pct):
    cov = np.column_stack([sub[f'cov_{r}'].notna().to_numpy() for r in runs])
    pct = np.column_stack([sub[f'pct_{r}'].to_numpy() for r in runs])
    valid_any = np.column_stack([(sub[f'pct_{r}'] >= min_pct).to_numpy() & cov[:, i] for i, r in enumerate(runs)])
    regions = {}
    for mode, keep in (('comparable', cov.all(axis=1)), ('unrestricted', np.ones(len(sub), bool))):
        v = valid_any[keep]
        bits = (v * (1 << np.arange(v.shape[1]))).sum(axis=1)
        regions[mode] = np.bincount(bits, minlength=1 << v.shape[1])
    keep = cov.all(axis=1)
    v = valid_any[keep]
    union = v.any(axis=1)
    j3 = v.all(axis=1).sum() / union.sum() if union.sum() else np.nan
    icc3 = float(icc_from_stats(icc_stats(pct[keep][union]), len(runs))) if union.sum() > 2 else np.nan
    return regions, int(keep.sum()), j3, icc3


# ---------------------------------------------------------------- drawing
def draw_venn(ax, counts, labels, title, note):
    """counts: bincount over membership bits (bit i = run i). Uses matplotlib_venn when installed."""
    if len(labels) == 3:
        sizes = tuple(int(counts[b]) for b in (1, 2, 3, 4, 5, 6, 7))      # 100 010 110 001 101 011 111
        totals = [int(sum(counts[b] for b in range(8) if b & (1 << i))) for i in range(3)]
        set_labels = [f'{l}\nTotal: {t:,}' for l, t in zip(labels, totals)]
        try:
            from matplotlib_venn import venn3
            v = venn3(subsets=sizes, set_labels=set_labels, ax=ax)
            for patch, color in zip(v.patches, VENN_PATCH_COLORS):
                if patch:
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
            for t in v.set_labels or []:
                if t:
                    t.set_fontsize(8)
            for t in v.subset_labels or []:
                if t:
                    t.set_fontsize(8)
        except ImportError:                        # plain, not area-proportional
            centers = [(-0.35, 0.2), (0.35, 0.2), (0.0, -0.4)]
            for (cx, cy), color in zip(centers, VENN_PATCH_COLORS[:2] + [VENN_PATCH_COLORS[3]]):
                ax.add_patch(Circle((cx, cy), 0.6, facecolor=color, alpha=0.4, edgecolor='k', lw=0.6))
            pos = {1: (-0.6, 0.35), 2: (0.6, 0.35), 3: (0, 0.45), 4: (0, -0.75), 5: (-0.33, -0.25),
                   6: (0.33, -0.25), 7: (0, 0.0)}
            for b, (px, py) in pos.items():
                ax.text(px, py, f'{int(counts[b]):,}', ha='center', va='center', fontsize=8)
            for (cx, cy), l in zip([(-0.75, 0.95), (0.75, 0.95), (0, -1.15)], set_labels):
                ax.text(cx, cy, l, ha='center', va='center', fontsize=8)
            ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.35, 1.25); ax.set_aspect('equal'); ax.axis('off')
    else:
        ax.text(0.5, 0.5, 'Venn drawn for 3 runs only', ha='center', transform=ax.transAxes)
        ax.axis('off')
    ax.set_title(f'{title}\n{note}', fontsize=9.5)


def heat_panel(ax, px, py, xl, yl, met, min_pct, bins, norm, cmap):
    h, xe, ye = np.histogram2d(px, py, bins=[bins, bins])
    im = ax.imshow(np.ma.masked_where(h == 0, h).T, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                   norm=norm, cmap=cmap, aspect='equal', interpolation='nearest')
    ax.plot([0, 100], [0, 100], color='k', ls='--', lw=0.8)
    ax.axhline(min_pct, color='r', ls='--', lw=0.8)
    ax.axvline(min_pct, color='r', ls='--', lw=0.8)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel(f'{xl} occupancy (%)', fontsize=9)
    ax.set_ylabel(f'{yl} occupancy (%)', fontsize=9)
    for s in ax.spines.values():
        s.set_edgecolor('black')
    return im, (f"n = {met['n_comparable']:,}   r = {met['r_all']:.3f}\n"
                f"ICC(2,1) = {met['icc']:.3f} (n = {met['n_icc']:,})   Jaccard = {met['jaccard']:.3f}")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--run', action='append', required=True, metavar='NAME=BED')
    ap.add_argument('--out', required=True)
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'(_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-)')
    ap.add_argument('--min-cov', type=int, default=20)
    ap.add_argument('--min-pct', type=float, default=20)
    ap.add_argument('--bootstrap', type=int, default=200, help='block-bootstrap replicates for CIs (0: none)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--style')
    ap.add_argument('--discordant-delta', type=float, default=50, help='list sites whose occupancy differs by >= this (points)')
    ap.add_argument('--discordant-top', type=int, default=200, help='at most this many per modification and pair')
    a = ap.parse_args()

    runs, paths = [], {}
    for item in a.run:
        name, _, path = item.partition('=')
        if not path or not os.path.exists(path):
            sys.exit(f'cmp_plots.py: --run {item}: no such file')
        runs.append(name); paths[name] = path
    if len(runs) < 2:
        sys.exit('cmp_plots.py: need at least two --run')
    label = {r: KNOWN_LABELS.get(r, r) for r in runs}
    short = {r: SHORT.get(r, r) for r in runs}
    pairs = [(x, y, why) for x, y, why in PREFERRED_PAIRS if x in runs and y in runs]
    if not pairs:
        pairs = [(runs[i], runs[j], '') for i in range(len(runs)) for j in range(i + 1, len(runs))]
    os.makedirs(f'{a.out}/figures', exist_ok=True); os.makedirs(f'{a.out}/tables', exist_ok=True)
    style = set_style(a.style)
    rng = np.random.default_rng(a.seed)
    contig_re = re.compile(a.contig_regex) if a.contig_regex else None
    mask = read_mask(a.mask_bed)

    # ---- load and join
    merged, load_rows = None, []
    for r in runs:
        df, n_all, n_masked = read_bedmethyl(paths[r], a.min_cov, contig_re, mask)
        load_rows.append({'run': r, 'rows_in_bed': n_all, f'rows_cov_ge_{a.min_cov}_kept': len(df),
                          'rows_masked_by_bed': n_masked})
        print(f'{r}: {n_all:,} rows; {len(df):,} with N_valid >= {a.min_cov} after masking', file=sys.stderr)
        df = df.rename(columns={'cov': f'cov_{r}', 'pct': f'pct_{r}', 'ndiff': f'ndiff_{r}'})
        merged = df if merged is None else merged.merge(df, on=['chrom', 'start', 'strand', 'code'], how='outer')
    pd.DataFrame(load_rows).to_csv(f'{a.out}/tables/loaded_rows.tsv', sep='\t', index=False)

    # ---- metrics per modification
    pair_rows, venn_rows, three_rows, store, disc_parts = [], [], [], {}, []
    for code, mod in MODS:
        sub = merged[merged['code'] == code]
        if sub.empty:
            print(f'{mod}: no rows', file=sys.stderr)
            continue
        regions, n_univ, j3, icc3 = three_way(sub, runs, a.min_pct)
        for mode, cnt in regions.items():
            row = {'modification': mod, 'mode': mode}
            for b in range(1, len(cnt)):
                row['&'.join(r for i, r in enumerate(runs) if b & (1 << i))] = int(cnt[b])
            venn_rows.append(row)
        three_rows.append({'modification': mod, 'n_comparable_all_runs': n_univ, 'jaccard_all_runs': j3,
                           'icc21_all_runs': icc3})
        store[mod] = {'regions': regions, 'n_univ': n_univ, 'j3': j3, 'icc3': icc3, 'pairs': []}
        for x, y, why in pairs:
            met, d, px, py = pair_metrics(sub, x, y, a.min_pct, a.bootstrap, rng)
            dd = discordant(sub, x, y, a.discordant_delta, a.discordant_top)
            if not dd.empty:
                dd.insert(0, 'run_y', y); dd.insert(0, 'run_x', x); dd.insert(0, 'modification', mod)
                disc_parts.append(dd.rename(columns=lambda c: c.replace(f'_{x}', '_x').replace(f'_{y}', '_y')))
            pair_rows.append({'modification': mod, 'run_x': x, 'run_y': y, 'comparison': why, **met})
            store[mod]['pairs'].append((x, y, why, met, px, py))
    pd.DataFrame(pair_rows).to_csv(f'{a.out}/tables/pairwise_metrics.tsv', sep='\t', index=False, float_format='%.5g')
    pd.DataFrame(venn_rows).to_csv(f'{a.out}/tables/venn_counts.tsv', sep='\t', index=False)
    if disc_parts:
        pd.concat(disc_parts, ignore_index=True).to_csv(f'{a.out}/tables/discordant_sites.tsv', sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(three_rows).to_csv(f'{a.out}/tables/all_runs_metrics.tsv', sep='\t', index=False, float_format='%.5g')
    thr = {r: read_thresholds(paths[r]) for r in runs}
    if any(thr.values()):
        pd.DataFrame([{'run': r, **(t or {})} for r, t in thr.items()]).to_csv(
            f'{a.out}/tables/thresholds.tsv', sep='\t', index=False)
    mods_present = [m for _, m in MODS if m in store]

    # ---- Venns
    for mode in ('comparable', 'unrestricted'):
        fig, axes = plt.subplots(2, 4, figsize=(18, 10.5))
        for ax, mod in zip(axes.flat, [m for _, m in MODS]):
            if mod not in store:
                ax.axis('off'); continue
            s = store[mod]
            note = (f"Jaccard = {s['j3']:.3f} | ICC(2,1) = {s['icc3']:.3f}\ncomparable sites {s['n_univ']:,}"
                    if mode == 'comparable' else "each run's valid sites, any coverage elsewhere")
            draw_venn(ax, s['regions'][mode], [short[r] for r in runs], mod, note)
        fig.suptitle(f"{a.sample}: valid sites (N_valid >= {a.min_cov}, >= {a.min_pct:g}% modified), "
                     f"{'comparable set: covered in every run' if mode == 'comparable' else 'unrestricted'}", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(f'{a.out}/figures/venn_{mode}.pdf'); plt.close(fig)

    # ---- heatmaps: one page (and one PDF) per modification, one panel per pair
    cmap = matplotlib.colormaps['viridis'].copy(); cmap.set_bad('white', 1.0)
    norm = mcolors.LogNorm(vmin=1, vmax=1e6)
    bins = np.arange(0, 105, 5)
    with PdfPages(f'{a.out}/figures/heatmaps.pdf') as pdf:
        for mod in mods_present:
            ps = store[mod]['pairs']
            fig, axes = plt.subplots(1, len(ps), figsize=(4.6 * len(ps) + 1.2, 5.4), squeeze=False, layout='constrained')
            im = None
            for ax, (x, y, why, met, px, py) in zip(axes[0], ps):
                im, mtxt = heat_panel(ax, px, py, short[x], short[y], met, a.min_pct, bins, norm, cmap)
                ax.set_title(f"{why or short[x] + ' vs ' + short[y]}\n{mtxt}", fontsize=8.5)
            fig.colorbar(im, ax=axes[0].tolist(), shrink=0.75, fraction=0.03, pad=0.02, label='Sites (log scale)')
            fig.suptitle(f'{a.sample} {mod}: occupancy at sites with N_valid >= {a.min_cov} in both runs', fontsize=12)
            pdf.savefig(fig); fig.savefig(f'{a.out}/figures/heatmap_{mod}.pdf'); plt.close(fig)

    # ---- summary: valid sites per run, Jaccard and ICC per pair, thresholds
    pr = pd.DataFrame(pair_rows)
    has_thr = any(thr.values())
    fig, axes = plt.subplots(1, 4 if has_thr else 3, figsize=(22 if has_thr else 17, 6.2))
    ax = axes[0]
    xs = np.arange(len(mods_present)); w = 0.8 / len(runs)
    for i, r in enumerate(runs):
        vals = [store[m]['regions']['unrestricted'][[b for b in range(1, 1 << len(runs)) if b & (1 << i)]].sum()
                for m in mods_present]
        ax.bar(xs + (i - (len(runs) - 1) / 2) * w, vals, w, color=RUN_COLORS[i], edgecolor='k', lw=0.4, label=short[r])
    ax.set_xticks(xs); ax.set_xticklabels(mods_present, rotation=45, ha='right')
    ax.set_yscale('log'); ax.set_ylabel('Valid sites'); ax.set_title('Valid sites per run'); ax.legend(fontsize=8, frameon=False)
    for ax, metric, ttl in ((axes[1], 'jaccard', 'Jaccard of valid sites'), (axes[2], 'icc', 'ICC(2,1) of occupancy')):
        for j, (x, y, why) in enumerate(pairs):
            d = pr[(pr.run_x == x) & (pr.run_y == y)].set_index('modification').reindex(mods_present)
            off = (j - (len(pairs) - 1) / 2) * 0.22
            yerr = np.vstack([d[metric] - d[f'{metric}_lo'], d[f'{metric}_hi'] - d[metric]])
            ax.errorbar(xs + off, d[metric], yerr=np.where(np.isfinite(yerr), yerr, 0), fmt='o', ms=5,
                        color=PAIR_COLORS[j % 3], capsize=2, lw=1, label=f'{short[x]} vs {short[y]}' + (f'\n({why})' if why else ''))
        ax.set_xticks(xs); ax.set_xticklabels(mods_present, rotation=45, ha='right')
        ax.set_ylim(0, 1.02); ax.set_title(ttl + (' (95% CI)' if a.bootstrap else ''))
        ax.legend(fontsize=7, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=1)
    if has_thr:
        ax = axes[3]; bases = ['A', 'C', 'G', 'T']; xb = np.arange(4)
        for i, r in enumerate(runs):
            t = thr[r] or {}
            ax.bar(xb + (i - (len(runs) - 1) / 2) * w, [t.get(b, np.nan) for b in bases], w,
                   color=RUN_COLORS[i], edgecolor='k', lw=0.4, label=short[r])
        ax.set_xticks(xb); ax.set_xticklabels(['A', 'C', 'G', 'U']); ax.set_ylim(0.5, 1.0)
        ax.set_ylabel('Call threshold'); ax.set_title('Per-sample thresholds (10th percentile)'); ax.legend(fontsize=8, frameon=False)
    for ax in axes:
        for s in ax.spines.values():
            s.set_edgecolor('black')
    fig.suptitle(f'{a.sample}: summary', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f'{a.out}/figures/summary.pdf'); plt.close(fig)

    # ---- legends
    runs_txt = '\n'.join(f'- **{short[r]}** (`{r}`): {label[r].replace(chr(10), " ")}; `{paths[r]}`' for r in runs)
    with open(f'{a.out}/legends.md', 'w') as fh:
        fh.write(f"""# {a.sample}: basecalling comparison (same reads in every run)

Runs:
{runs_txt}

Settings: valid site = N_valid_cov >= {a.min_cov} and raw modkit percent >= {a.min_pct:g}; no gIVT correction.
Masked: contigs matching `{a.contig_regex}`{'; intervals in `' + a.mask_bed + '`' if a.mask_bed else ''}.
Style: {style}. CIs: {a.bootstrap} block-bootstrap replicates over {BLOCK // 1000} kb genomic blocks{' (none)' if not a.bootstrap else ''}.

**venn_comparable.pdf**: per modification, valid sites among sites covered (N_valid >= {a.min_cov}) in every run, so a
difference is a different call, not different coverage. Below each: Jaccard (sites valid in all runs / valid in any)
and ICC(2,1) over all runs, on comparable sites valid in at least one run.

**venn_unrestricted.pdf**: each run's valid sites, whatever the other runs' coverage at those sites.

**heatmaps.pdf** (and heatmap_<mod>.pdf): per modification, one panel per pair of runs; every site with N_valid >= {a.min_cov}
in both, binned at 5% occupancy, log colour scale; dashed y = x, red lines at {a.min_pct:g}%. Text: n sites, Pearson r
over all of them, ICC(2,1) over sites valid in either run, Jaccard of valid sites.

**summary.pdf**: valid sites per run (log scale); Jaccard and ICC(2,1) per pair with 95% block-bootstrap CIs; the
per-sample call thresholds (A, C, G, U) each run's pileup used.

Tables (tables/): pairwise_metrics.tsv, all_runs_metrics.tsv, venn_counts.tsv, thresholds.tsv, loaded_rows.tsv,
discordant_sites.tsv (per modification and pair, up to {a.discordant_top} sites whose occupancy differs by >= {a.discordant_delta:g}
points, with each run's coverage, occupancy and mismatch fraction N_diff / (N_valid + N_diff), 1-based position).
""")
    print(f'wrote {a.out}/figures, {a.out}/tables, {a.out}/legends.md', file=sys.stderr)


if __name__ == '__main__':
    main()
