#!/usr/bin/env python3
"""make_report_figures.py: figures for the v6 vs v5.2 report (Fig1-Fig9, FigS1).

    python make_report_figures.py --base /scratch/stein.an/basecall_cmp --out <report>/figures

Fig1 answers one question (does the Dorado version matter, same v5.2 models?); every other figure compares v6.0.0 with
v5.2.0, both on Dorado 2.1.2, and never shows the Dorado 1.3.0 run. Reads only analysis outputs:
  basecall_qc/GM12878/tables/per_read_<run>.tsv.gz      Fig2        confusion/<s>/tables/*.tsv            Fig1A, Fig3A-B
  plots/<s>_<arm>/tables/{pairwise_metrics,discordant_sites}.tsv, plots/HEK293_replicates_<run>_<arm>/   Fig1B, Fig3C-D, Fig8
  ortho/<s>_<arm>/tables/{ci,pr_curves,summary}.tsv     Fig1C, Fig5, Fig6, Fig7A
  truthsets/<s>_<arm>_tol<t>/tables/{pairs,sets}.tsv    Fig4, FigS1 motif/<s>_<arm>/tables/motif.tsv     Fig7B-C
  sampleprobs_050/...thresholds.tsv                      Fig9
Each figure is built on its own; a missing input skips that figure. Writes Fig*.pdf/png and figure_numbers.md.
"""
import argparse
import os
import re
import sys
import traceback

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

HEK = ['08_01_25_RNA004_HEK293', '09_09_25_RNA004_HEK293', '10_06_25_RNA004_HEK293_CP']
SAMPLES = ['GM12878'] + HEK
SLAB = {'GM12878': 'GM12878', HEK[0]: '08/01', HEK[1]: '09/09', HEK[2]: '10/06'}
V6, V52, V52OLD = 'A_d212_v6', 'C_d212_v52', 'B_d130_v52'
COL = {V6: '#d95f02', V52: '#1b9e77'}
LAB = {V6: 'v6.0.0', V52: 'v5.2.0'}
ARMS = {'nofilter': 'no threshold', 'lab050': 'lab rule'}
MODS = ['m6A', 'Inosine', '2OmeA', 'm5C', '2OmeC', 'Psi', '2OmeU', '2OmeG']
MLAB = {'m6A': 'm$^6$A', 'Inosine': 'Inosine', '2OmeA': "2$'$-OMe\nA", 'm5C': 'm$^5$C', '2OmeC': "2$'$-OMe\nC", 'Psi': r'$\Psi$', '2OmeU': "2$'$-OMe\nU",
        '2OmeG': "2$'$-OMe\nG"}
MAIN = [('m6A', 'GLORI1.0∩2.0'), ('Inosine', 'SLIC-seq_2of3'), ('Psi', 'BID∩PRAISE'), ('m5C', 'UBS-seq')]
SETLAB = {'GLORI1.0∩2.0': r'GLORI 1.0 $\cap$ 2.0', 'SLIC-seq_2of3': 'SLIC-seq, 2+ replicates', 'BID∩PRAISE': r'BID-seq $\cap$ PRAISE',
          'UBS-seq': 'UBS-seq'}
# every truth set in the paired accuracy plot (Fig6), by modification; DART-seq is left out (no base-resolution agreement)
PAIRED = {'m6A': ['GLORI1.0∩2.0', 'm6A-SAC-seq', 'miCLIP-CITS', 'miCLIP-CIMS', 'm6A_2+methods'], 'Inosine': ['SLIC-seq_2of3'],
          'Psi': ['BID∩PRAISE', 'PRAISE_PUS-dependent'], 'm5C': ['UBS-seq', 'UBS-seq_NSUN2-dependent']}
MODCOL = {'m6A': '#3b6ea8', 'Inosine': '#c2508a', 'Psi': '#6a9a3a', 'm5C': '#8c6bb1'}
# concordance matrices: display names; anything not listed (DART-seq, the Dorado 1.3.0 run, derived subsets) is left out
MATRIX = {'m6A': ['GLORI1.0', 'GLORI2.0', 'GLORI1.0_2022', 'm6A-SAC-seq', 'miCLIP-CITS', 'miCLIP-CIMS', 'nanopore v6', 'nanopore v5.2 (2.1.2)'],
          'Psi': ['BID-seq', 'PRAISE', 'nanopore v6', 'nanopore v5.2 (2.1.2)'],
          'Inosine': ['SLIC-seq', 'nanopore v6', 'nanopore v5.2 (2.1.2)'], 'm5C': ['UBS-seq', 'nanopore v6', 'nanopore v5.2 (2.1.2)']}
DISP = {'GLORI1.0': 'GLORI 1.0, 2025', 'GLORI2.0': 'GLORI 2.0, 2025', 'GLORI1.0_2022': 'GLORI 1.0, 2022',
        'miCLIP-CITS': 'miCLIP (SySy)', 'miCLIP-CIMS': 'miCLIP (Abcam)', 'nanopore v6': 'v6.0.0', 'nanopore v5.2 (2.1.2)': 'v5.2.0'}
W = 7.2
missing, numbers = [], []


def read(path, **kw):
    if not os.path.exists(path):
        missing.append(path)
        return None
    return pd.read_csv(path, sep='\t', **kw)


def label(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')


def save(fig, out, name):
    fig.savefig(f'{out}/{name}.pdf'); fig.savefig(f'{out}/{name}.png', dpi=300); plt.close(fig)
    print(f'  {name}', file=sys.stderr)


def span(v):
    return [np.clip(v['mean'] - v['min'], 0, None), np.clip(v['max'] - v['mean'], 0, None)]


def thresholds(path):
    if not os.path.exists(path):
        missing.append(path)
        return {}
    out = {}
    for line in open(path):
        t = line.split()
        if t and len(t[0]) == 1 and t[0] in 'ACGT':
            try:
                out['U' if t[0] == 'T' else t[0]] = float(t[-1])
            except ValueError:
                pass
    return out


def call_ratio(s, x, y):
    """calls per modification, run x divided by run y, from the read-level confusion table (identical reads)."""
    t = read(f'{B}/confusion/{s}/tables/confusion_{x}_vs_{y}.tsv')
    if t is None:
        return {}
    return {m: t.loc[t['call_x'] == m, 'n'].sum() / max(t.loc[t['call_y'] == m, 'n'].sum(), 1) for m in MODS}


def ratio_bars(ax, x, y, color):
    r = pd.DataFrame([dict(call_ratio(s, x, y), sample=s) for s in SAMPLES]).set_index('sample')
    g = np.log2(r[MODS]).agg(['mean', 'min', 'max']).T
    xs = np.arange(len(MODS))
    ax.bar(xs, g['mean'], 0.6, color=color)
    ax.errorbar(xs, g['mean'], yerr=span(g), fmt='none', ecolor='k', lw=0.7, capsize=2)
    ax.axhline(0, color='k', lw=0.5); ax.set_xticks(xs); ax.set_xticklabels([MLAB[m] for m in MODS], fontsize=6.3)
    ax.set_ylim(-1.3, 1.3); ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticklabels(['0.5×', '0.71×', '1×', '1.41×', '2×'])
    return r


def icc(pairs, rx, ry):
    rows = []
    for arm in ARMS:
        for s in pairs:
            t = read(f'{B}/plots/{s}_{arm}/tables/pairwise_metrics.tsv')
            if t is not None:
                rows.append(t[(t['run_x'] == rx) & (t['run_y'] == ry)][['modification', 'icc']].assign(arm=arm, sample=s))
    return pd.concat(rows) if rows else pd.DataFrame(columns=['modification', 'icc', 'arm', 'sample'])


# ================================================================ Fig 1: does the Dorado version matter?
def fig1(out):
    fig = plt.figure(figsize=(W, 5.0), layout='constrained')
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'A')
    r = ratio_bars(ax, V52, V52OLD, '0.45')
    ax.set_ylabel('calls per read, 2.1.2 ÷ 1.3.0'); ax.set_title('Calls per read, same reads', fontsize=8)
    numbers.append(('Fig1A call ratio Dorado 2.1.2 / 1.3.0 (min-max)', r[MODS].agg(['min', 'max']).round(3).to_string()))

    ax = fig.add_subplot(gs[0, 1]); label(ax, 'B')
    t = icc(SAMPLES, V52, V52OLD)
    xs = np.arange(len(MODS)); w = 0.38
    for i, (arm, shade) in enumerate((('nofilter', '0.35'), ('lab050', '0.7'))):
        g = t[t['arm'] == arm].groupby('modification')['icc'].agg(['mean', 'min', 'max']).reindex(MODS)
        ax.bar(xs + (i - 0.5) * w, g['mean'], w, color=shade, label=ARMS[arm])
        ax.errorbar(xs + (i - 0.5) * w, g['mean'], yerr=span(g), fmt='none', ecolor='k', lw=0.6, capsize=1.5)
    ax.set_xticks(xs); ax.set_xticklabels([MLAB[m] for m in MODS]); ax.set_ylim(0, 1.25); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel('site ICC(2,1)'); ax.set_title('Site agreement, 2.1.2 against 1.3.0', fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc='upper center', ncol=2)

    ax = fig.add_subplot(gs[1, :]); label(ax, 'C', x=-0.06)
    rows = []
    for arm in ARMS:
        for s in HEK:
            ci = read(f'{B}/ortho/{s}_{arm}/tables/ci.tsv')
            if ci is not None:
                rows.append(ci[ci['comparison'] == f'{V52} - {V52OLD}'].assign(arm=arm, sample=s))
    if rows:
        v = pd.concat(rows); y = 0; ticks = []
        for mod, st in MAIN:
            for metric, mlab in (('precision_at_20', 'precision at 20%'), ('auprc', 'AUPRC')):
                for a_i, (arm, mk, c) in enumerate((('nofilter', 'o', '0.3'), ('lab050', 's', '0.6'))):
                    for j, s in enumerate(HEK):
                        q = v[(v['set'] == st) & (v['metric'] == metric) & (v['arm'] == arm) & (v['sample'] == s)]
                        if len(q):
                            q = q.iloc[0]
                            ax.errorbar(q['estimate'], y + a_i * 0.36 + j * 0.1, xerr=[[max(q['estimate'] - q['ci_low'], 0)], [max(q['ci_high'] - q['estimate'], 0)]],
                                        fmt=mk, ms=3, color=c, lw=0.7)
                ticks.append((y + 0.28, f'{MLAB[mod]}, {mlab}')); y += 1
        ax.axvline(0, color='k', lw=0.5); ax.set_xlim(-0.1, 0.1)
        ax.set_yticks([p for p, _ in ticks]); ax.set_yticklabels([l for _, l in ticks], fontsize=7); ax.invert_yaxis()
        ax.set_xlabel('Dorado 2.1.2 minus 1.3.0 (95% paired block bootstrap interval)')
        ax.plot([], [], 'o', color='0.3', ms=3.5, label='no threshold'); ax.plot([], [], 's', color='0.6', ms=3.5, label='lab rule')
        ax.legend(frameon=False, fontsize=7, loc='lower right', title='three HEK293 runs each', title_fontsize=6.5)
        numbers.append(('Fig1C max |difference|, Dorado 2.1.2 - 1.3.0',
                        v[v['set'].isin([s_ for _, s_ in MAIN]) & v['metric'].isin(['precision_at_20', 'auprc'])]['estimate'].abs().max()))
    save(fig, out, 'Fig1')


# ================================================================ Fig 2: basecalling, v6 against v5.2
def fig2(out, s='GM12878'):
    t = {r: read(f'{B}/basecall_qc/{s}/tables/per_read_{r}.tsv.gz') for r in (V6, V52)}
    if any(v is None for v in t.values()):
        return
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.7), layout='constrained')
    ax = axes[0]; label(ax, 'A')
    bins = np.linspace(0.90, 1.0, 101)
    for r in (V6, V52):
        ax.hist(t[r]['identity'], bins=bins, histtype='step', color=COL[r], lw=1.1, density=True, label=LAB[r])
    ax.set_xlabel('read identity'); ax.set_ylabel('density'); ax.set_xlim(0.90, 1.0)
    ax.legend(frameon=False, fontsize=7, loc='upper left'); ax.set_title(s, fontsize=8)
    ax = axes[1]; label(ax, 'B')
    m = t[V6][['read_id', 'identity']].merge(t[V52][['read_id', 'identity']], on='read_id', suffixes=('_v6', '_v52'))
    d = (m['identity_v6'] - m['identity_v52']) * 100
    inside = d[d.abs() <= 2]
    ax.hist(inside, bins=np.linspace(-2, 2, 81), color='0.55')
    ax.axvline(0, color='k', lw=0.5)
    better, worse, same = (d > 0).mean() * 100, (d < 0).mean() * 100, (d == 0).mean() * 100
    ax.text(0.98, 0.96, f'v6 better: {better:.0f}%\nv6 worse: {worse:.0f}%\nidentical: {same:.0f}%\n'
                        f'outside ±2: {100 * (1 - len(inside) / len(d)):.1f}%', transform=ax.transAxes, ha='right', va='top', fontsize=7)
    ax.set_xlabel('identity, v6 minus v5.2 (points, same read)'); ax.set_ylabel('reads'); ax.set_xlim(-2, 2)
    numbers.append(('Fig2B shares (%)', f'better {better:.1f}, worse {worse:.1f}, identical {same:.1f}, median {d.median():+.3f}'))
    save(fig, out, 'Fig2')


# ================================================================ Fig 3: how the calls change
def fig3(out):
    fig, axes = plt.subplots(2, 2, figsize=(W, 5.0), layout='constrained')
    ax = axes[0, 0]; label(ax, 'A')
    r = ratio_bars(ax, V6, V52, COL[V6])
    ax.set_ylabel('calls per read, v6 ÷ v5.2'); ax.set_title('Calls per read, same reads', fontsize=8)
    numbers.append(('Fig3A call ratio v6 / v5.2 (min-max)', r[MODS].agg(['min', 'max']).round(3).to_string()))
    ax = axes[0, 1]; label(ax, 'B')
    c = read(f'{B}/confusion/GM12878/tables/agreement_by_confidence.tsv')
    if c is not None:
        c = c[c['pair'] == f'{V6} vs {V52}']
        for base, col in zip('ACGT', ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00']):
            d = c[c['base'] == base]
            mid = d['call_prob_x_bin'].map(lambda b: np.mean([float(v) for v in re.findall(r'[\d.]+', b)][:2]))
            ax.plot(mid, d['identical_share'] * 100, '-o', ms=2.5, lw=1, color=col, label='U' if base == 'T' else base)
    ax.set_xlabel('v6 call probability'); ax.set_ylabel('v5.2 makes the same call (%)'); ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=7, title='base', title_fontsize=7, loc='lower right'); ax.set_title('Agreement by confidence, GM12878', fontsize=8)
    same = icc(SAMPLES, V6, V52)
    reps = []
    for rr in (V6, V52):
        for arm in ARMS:
            t = read(f'{B}/plots/HEK293_replicates_{rr}_{arm}/tables/pairwise_metrics.tsv')
            if t is not None:
                reps.append(t.assign(model=rr, arm=arm))
    reps = pd.concat(reps) if reps else pd.DataFrame(columns=['modification', 'icc', 'model', 'arm'])
    xs = np.arange(len(MODS)); w = 0.27
    for ax, arm, lab_ in ((axes[1, 0], 'nofilter', 'C'), (axes[1, 1], 'lab050', 'D')):
        label(ax, lab_)
        series = [('v6 vs v5.2, same reads', same[same['arm'] == arm].groupby('modification')['icc'].mean(), '0.25'),
                  ('v6, run vs run', reps[(reps['arm'] == arm) & (reps['model'] == V6)].groupby('modification')['icc'].mean(), COL[V6]),
                  ('v5.2, run vs run', reps[(reps['arm'] == arm) & (reps['model'] == V52)].groupby('modification')['icc'].mean(), COL[V52])]
        for i, (name, v, col) in enumerate(series):
            ax.bar(xs + (i - 1) * w, v.reindex(MODS), w, color=col, label=name)
        ax.axhline(0, color='k', lw=0.5); ax.set_xticks(xs); ax.set_xticklabels([MLAB[m] for m in MODS]); ax.set_ylim(-0.3, 1.05)
        ax.set_ylabel('site ICC(2,1)'); ax.set_title(f'Site agreement, {ARMS[arm]}', fontsize=8)
        if arm == 'nofilter':
            ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    save(fig, out, 'Fig3')


# ================================================================ Fig 4 / S1: how well the maps agree
def matrix(ax, p, sz, mod, title, fs):
    names = [n for n in MATRIX[mod] if n in set(p['set_a']) | set(p['set_b'])]
    M = np.full((len(names), len(names)), np.nan)
    for _, r in p[p['set_a'].isin(names) & p['set_b'].isin(names)].iterrows():
        i, j = names.index(r['set_a']), names.index(r['set_b'])
        M[i, j], M[j, i] = r['share_a_in_b'], r['share_b_in_a']
    im = ax.imshow(M * 100, cmap='viridis', vmin=0, vmax=100)
    for i in range(len(names)):
        for j in range(len(names)):
            if i != j and np.isfinite(M[i, j]):
                ax.text(j, i, f'{100 * M[i, j]:.0f}', ha='center', va='center', fontsize=fs, color='white' if M[i, j] < 0.6 else 'black')
    nsz = sz.set_index('set')['sites']
    ax.set_xticks(range(len(names))); ax.set_xticklabels([DISP.get(n, n) for n in names], rotation=40, ha='right', fontsize=fs)
    ax.set_yticks(range(len(names))); ax.set_yticklabels([f'{DISP.get(n, n)} ({int(nsz.get(n, 0)):,})' for n in names], fontsize=fs)
    for k in range(len(names) + 1):                      # thin white grid between cells
        ax.axhline(k - 0.5, color='white', lw=0.8); ax.axvline(k - 0.5, color='white', lw=0.8)
    ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_title(title, fontsize=9)
    return im


def fig4(out, tol=0, name='Fig4', s=HEK[0]):
    d = f'{B}/truthsets/{s}_nofilter_tol{tol}/tables'
    p, sz = read(f'{d}/pairs.tsv'), read(f'{d}/sets.tsv')
    if p is None or sz is None:
        return
    p, sz = p[p['universe'] == 'detectable'], sz[sz['universe'] == 'detectable']
    fig = plt.figure(figsize=(W, 8.4), layout='constrained')
    gs = fig.add_gridspec(2, 3, height_ratios=[1.75, 1])
    ims = []
    for spec, mod, lab_, fs in ((gs[0, :], 'm6A', 'A', 7.5), (gs[1, 0], 'Psi', 'B', 7), (gs[1, 1], 'Inosine', 'C', 7), (gs[1, 2], 'm5C', 'D', 7)):
        ax = fig.add_subplot(spec); label(ax, lab_, x=-0.3 if mod != 'm6A' else -0.02)
        q, z = p[p['mod'] == mod], sz[sz['mod'] == mod]
        if len(q):
            ims.append(matrix(ax, q, z, mod, MLAB[mod], fs))
    if ims:
        fig.colorbar(ims[0], ax=fig.axes, shrink=0.35, label="% of the row's sites found in the column's set")
    save(fig, out, name)


# ================================================================ Fig 5: precision against recall
def fig5(out, s=HEK[0]):
    curves = {arm: read(f'{B}/ortho/{s}_{arm}/tables/pr_curves.tsv') for arm in ARMS}
    if all(c is None for c in curves.values()):
        return
    fig, axes = plt.subplots(2, 2, figsize=(W, 5.6), layout='constrained')
    for ax, lab_, (mod, st) in zip(axes.flat, 'ABCD', MAIN):
        label(ax, lab_)
        for arm, ls in (('nofilter', '-'), ('lab050', '--')):
            c = curves[arm]
            if c is None:
                continue
            for r in (V6, V52):
                g = c[(c['set'] == st) & (c['run'] == r) & (c['n_called'] >= 20)].sort_values('cutoff_pct')
                ax.plot(g['recall'], g['precision'], ls, color=COL[r], lw=1.2)
                p20 = g[g['cutoff_pct'] == 20]
                ax.plot(p20['recall'], p20['precision'], 'o', color=COL[r], ms=3.5, mec='k', mew=0.4)
        ax.set_title(f'{MLAB[mod]} against {SETLAB[st]}', fontsize=8); ax.set_xlabel('recall'); ax.set_ylabel('precision')
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    h = [plt.Line2D([], [], color=COL[V6], lw=1.5), plt.Line2D([], [], color=COL[V52], lw=1.5),
         plt.Line2D([], [], color='0.3', lw=1.2), plt.Line2D([], [], color='0.3', lw=1.2, ls='--')]
    axes[0, 0].legend(h, ['v6.0.0', 'v5.2.0', 'no threshold', 'lab rule'], frameon=False, fontsize=7, loc='lower left')
    save(fig, out, 'Fig5')


# ================================================================ Fig 6: v6 against v5.2 on every map
def fig6(out):
    rows = []
    for arm in ARMS:
        for s in HEK:
            sm = read(f'{B}/ortho/{s}_{arm}/tables/summary.tsv')
            if sm is not None:
                rows.append(sm.assign(arm=arm, sample=s))
    if not rows:
        return
    sm = pd.concat(rows)
    fig, axes = plt.subplots(1, 2, figsize=(W, 3.5), layout='constrained')
    for ax, lab_, metric, title in ((axes[0], 'A', 'precision_at_20', 'Precision at 20% occupancy'), (axes[1], 'B', 'auprc', 'AUPRC')):
        label(ax, lab_)
        vals = []
        for mod, sets_ in PAIRED.items():
            for st in sets_:
                for arm, mk in (('nofilter', 'o'), ('lab050', 's')):
                    q = sm[(sm['set'] == st) & (sm['arm'] == arm)].pivot_table(index='sample', columns='run', values=metric)
                    if V6 in q and V52 in q:
                        ax.scatter(q[V52], q[V6], marker=mk, s=16, color=MODCOL[mod], edgecolor='white', lw=0.4, zorder=3)
                        vals += list(q[V52]) + list(q[V6])
        vals = [v for v in vals if v > 0]
        if vals:
            lo, hi = min(vals) / 1.6, min(1, max(vals) * 1.6)
            ax.plot([lo, hi], [lo, hi], color='k', lw=0.6, zorder=1)
            ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.text(0.04, 0.96, 'v6 better', transform=ax.transAxes, va='top', fontsize=7, color='0.35')
        ax.text(0.96, 0.04, 'v5.2 better', transform=ax.transAxes, ha='right', fontsize=7, color='0.35')
        fmt = (lambda v, _: f'{100 * v:g}%') if metric == 'precision_at_20' else (lambda v, _: f'{v:g}')
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=(1.0, 3.0)))
            axis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt)); axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel('v5.2.0'); ax.set_ylabel('v6.0.0'); ax.set_title(title, fontsize=8); ax.set_aspect('equal')
    h = [plt.Line2D([], [], marker='o', ls='', color=c, ms=5) for c in MODCOL.values()] + \
        [plt.Line2D([], [], marker='o', ls='', color='0.5', ms=5), plt.Line2D([], [], marker='s', ls='', color='0.5', ms=5)]
    axes[1].legend(h, [MLAB[m] for m in MODCOL] + ['no threshold', 'lab rule'], frameon=False, fontsize=7, loc='center left', bbox_to_anchor=(1.02, 0.5))
    save(fig, out, 'Fig6')


# ================================================================ Fig 7: m6A
def fig7(out):
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.9), layout='constrained', gridspec_kw={'width_ratios': [1.15, 1, 1]})
    ax = axes[0]; label(ax, 'A', x=-0.22)
    rows = []
    for arm in ARMS:
        for s in HEK:
            c = read(f'{B}/ortho/{s}_{arm}/tables/pr_curves.tsv')
            if c is None:
                continue
            g = c[(c['set'] == 'GLORI1.0∩2.0') & (c['mod'] == 'm6A')]
            ref = g[(g['run'] == V52) & (g['cutoff_pct'] == 20)]
            ok = g[(g['run'] == V6)].sort_values('cutoff_pct')
            ok = ok[ok['precision'] >= ref['precision'].iloc[0]] if len(ref) else ok.iloc[0:0]
            if len(ref) and len(ok):
                rows.append({'arm': arm, 'sample': s, 'v52': ref['recall'].iloc[0], 'v6': ok['recall'].iloc[0], 'cut': ok['cutoff_pct'].iloc[0]})
    if rows:
        op = pd.DataFrame(rows); x = np.arange(len(op))
        ax.bar(x - 0.2, op['v52'], 0.4, color=COL[V52], label='v5.2 at 20%')
        ax.bar(x + 0.2, op['v6'], 0.4, color=COL[V6], label="v6 at v5.2's precision")
        for xi, cut in zip(x, op['cut']):
            ax.text(xi + 0.2, 0.03, f'{cut:.0f}%', rotation=90, ha='center', va='bottom', fontsize=6, color='white')
        ax.set_xticks(x); ax.set_xticklabels([f'{SLAB[s]}\n{"none" if a == "nofilter" else "lab"}' for a, s in zip(op['arm'], op['sample'])], fontsize=6.5)
        ax.set_ylim(0, 1); ax.set_ylabel(r'recall, GLORI 1.0 $\cap$ 2.0'); ax.legend(frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, 1.22))
    mot = [read(f'{B}/motif/{s}_nofilter/tables/motif.tsv') for s in HEK]
    mot = [m for m in mot if m is not None]
    classes = ['both', 'v6 only', 'v5.2 only', 'background (all covered A)']
    clab = ['both', 'v6\nonly', 'v5.2\nonly', 'all A']
    if mot:
        m = pd.concat(mot)
        for ax, lab_, cols, title in ((axes[1], 'B', [('drach_share', 'DRACH', '#4d4d4d'), ('glori_share', r'GLORI 1.0 $\cap$ 2.0', '#3b6ea8')], 'Motif and GLORI'),
                                      (axes[2], 'C', [('share_m6A-SAC-seq', 'm6A-SAC-seq', '#7fa7d9'), ('share_miCLIP (either antibody)', 'miCLIP', '#e3a33b'),
                                                      ('share_2+ methods', '2+ methods', '#b04a4a')], 'Other maps')):
            label(ax, lab_, x=-0.2)
            cols = [c for c in cols if c[0] in m.columns]
            w = 0.8 / max(len(cols), 1); xx = np.arange(len(classes))
            for j, (c, l, col) in enumerate(cols):
                g = (m.groupby('class')[c].agg(['mean', 'min', 'max']).reindex(classes) * 100)
                ax.bar(xx + (j - (len(cols) - 1) / 2) * w, g['mean'], w, color=col, label=l)
                ax.errorbar(xx + (j - (len(cols) - 1) / 2) * w, g['mean'], yerr=span(g), fmt='none', ecolor='k', lw=0.5, capsize=1.2)
            ax.set_xticks(xx); ax.set_xticklabels(clab, fontsize=7); ax.set_ylabel('% of sites'); ax.set_title(title, fontsize=8)
            ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    save(fig, out, 'Fig7')


# ================================================================ Fig 8: the STMN1 cluster
def fig8(out):
    d = read(f'{B}/plots/GM12878_nofilter/tables/discordant_sites.tsv')
    if d is None:
        return
    d = d[(d['modification'] == 'Inosine') & (d['run_x'] == V6) & (d['run_y'] == V52) & (d['chrom'] == 'chr1')
          & d['pos_1based'].between(25_900_500, 25_904_500)]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(3.4, 3.0), layout='constrained')
    for _, r in d.iterrows():
        ax.plot([0, 1], [r['pct_y'], r['pct_x']], color='0.6', lw=0.7, zorder=1)
    ax.scatter(np.zeros(len(d)), d['pct_y'], color=COL[V52], s=14, zorder=2)
    ax.scatter(np.ones(len(d)), d['pct_x'], color=COL[V6], s=14, zorder=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['v5.2.0', 'v6.0.0']); ax.set_xlim(-0.35, 1.35); ax.set_ylim(-3, 100)
    ax.set_ylabel('inosine occupancy (%)'); ax.set_title(f'STMN1, {len(d)} positions (GM12878)', fontsize=8)
    numbers.append(('Fig8 STMN1', f'{len(d)} positions; v5.2 {d["pct_y"].min():.0f}-{d["pct_y"].max():.0f}%; v6 max {d["pct_x"].max():.1f}%; '
                                  f'depth v6 median {d["cov_x"].median():.0f}, v5.2 {d["cov_y"].min()}-{d["cov_y"].max()}'))
    save(fig, out, 'Fig8')


# ================================================================ Fig 9: lab rule thresholds
def fig9(out):
    fig, ax = plt.subplots(figsize=(3.8, 2.8), layout='constrained')
    bases = ['A', 'C', 'G', 'U']
    for j, r in enumerate((V6, V52)):
        for i, b in enumerate(bases):
            v = [thresholds(f'{B}/sampleprobs_050/{r}/thresholds.tsv' if s == 'GM12878' else f'{B}/sampleprobs_050/{r}/{s}/thresholds.tsv').get(b, np.nan)
                 for s in SAMPLES]
            ax.scatter(np.full(len(v), i + (j - 0.5) * 0.3), v, s=16, color=COL[r], label=LAB[r] if i == 0 else None, zorder=2)
    ax.set_xticks(range(4)); ax.set_xticklabels(bases); ax.set_xlabel('base'); ax.set_ylabel('threshold (10th percentile)')
    ax.set_ylim(0.5, 1.0); ax.legend(frameon=False, fontsize=7, loc='lower right')
    save(fig, out, 'Fig9')


def main():
    global B
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--base', default='/scratch/stein.an/basecall_cmp')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    B = a.base
    os.makedirs(a.out, exist_ok=True)
    plt.rcParams.update({'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
    jobs = [('Fig1', lambda: fig1(a.out)), ('Fig2', lambda: fig2(a.out)), ('Fig3', lambda: fig3(a.out)), ('Fig4', lambda: fig4(a.out)),
            ('FigS1', lambda: fig4(a.out, tol=1, name='FigS1')), ('Fig5', lambda: fig5(a.out)), ('Fig6', lambda: fig6(a.out)),
            ('Fig7', lambda: fig7(a.out)), ('Fig8', lambda: fig8(a.out)), ('Fig9', lambda: fig9(a.out))]
    failed = []
    for name, f in jobs:
        try:
            f()
        except Exception:
            failed.append(name); print(f'{name} FAILED:\n{traceback.format_exc()}', file=sys.stderr)
    with open(f'{a.out}/figure_numbers.md', 'w') as fh:
        for title, body in numbers:
            fh.write(f'## {title}\n\n```\n{body}\n```\n\n')
        fh.write('## Missing inputs\n\n' + ('\n'.join(f'- {m}' for m in sorted(set(missing))) or 'none') + '\n')
    print(f'wrote {a.out}; failed: {", ".join(failed) or "none"}; missing inputs: {len(set(missing))}', file=sys.stderr)


B = None
if __name__ == '__main__':
    main()
