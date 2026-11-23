#!/usr/bin/env python3
"""make_short_figures.py: the three figures of the short v6 vs v5.2 report (no probability filter throughout).

    python make_short_figures.py --base /scratch/stein.an/basecall_cmp --out <report>/figures

Needs make_report_figures.py in the same directory (shared helpers, colours and input paths).
  Fig1  v6 against v5.2: calls per read, site agreement against run-to-run agreement
  Fig2  v6 against v5.2 on every map (precision), precision against recall for m6A and inosine
  Fig3  where the differences come from: m6A sites by class, the STMN1 inosine cluster
Only the no-threshold calls are used; the Dorado 1.3.0 run is not shown.
"""
import argparse
import os
import sys
import traceback

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_report_figures as mr  # noqa: E402

V6, V52, COL = mr.V6, mr.V52, mr.COL
XS = np.arange(len(mr.MODS))


def grouped(ax, series, ylim, ylabel):
    """bars per modification for several series of mean values, with whiskers when min/max are given."""
    w = 0.8 / len(series)
    for i, (name, g, col) in enumerate(series):
        x = XS + (i - (len(series) - 1) / 2) * w
        g = g.reindex(mr.MODS)
        ax.bar(x, g['mean'], w, color=col, label=name)
        if 'min' in g:
            ax.errorbar(x, g['mean'], yerr=mr.span(g), fmt='none', ecolor='k', lw=0.5, capsize=1.2)
    ax.axhline(0, color='k', lw=0.5); ax.set_xticks(XS); ax.set_xticklabels([mr.MLAB[m] for m in mr.MODS], fontsize=6.3)
    ax.set_ylim(*ylim); ax.set_ylabel(ylabel)


def fig2(out):
    fig, axes = plt.subplots(1, 2, figsize=(mr.W, 2.5), layout='constrained')
    ax = axes[0]; mr.label(ax, 'A')
    mr.ratio_bars(ax, V6, V52, COL[V6]); ax.set_ylabel('calls per read, v6 ÷ v5.2')
    ax.set_title('Calls per read, same reads', fontsize=8)
    ax = axes[1]; mr.label(ax, 'B')
    same = mr.icc(mr.SAMPLES, V6, V52)
    reps = [mr.read(f'{mr.B}/plots/HEK293_replicates_{r}_nofilter/tables/pairwise_metrics.tsv') for r in (V6, V52)]
    series = [('v6 vs v5.2, same reads', same[same['arm'] == 'nofilter'].groupby('modification')['icc'].agg(['mean']), '0.25')]
    for r, t in zip((V6, V52), reps):
        if t is not None:
            series.append(({V6: 'v6', V52: 'v5.2'}[r] + ', run vs run', t.groupby('modification')['icc'].agg(['mean']), COL[r]))
    grouped(ax, series, (-0.3, 1.05), 'site ICC(2,1)')
    ax.legend(frameon=False, fontsize=6.5, loc='upper right'); ax.set_title('Site agreement, no threshold', fontsize=8)
    mr.save(fig, out, 'Fig1')


def paired(ax, sm, metric):
    vals = []
    for mod, sets_ in mr.PAIRED.items():
        for st in sets_:
            for arm, mk in (('nofilter', 'o'),):
                q = sm[(sm['set'] == st) & (sm['arm'] == arm)].pivot_table(index='sample', columns='run', values=metric)
                if V6 in q and V52 in q:
                    ax.scatter(q[V52], q[V6], marker=mk, s=15, color=mr.MODCOL[mod], edgecolor='white', lw=0.4, zorder=3)
                    vals += list(q[V52]) + list(q[V6])
    vals = [v for v in vals if v > 0]
    if vals:
        lo, hi = min(vals) / 1.6, min(1, max(vals) * 1.6)
        ax.plot([lo, hi], [lo, hi], color='k', lw=0.6, zorder=1)
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=(1.0, 3.0)))
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{100 * v:g}%'))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.text(0.05, 0.95, 'v6 better', transform=ax.transAxes, va='top', fontsize=6.5, color='0.35')
    ax.text(0.95, 0.05, 'v5.2 better', transform=ax.transAxes, ha='right', fontsize=6.5, color='0.35')
    ax.set_xlabel('v5.2.0'); ax.set_ylabel('v6.0.0'); ax.set_aspect('equal')


def pr(ax, curves, st, title):
    for arm, ls in (('nofilter', '-'),):
        c = curves.get(arm)
        if c is None:
            continue
        for r in (V6, V52):
            g = c[(c['set'] == st) & (c['run'] == r) & (c['n_called'] >= 20)].sort_values('cutoff_pct')
            ax.plot(g['recall'], g['precision'], ls, color=COL[r], lw=1.1)
            p20 = g[g['cutoff_pct'] == 20]
            ax.plot(p20['recall'], p20['precision'], 'o', color=COL[r], ms=3.2, mec='k', mew=0.4)
    ax.set_xlim(left=0); ax.set_ylim(bottom=0); ax.set_xlabel('recall'); ax.set_ylabel('precision'); ax.set_title(title, fontsize=8)


def fig3(out, s=mr.HEK[0]):
    rows = []
    for arm in ('nofilter',):
        for smp in mr.HEK:
            t = mr.read(f'{mr.B}/ortho/{smp}_{arm}/tables/summary.tsv')
            if t is not None:
                rows.append(t.assign(arm=arm, sample=smp))
    curves = {'nofilter': mr.read(f'{mr.B}/ortho/{s}_nofilter/tables/pr_curves.tsv')}
    fig, axes = plt.subplots(1, 3, figsize=(mr.W, 2.75), layout='constrained', gridspec_kw={'width_ratios': [1.05, 1, 1]})
    mr.label(axes[0], 'A', x=-0.3)
    if rows:
        paired(axes[0], pd.concat(rows), 'precision_at_20')
        axes[0].set_title('Precision at 20%, every map', fontsize=8)
        h = [plt.Line2D([], [], marker='o', ls='', color=c, ms=4) for c in mr.MODCOL.values()]
        axes[0].legend(h, [mr.MLAB[m] for m in mr.MODCOL], frameon=False, fontsize=6, loc='lower right', bbox_to_anchor=(1.0, 0.1),
                       handletextpad=0.1, borderpad=0.2)
    for ax, lab_, st, title in ((axes[1], 'B', 'GLORI1.0∩2.0', r'm$^6$A against GLORI 1.0 $\cap$ 2.0'),
                                (axes[2], 'C', 'SLIC-seq_2of3', 'Inosine against SLIC-seq')):
        mr.label(ax, lab_, x=-0.2)
        pr(ax, curves, st, title)
    h = [plt.Line2D([], [], color=COL[V6], lw=1.4), plt.Line2D([], [], color=COL[V52], lw=1.4)]
    axes[1].legend(h, ['v6.0.0', 'v5.2.0'], frameon=False, fontsize=7, loc='lower left')
    mr.save(fig, out, 'Fig2')


def fig4(out):
    fig, axes = plt.subplots(1, 2, figsize=(mr.W, 2.7), layout='constrained', gridspec_kw={'width_ratios': [2.1, 1]})
    ax = axes[0]; mr.label(ax, 'A', x=-0.08)
    mot = [mr.read(f'{mr.B}/motif/{s}_nofilter/tables/motif.tsv') for s in mr.HEK]
    mot = [m for m in mot if m is not None]
    classes = ['both', 'v6 only', 'v5.2 only', 'background (all covered A)']
    if mot:
        m = pd.concat(mot)
        cols = [('drach_share', 'DRACH motif', '#4d4d4d'), ('glori_share', r'GLORI 1.0 $\cap$ 2.0', '#3b6ea8'),
                ('share_m6A-SAC-seq', 'm6A-SAC-seq', '#7fa7d9'), ('share_miCLIP (either antibody)', 'miCLIP', '#e3a33b')]
        cols = [c for c in cols if c[0] in m.columns]
        w = 0.8 / len(cols); xx = np.arange(len(classes))
        for j, (c, l, col) in enumerate(cols):
            g = m.groupby('class')[c].agg(['mean', 'min', 'max']).reindex(classes) * 100
            x = xx + (j - (len(cols) - 1) / 2) * w
            ax.bar(x, g['mean'].clip(lower=0.05), w, color=col, label=l)
            ax.errorbar(x, g['mean'].clip(lower=0.05), yerr=mr.span(g), fmt='none', ecolor='k', lw=0.5, capsize=1.2)
        ax.set_yscale('log'); ax.set_ylim(0.05, 100)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}%'))
        ax.set_xticks(xx); ax.set_xticklabels(['both models', 'v6 only', 'v5.2 only', 'all covered A'])
        ax.set_ylabel('share of sites'); ax.legend(frameon=False, fontsize=6.5, ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.13))
        ax.set_title(r'm$^6$A sites, no threshold', fontsize=8, pad=16)
    ax = axes[1]; mr.label(ax, 'B', x=-0.25)
    d = mr.read(f'{mr.B}/plots/GM12878_nofilter/tables/discordant_sites.tsv')
    if d is not None:
        d = d[(d['modification'] == 'Inosine') & (d['run_x'] == V6) & (d['run_y'] == V52) & (d['chrom'] == 'chr1')
              & d['pos_1based'].between(25_900_500, 25_904_500)]
        for _, r in d.iterrows():
            ax.plot([0, 1], [r['pct_y'], r['pct_x']], color='0.6', lw=0.6, zorder=1)
        ax.scatter(np.zeros(len(d)), d['pct_y'], color=COL[V52], s=12, zorder=2)
        ax.scatter(np.ones(len(d)), d['pct_x'], color=COL[V6], s=12, zorder=2)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['v5.2.0', 'v6.0.0']); ax.set_xlim(-0.35, 1.35); ax.set_ylim(-3, 100)
        ax.set_ylabel('inosine occupancy (%)'); ax.set_title(f'STMN1, {len(d)} positions', fontsize=8)
    mr.save(fig, out, 'Fig3')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--base', default='/scratch/stein.an/basecall_cmp')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    mr.B = a.base
    os.makedirs(a.out, exist_ok=True)
    plt.rcParams.update({'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.linewidth': 0.6})
    failed = []
    for name, f in (('Fig1', fig2), ('Fig2', fig3), ('Fig3', fig4)):
        try:
            f(a.out)
        except Exception:
            failed.append(name); print(f'{name} FAILED:\n{traceback.format_exc()}', file=sys.stderr)
    print(f'wrote {a.out}; failed: {", ".join(failed) or "none"}; missing inputs: {len(set(mr.missing))}', file=sys.stderr)
    for m in sorted(set(mr.missing)):
        print(f'  missing: {m}', file=sys.stderr)


if __name__ == '__main__':
    main()
