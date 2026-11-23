#!/usr/bin/env python3
"""make_article_figures.py: multi-panel figures (Fig1-Fig7, FigS1) and placeholder numbers for the v6-vs-v5.2 article.

    python make_article_figures.py --base /scratch/stein.an/basecall_cmp --out <article>/figures [--style lab.mplstyle]

Reads only analysis outputs (no BAMs, no pileups):
  basecall_qc/<s>/tables/per_read_<run>.tsv.gz        Fig1 B-D       confusion/<s>/tables/*.tsv        Fig2 A-B
  plots/<s>_<arm>/tables/{pairwise_metrics,discordant_sites}.tsv, plots/HEK293_replicates_<run>_<arm>/...   Fig2 C-D, Fig6 B
  truthsets/<s>_<arm>_tol<t>/tables/{pairs,sets}.tsv  Fig3, FigS1    ortho/<s>_<arm>/tables/{ci,pr_curves,summary,orthogonal_qc}.tsv  Fig4-6
  motif/<s>_<arm>/tables/motif.tsv                    Fig5 B-C       tower_ubam/tower_basecall.log     Fig1 E
  sampleprobs_050/..., pileup064/<run>/<s>.thresholds.tsv, sp064_diag/nosampling/thresholds.tsv          Fig7
Writes <out>/Fig*.pdf (vector) and .png, and <out>/article_numbers.md with the values for the article's placeholders.
Each figure is built independently: a missing input is reported and that figure skipped (the article shows its placeholder).
"""
import argparse
import os
import re
import sys
import traceback

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

HEK = ['08_01_25_RNA004_HEK293', '09_09_25_RNA004_HEK293', '10_06_25_RNA004_HEK293_CP']
SAMPLES = ['GM12878'] + HEK
SLAB = {'GM12878': 'GM12878', HEK[0]: 'HEK293 08/01', HEK[1]: 'HEK293 09/09', HEK[2]: 'HEK293 10/06'}
RUNS = ['A_d212_v6', 'C_d212_v52', 'B_d130_v52']
RLAB = {'A_d212_v6': 'v6 / Dorado 2.1.2', 'C_d212_v52': 'v5.2 / Dorado 2.1.2', 'B_d130_v52': 'v5.2 / Dorado 1.3.0'}
COL = {'A_d212_v6': '#d95f02', 'C_d212_v52': '#1b9e77', 'B_d130_v52': '#7570b3'}
ARMS = {'nofilter': 'no threshold', 'lab050': 'lab rule'}
MODS = ['m6A', 'Inosine', '2OmeA', 'm5C', '2OmeC', 'Psi', '2OmeU', '2OmeG']
MLAB = {'m6A': 'm$^6$A', 'Inosine': 'Inosine', '2OmeA': 'Am', 'm5C': 'm$^5$C', '2OmeC': 'Cm', 'Psi': r'$\Psi$', '2OmeU': 'Um',
        '2OmeG': 'Gm'}
CODE = {'a': 'm6A', '17596': 'Inosine', '69426': '2OmeA', 'm': 'm5C', '19228': '2OmeC', '17802': 'Psi', '19227': '2OmeU',
        '19229': '2OmeG'}
KEYSETS = [('m6A', 'GLORI1.0∩2.0'), ('Inosine', 'SLIC-seq_2of3'), ('Psi', 'BID∩PRAISE'), ('m5C', 'UBS-seq')]
DERIVED = ('∩', '_2of3', '2+methods', 'PUS-dependent', 'NSUN2-dependent')
W = 7.2                                                   # two-column width, inches
missing, numbers = [], []


def disp(n):
    """Set names for display: the default font lacks the intersection glyph."""
    return str(n).replace('∩', r' $\cap$ ').replace('_2of3', ', 2+ reps')


def span(g):
    """Error-bar half-widths from mean/min/max, clipped at zero against rounding."""
    return [np.clip(g['mean'] - g['min'], 0, None), np.clip(g['max'] - g['mean'], 0, None)]


def read(path, **kw):
    if not os.path.exists(path):
        missing.append(path)
        return None
    return pd.read_csv(path, sep='\t', **kw)


def label(ax, s, x=-0.14, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left')


def save(fig, out, name):
    fig.savefig(f'{out}/{name}.pdf'); fig.savefig(f'{out}/{name}.png', dpi=300); plt.close(fig)
    print(f'  {name}', file=sys.stderr)


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


# ================================================================ Fig 1: design and basecalling
def fig1(B, out, s='GM12878'):
    fig = plt.figure(figsize=(W, 4.9), layout='constrained')
    gs = fig.add_gridspec(2, 4, height_ratios=[0.9, 1.1])
    ax = fig.add_subplot(gs[0, :]); label(ax, 'A', x=-0.02); ax.axis('off'); ax.set_xlim(0, 10); ax.set_ylim(0, 3)
    def box(x, y, w, h, txt, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.08', fc=fc, ec='k', lw=0.6))
        ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=6.5)
    box(0.1, 1.0, 1.9, 1.0, 'identical pod5 subsets\nGM12878 + 3 HEK293 runs\n1.5-2.1 M reads each', '#f0f0f0')
    for i, (r, txt) in enumerate([('A_d212_v6', 'A: v6.0.0 SUP + v2/v1 mod models\nDorado 2.1.2'),
                                  ('C_d212_v52', 'C: v5.2.0 SUP + v1 mod models\nDorado 2.1.2'),
                                  ('B_d130_v52', 'B: v5.2.0 SUP + v1 mod models\nDorado 1.3.0 (production)')]):
        y = 2.15 - i * 0.95
        box(2.9, y, 2.6, 0.7, txt, matplotlib.colors.to_rgba(COL[r], 0.25))
        ax.annotate('', xy=(2.88, y + 0.35), xytext=(2.02, 1.5), arrowprops=dict(arrowstyle='->', lw=0.6))
    box(6.3, 1.0, 1.6, 1.0, 'minimap2 splice\nmodkit 0.6.4 pileup\n2 calling rules', '#f0f0f0')
    for i in range(3):
        ax.annotate('', xy=(6.28, 1.5), xytext=(5.52, 2.5 - i * 0.95), arrowprops=dict(arrowstyle='->', lw=0.6))
    box(8.4, 1.0, 1.5, 1.0, 'A vs C: models\nC vs B: Dorado\nvs orthogonal maps', '#f0f0f0')
    ax.annotate('', xy=(8.38, 1.5), xytext=(7.92, 1.5), arrowprops=dict(arrowstyle='->', lw=0.6))

    t = {r: read(f'{B}/basecall_qc/{s}/tables/per_read_{r}.tsv.gz') for r in RUNS}
    ax = fig.add_subplot(gs[1, 0]); label(ax, 'B')
    bins = np.linspace(0.93, 1.0, 71)
    for r in RUNS:
        if t[r] is not None:
            ax.hist(t[r]['identity'], bins=bins, histtype='step', color=COL[r], lw=0.9, density=True, label=RLAB[r])
    ax.set_xlabel('read identity'); ax.set_ylabel('density'); ax.legend(frameon=False, fontsize=4.8, loc='upper left')
    ax.set_title(SLAB[s], fontsize=7)
    ax = fig.add_subplot(gs[1, 1]); label(ax, 'C')
    for x, y, ls in (('A_d212_v6', 'C_d212_v52', '-'), ('C_d212_v52', 'B_d130_v52', '--')):
        if t[x] is not None and t[y] is not None:
            m = t[x][['read_id', 'identity']].merge(t[y][['read_id', 'identity']], on='read_id')
            d = (m['identity_x'] - m['identity_y']) * 100
            ax.hist(d.clip(-1.5, 1.5), bins=np.linspace(-1.5, 1.5, 61), histtype='step', color=COL[x], ls=ls, lw=0.9,
                    label=f'{RLAB[x].split(" /")[0] if x == "A_d212_v6" else "Dorado 2.1.2"} - '
                          f'{RLAB[y].split(" /")[0] if x == "A_d212_v6" else "1.3.0"} (median {d.median():+.2f})')
    ax.axvline(0, color='k', lw=0.4); ax.set_xlabel('identity difference, same read (points)'); ax.set_ylabel('reads')
    ax.legend(frameon=False, fontsize=4.8, loc='upper left')
    ax = fig.add_subplot(gs[1, 2]); label(ax, 'D')
    if t['A_d212_v6'] is not None and t['C_d212_v52'] is not None:
        m = t['A_d212_v6'][['read_id', 'polya']].merge(t['C_d212_v52'][['read_id', 'polya']], on='read_id').dropna()
        ax.hexbin(m['polya_y'], m['polya_x'], gridsize=60, extent=(0, 250, 0, 250), bins='log', cmap='Greys', mincnt=1)
        ax.plot([0, 250], [0, 250], color=COL['A_d212_v6'], lw=0.6)
        ax.text(0.04, 0.94, f'r = {m["polya_x"].corr(m["polya_y"]):.3f}\nmedian diff {np.median(m["polya_x"] - m["polya_y"]):+.0f} nt',
                transform=ax.transAxes, va='top', fontsize=6)
    ax.set_xlabel('poly(A), v5.2 (nt)'); ax.set_ylabel('poly(A), v6 (nt)'); ax.set_xlim(0, 250); ax.set_ylim(0, 250)
    ax = fig.add_subplot(gs[1, 3]); label(ax, 'E')
    logf = f'{B}/tower_ubam/tower_basecall.log'
    if os.path.exists(logf):
        rt = pd.DataFrame([m_.groups() for m_ in re.finditer(r'done\s+(\S+)\s+(\S+)\s+in\s+(\d+)\s+min', open(logf).read())],
                          columns=['run', 'sample', 'minutes']).astype({'minutes': int})
        x = np.arange(len(SAMPLES)); w = 0.26
        for i, r in enumerate(RUNS):
            ax.bar(x + (i - 1) * w, rt[rt['run'] == r].set_index('sample')['minutes'].reindex(SAMPLES), w, color=COL[r])
        ax.set_xticks(x); ax.set_xticklabels([SLAB[s_].replace('HEK293 ', '') for s_ in SAMPLES], fontsize=5.5, rotation=30, ha='right')
        ax.set_ylabel('minutes (4 x A100)')
        numbers.append(('runtimes (min)', rt.pivot_table(index='sample', columns='run', values='minutes').to_string()))
    else:
        missing.append(logf)
    save(fig, out, 'Fig1')


# ================================================================ Fig 2: redistribution
def call_counts(path, mod):
    t = read(path)
    if t is None:
        return np.nan, np.nan
    return t.loc[t['call_x'] == mod, 'n'].sum(), t.loc[t['call_y'] == mod, 'n'].sum()


def fig2(B, out):
    fig, axes = plt.subplots(2, 2, figsize=(W, 5.2), layout='constrained')
    ax = axes[0, 0]; label(ax, 'A')
    rows = []
    for s in SAMPLES:
        for x, y, kind in (('A_d212_v6', 'C_d212_v52', 'models'), ('C_d212_v52', 'B_d130_v52', 'Dorado')):
            for mod in MODS:
                nx, ny = call_counts(f'{B}/confusion/{s}/tables/confusion_{x}_vs_{y}.tsv', mod)
                rows.append({'sample': s, 'kind': kind, 'mod': mod, 'ratio': nx / ny if ny else np.nan})
    r = pd.DataFrame(rows)
    xs = np.arange(len(MODS))
    for kind, col, w, off in (('models', COL['A_d212_v6'], 0.5, 0), ('Dorado', '0.3', 0.5, 0)):
        g = r[r['kind'] == kind].groupby('mod')['ratio'].agg(['mean', 'min', 'max']).reindex(MODS)
        lg = np.log2(g)
        if kind == 'models':
            ax.bar(xs, lg['mean'], w, color=col, label='v6 / v5.2 (same Dorado)')
            ax.errorbar(xs, lg['mean'], yerr=span(lg), fmt='none', ecolor='k', lw=0.6, capsize=1.5)
        else:
            ax.plot(xs, lg['mean'], '_', color=col, ms=12, mew=1.5, label='Dorado 2.1.2 / 1.3.0 (same models)')
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(xs); ax.set_xticklabels([MLAB[m] for m in MODS])
    ax.set_ylabel('log$_2$ calls per read, ratio')
    ax.legend(frameon=False, fontsize=5.5, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2)
    numbers.append(('per-read call ratio, v6/v5.2 (min-max over samples)',
                    r[r['kind'] == 'models'].groupby('mod')['ratio'].agg(['min', 'max']).reindex(MODS).round(2).to_string()))

    ax = axes[0, 1]; label(ax, 'B')
    c = read(f'{B}/confusion/GM12878/tables/agreement_by_confidence.tsv')
    if c is not None:
        c = c[c['pair'] == 'A_d212_v6 vs C_d212_v52']
        for base, col in zip('ACGT', ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00']):
            d = c[c['base'] == base]
            mid = d['call_prob_x_bin'].map(lambda b: np.mean([float(v) for v in re.findall(r'[\d.]+', b)][:2]))
            ax.plot(mid, d['identical_share'] * 100, '-o', ms=2.5, lw=0.9, color=col, label='U' if base == 'T' else base)
    ax.set_xlabel('v6 call probability'); ax.set_ylabel('same call by v5.2 (%)'); ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=6, title='reference base', title_fontsize=6); ax.set_title('GM12878', fontsize=7)

    for ax, arm, lab in ((axes[1, 0], 'nofilter', 'C'), (axes[1, 1], 'lab050', 'D')):
        label(ax, lab)
        same = []
        for s in HEK:
            t = read(f'{B}/plots/{s}_{arm}/tables/pairwise_metrics.tsv')
            if t is not None:
                same.append(t)
        reps = []
        for rr in ('A_d212_v6', 'C_d212_v52'):
            t = read(f'{B}/plots/HEK293_replicates_{rr}_{arm}/tables/pairwise_metrics.tsv')
            if t is not None:
                reps.append(t.assign(model=rr))
        if not same:
            continue
        x = pd.concat(same)
        series = {'v6 vs v5.2, same reads': (x[(x['run_x'] == 'A_d212_v6') & (x['run_y'] == 'C_d212_v52')].groupby('modification')['icc'].mean(), COL['A_d212_v6']),
                  'Dorado 2.1.2 vs 1.3.0, same reads': (x[(x['run_x'] == 'C_d212_v52') & (x['run_y'] == 'B_d130_v52')].groupby('modification')['icc'].mean(), '0.55')}
        if reps:
            g = pd.concat(reps).groupby(['model', 'modification'])['icc'].mean()
            series['v6, across HEK293 samples'] = (g.loc['A_d212_v6'] if 'A_d212_v6' in g.index.get_level_values(0) else pd.Series(dtype=float), '#fdb863')
            series['v5.2, across HEK293 samples'] = (g.loc['C_d212_v52'] if 'C_d212_v52' in g.index.get_level_values(0) else pd.Series(dtype=float), '#80cdc1')
        w = 0.8 / len(series)
        for i, (name, (v, col)) in enumerate(series.items()):
            ax.bar(xs + (i - (len(series) - 1) / 2) * w, v.reindex(MODS), w, color=col, label=name)
        ax.axhline(0, color='k', lw=0.4); ax.set_xticks(xs); ax.set_xticklabels([MLAB[m] for m in MODS]); ax.set_ylim(-0.3, 1.05)
        ax.set_ylabel('site ICC(2,1)'); ax.set_title(f'HEK293, {ARMS[arm]}', fontsize=7)
        if arm == 'nofilter':
            ax.legend(frameon=False, fontsize=5.5, loc='upper right', ncol=1)
    save(fig, out, 'Fig2')


# ================================================================ Fig 3 / S1: truth-set concordance
def matrix(ax, p, sz, mod, title):
    p = p[(p['mod'] == mod) & (p['universe'] == 'detectable')]
    names = [n for n in dict.fromkeys(list(p['set_a']) + list(p['set_b'])) if not any(k in n for k in DERIVED)]
    names = [n for n in names if not n.startswith('nanopore')] + [n for n in names if n.startswith('nanopore')]
    p = p[p['set_a'].isin(names) & p['set_b'].isin(names)]
    M = np.full((len(names), len(names)), np.nan)
    for _, r in p.iterrows():
        i, j = names.index(r['set_a']), names.index(r['set_b'])
        M[i, j], M[j, i] = r['share_a_in_b'], r['share_b_in_a']
    im = ax.imshow(M * 100, cmap='viridis', vmin=0, vmax=100)
    for i in range(len(names)):
        for j in range(len(names)):
            if i != j and np.isfinite(M[i, j]):
                ax.text(j, i, f'{100 * M[i, j]:.0f}', ha='center', va='center', fontsize=5.5, color='white' if M[i, j] < 0.6 else 'black')
    short = lambda n: n.replace('nanopore ', '').replace(' (2.1.2)', '').replace(' (1.3.0)', '/1.3.0')
    sz_ = sz[(sz['mod'] == mod) & (sz['universe'] == 'detectable')].set_index('set')['sites']
    ax.set_xticks(range(len(names))); ax.set_xticklabels([short(n) for n in names], rotation=45, ha='right', fontsize=5.5)
    ax.set_yticks(range(len(names))); ax.set_yticklabels([f'{short(n)} ({int(sz_.get(n, 0)):,})' for n in names], fontsize=5.5)
    ax.set_title(title, fontsize=7)
    return im


def fig3(B, out, s=HEK[0], arm='nofilter', tol=0, name='Fig3'):
    d = f'{B}/truthsets/{s}_{arm}_tol{tol}/tables'
    p, sz = read(f'{d}/pairs.tsv'), read(f'{d}/sets.tsv')
    if p is None or sz is None:
        return
    fig, axes = plt.subplots(2, 2, figsize=(W, 6.6), layout='constrained')
    im = None
    for ax, lab, mod in zip(axes.flat, 'ABCD', ['m6A', 'Psi', 'Inosine', 'm5C']):
        label(ax, lab, x=-0.45)
        if (p['mod'] == mod).any():
            im = matrix(ax, p, sz, mod, MLAB[mod])
        else:
            ax.axis('off')
    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.5, label='% of row set found in column set')
    save(fig, out, name)
    if tol == 0:
        for mod in ['m6A', 'Psi', 'Inosine', 'm5C']:
            x = p[(p['mod'] == mod) & (p['universe'] == 'detectable')]
            numbers.append((f'truth-set concordance, {mod}, {SLAB[s]}, {ARMS[arm]}, detectable positions, exact',
                            x[['set_a', 'set_b', 'sites_a', 'sites_b', 'share_a_in_b', 'share_b_in_a', 'jaccard', 'levels_r']]
                            .round(3).to_string(index=False)))


# ================================================================ Fig 4: accuracy
def fig4(B, out, s=HEK[0]):
    fig = plt.figure(figsize=(W, 6.4), layout='constrained')
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.25])
    for k, (mod, st) in enumerate(KEYSETS):
        ax = fig.add_subplot(gs[0, k])
        if k == 0:
            label(ax, 'A', x=-0.35)
        for arm, ls in (('nofilter', '-'), ('lab050', '--')):
            c = read(f'{B}/ortho/{s}_{arm}/tables/pr_curves.tsv')
            if c is None:
                continue
            for r in RUNS:
                g = c[(c['set'] == st) & (c['mod'] == mod) & (c['run'] == r)].sort_values('cutoff_pct')
                ax.plot(g['recall'], g['precision'], ls, color=COL[r], lw=0.9, label=f'{RLAB[r]}' if arm == 'nofilter' else None)
                p20 = g[g['cutoff_pct'] == 20]
                ax.plot(p20['recall'], p20['precision'], 'o', color=COL[r], ms=2.5, mec='k', mew=0.3)
        ax.set_title(f'{MLAB[mod]} ({disp(st)})', fontsize=6.5)
        ax.set_xlabel('recall'); ax.set_ylim(bottom=0)
        if k == 0:
            ax.set_ylabel('precision (lower bound)'); ax.legend(frameon=False, fontsize=5)
    ax = fig.add_subplot(gs[1, :]); label(ax, 'B', x=-0.06)
    rows = []
    for arm in ARMS:
        for s_ in HEK:
            ci = read(f'{B}/ortho/{s_}_{arm}/tables/ci.tsv')
            if ci is not None:
                rows.append(ci[ci['comparison'] == 'A_d212_v6 - C_d212_v52'].assign(arm=arm, sample=s_))
    if rows:
        v = pd.concat(rows)
        ylab, y = [], 0
        for mod, st in KEYSETS:
            for metric, mlab in (('precision_at_20', 'precision'), ('auprc', 'AUPRC'), ('recall_at_20', 'recall')):
                for arm, mk in (('nofilter', 'o'), ('lab050', 's')):
                    for j, s_ in enumerate(HEK):
                        q = v[(v['mod'] == mod) & (v['set'] == st) & (v['metric'] == metric) & (v['arm'] == arm) & (v['sample'] == s_)]
                        if q.empty:
                            continue
                        q = q.iloc[0]
                        ax.errorbar(q['estimate'], y + j * 0.12 + (0.4 if arm == 'lab050' else 0), xerr=[[max(q['estimate'] - q['ci_low'], 0)], [max(q['ci_high'] - q['estimate'], 0)]],
                                    fmt=mk, ms=2.5, color=COL['A_d212_v6'] if arm == 'nofilter' else '0.35', lw=0.6, capsize=0)
                ylab.append((y + 0.32, f'{MLAB[mod]} {mlab}')); y += 1
        ax.axvline(0, color='k', lw=0.4); ax.set_yticks([t for t, _ in ylab]); ax.set_yticklabels([l for _, l in ylab], fontsize=6)
        ax.invert_yaxis(); ax.set_xlabel('v6 minus v5.2 (same Dorado), 95% paired block-bootstrap interval')
        ax.plot([], [], 'o', color=COL['A_d212_v6'], ms=3, label='no threshold'); ax.plot([], [], 's', color='0.35', ms=3, label='lab rule')
        ax.legend(frameon=False, fontsize=6, loc='lower right', title='three HEK293 samples per rule', title_fontsize=6)
    save(fig, out, 'Fig4')


# ================================================================ Fig 5: m6A
def fig5(B, out):
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.9), layout='constrained', gridspec_kw={'width_ratios': [1, 1.2, 1.5]})
    ax = axes[0]; label(ax, 'A', x=-0.3)
    rows = []
    for arm in ARMS:
        for s in HEK:
            c = read(f'{B}/ortho/{s}_{arm}/tables/pr_curves.tsv')
            if c is None:
                continue
            g = c[(c['set'] == 'GLORI1.0∩2.0') & (c['mod'] == 'm6A')]
            ref = g[(g['run'] == 'C_d212_v52') & (g['cutoff_pct'] == 20)]
            v6 = g[g['run'] == 'A_d212_v6'].sort_values('cutoff_pct')
            if ref.empty or v6.empty:
                continue
            ok = v6[v6['precision'] >= ref['precision'].iloc[0]]
            rows.append({'arm': arm, 'sample': s, 'v52_recall_20': ref['recall'].iloc[0],
                         'v6_cutoff': ok['cutoff_pct'].iloc[0] if len(ok) else np.nan, 'v6_recall_matched': ok['recall'].iloc[0] if len(ok) else np.nan})
    if rows:
        op = pd.DataFrame(rows); x = np.arange(len(op))
        ax.bar(x - 0.2, op['v52_recall_20'], 0.4, color=COL['C_d212_v52'], label='v5.2 at 20%')
        ax.bar(x + 0.2, op['v6_recall_matched'], 0.4, color=COL['A_d212_v6'], label='v6 at matched precision')
        for xi, c in zip(x, op['v6_cutoff']):
            ax.text(xi + 0.2, 0.02, f'{c:.0f}%', rotation=90, ha='center', va='bottom', fontsize=5, color='white')
        ax.set_xticks(x); ax.set_xticklabels([f'{SLAB[s][7:12]}\n{"none" if a == "nofilter" else "lab"}' for a, s in zip(op['arm'], op['sample'])], fontsize=5)
        ax.set_ylabel(r'recall, GLORI 1.0 $\cap$ 2.0'); ax.set_ylim(0, 1); ax.legend(frameon=False, fontsize=5, loc='upper center', bbox_to_anchor=(0.5, 1.2))
        numbers.append(('m6A operating points', op.round(3).to_string(index=False)))
    mot = []
    for arm in ARMS:
        for s in HEK:
            t = read(f'{B}/motif/{s}_{arm}/tables/motif.tsv')
            if t is not None:
                mot.append(t.assign(arm=arm))
    if mot:
        m = pd.concat(mot)
        classes = ['both', 'v6 only', 'v5.2 only', 'background (all covered A)']
        clab = ['both', 'v6 only', 'v5.2 only', 'all A']
        for ax, lab, cols, title in ((axes[1], 'B', [('drach_share', 'DRACH'), ('glori_share', r'GLORI 1.0 $\cap$ 2.0')], 'context and GLORI support'),
                                     (axes[2], 'C', [(c, c.replace('share_', '')) for c in m.columns if c.startswith('share_') and 'GLORI' not in c],
                                      'support by other chemistries')):
            label(ax, lab, x=-0.2)
            cols = [(c, l) for c, l in cols if c in m.columns and m[c].notna().any()]
            if not cols:
                ax.text(0.5, 0.5, 'pending', ha='center', transform=ax.transAxes); continue
            k = len(cols) * 2; w = 0.8 / k; xx = np.arange(len(classes)); j = 0
            pal = ['#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#fdbf6f', '#e31a1c', '#cab2d6']
            for ci_, (c, l) in enumerate(cols):
                for arm, alpha in (('nofilter', 1.0), ('lab050', 0.45)):
                    g = m[m['arm'] == arm].groupby('class')[c].agg(['mean', 'min', 'max']).reindex(classes) * 100
                    ax.bar(xx + (j - (k - 1) / 2) * w, g['mean'], w, color=pal[ci_ % len(pal)], alpha=alpha,
                           label=f'{l}' if arm == 'nofilter' else None)
                    ax.errorbar(xx + (j - (k - 1) / 2) * w, g['mean'], yerr=span(g), fmt='none', ecolor='k', lw=0.4)
                    j += 1
            ax.set_xticks(xx); ax.set_xticklabels(clab, fontsize=6); ax.set_ylabel('% of sites'); ax.set_title(title, fontsize=7)
            ax.legend(frameon=False, fontsize=5, title='solid: no threshold; light: lab rule', title_fontsize=5)
        numbers.append(('m6A class support (mean over HEK293 samples, %)',
                        (m.groupby(['arm', 'class'])[[c for c in m.columns if c.endswith('_share') or c.startswith('share_')]].mean() * 100).round(1).to_string()))
    save(fig, out, 'Fig5')


# ================================================================ Fig 6: inosine
def fig6(B, out):
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.9), layout='constrained', gridspec_kw={'width_ratios': [1.1, 1.4]})
    ax = axes[0]; label(ax, 'A', x=-0.18)
    rows = []
    for arm in ARMS:
        for s in HEK:
            ci = read(f'{B}/ortho/{s}_{arm}/tables/ci.tsv')
            if ci is not None:
                rows.append(ci[(ci['set'] == 'SLIC-seq_2of3') & ci['comparison'].isin(RUNS)].assign(arm=arm, sample=s))
    if rows:
        v = pd.concat(rows); xt = []
        for k, (metric, mlab) in enumerate((('precision_at_20', 'precision'), ('recall_at_20', 'recall'), ('auprc', 'AUPRC'))):
            for a_i, arm in enumerate(ARMS):
                for r_i, r in enumerate(RUNS):
                    q = v[(v['metric'] == metric) & (v['arm'] == arm) & (v['comparison'] == r)]
                    xpos = k * 2.6 + a_i * 1.2 + r_i * 0.3
                    ax.errorbar(np.full(len(q), xpos), q['estimate'], yerr=[np.clip(q['estimate'] - q['ci_low'], 0, None), np.clip(q['ci_high'] - q['estimate'], 0, None)],
                                fmt='o' if arm == 'nofilter' else 's', ms=2.3, color=COL[r], lw=0.5, label=RLAB[r] if k == 0 and a_i == 0 else None)
                xt.append((k * 2.6 + a_i * 1.2 + 0.3, f'{mlab}\n{"none" if arm == "nofilter" else "lab"}'))
        ax.set_xticks([t for t, _ in xt]); ax.set_xticklabels([l for _, l in xt], fontsize=5.5); ax.set_ylim(bottom=0)
        ax.set_ylabel('value (three HEK293 samples)'); ax.legend(frameon=False, fontsize=5); ax.set_title('Inosine vs SLIC-seq (2+ replicates)', fontsize=7)
    ax = axes[1]; label(ax, 'B', x=-0.12)
    d = read(f'{B}/plots/GM12878_nofilter/tables/discordant_sites.tsv')
    if d is not None:
        d = d[(d['modification'] == 'Inosine') & (d['run_x'] == 'A_d212_v6') & (d['run_y'] == 'C_d212_v52') & (d['chrom'] == 'chr1')
              & d['pos_1based'].between(25_900_500, 25_904_500)].sort_values('pos_1based')
        if len(d):
            ax.vlines(d['pos_1based'], 0, d['pct_y'], color=COL['C_d212_v52'], lw=1.2, label='v5.2 occupancy')
            ax.plot(d['pos_1based'], d['pct_x'], 'o', ms=2.5, color=COL['A_d212_v6'], label='v6 occupancy')
            ax.set_ylim(0, 105); ax.set_ylabel('inosine occupancy (%)'); ax.set_xlabel('chr1 position (STMN1 exons, minus strand)')
            a2 = ax.twinx(); a2.plot(d['pos_1based'], d['cov_x'], '-', color=COL['A_d212_v6'], lw=0.5, alpha=0.6)
            a2.plot(d['pos_1based'], d['cov_y'], '-', color=COL['C_d212_v52'], lw=0.5, alpha=0.6); a2.set_ylabel('valid reads (lines)', fontsize=6)
            a2.tick_params(labelsize=6); a2.spines['top'].set_visible(False)
            ax.legend(frameon=False, fontsize=5.5, loc='lower left', bbox_to_anchor=(0, 1.08), ncol=2); ax.ticklabel_format(axis='x', style='plain', useOffset=False)
            ax.set_title(f'GM12878, {len(d)} positions with |difference| $\\geq$ 50 points', fontsize=7)
            numbers.append(('STMN1 cluster', f'{len(d)} positions; v5.2 pct {d["pct_y"].min():.0f}-{d["pct_y"].max():.0f}; v6 pct max {d["pct_x"].max():.1f}; '
                                             f'valid reads v6 {d["cov_x"].median():.0f}, v5.2 {d["cov_y"].min()}-{d["cov_y"].max()}'))
    save(fig, out, 'Fig6')


# ================================================================ Fig 7: thresholds
def fig7(B, out):
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.6), layout='constrained', gridspec_kw={'width_ratios': [1.6, 1]})
    ax = axes[0]; label(ax, 'A', x=-0.1)
    bases = ['A', 'C', 'G', 'U']
    for i, b in enumerate(bases):
        for j, r in enumerate(RUNS):
            vals = [thresholds(f'{B}/sampleprobs_050/{r}/thresholds.tsv' if s == 'GM12878' else f'{B}/sampleprobs_050/{r}/{s}/thresholds.tsv').get(b, np.nan)
                    for s in SAMPLES]
            ax.plot(np.full(len(vals), i + (j - 1) * 0.22), vals, 'o', ms=3, color=COL[r], label=RLAB[r] if i == 0 else None, alpha=0.85)
    ax.set_xticks(range(4)); ax.set_xticklabels(bases); ax.set_xlabel('canonical base'); ax.set_ylabel('10th-percentile threshold')
    ax.set_ylim(0.5, 1.01); ax.legend(frameon=False, fontsize=5.5, loc='lower right'); ax.set_title('modkit 0.5.0, four samples', fontsize=7)
    ax = axes[1]; label(ax, 'B', x=-0.18)
    est = {'0.5.0': thresholds(f'{B}/sampleprobs_050/B_d130_v52/thresholds.tsv'),
           '0.6.4, reference positions': thresholds(f'{B}/pileup064/B_d130_v52/GM12878.thresholds.tsv'),
           '0.6.4, all calls': thresholds(f'{B}/sp064_diag/nosampling/thresholds.tsv')}
    w = 0.26
    for k, (name, t) in enumerate(est.items()):
        ax.bar(np.arange(4) + (k - 1) * w, [t.get(b, np.nan) for b in bases], w, color=['0.25', '#e31a1c', '0.65'][k], label=name)
    ax.set_xticks(range(4)); ax.set_xticklabels(bases); ax.set_ylim(0.5, 1.13); ax.set_ylabel('threshold')
    ax.legend(frameon=False, fontsize=5, loc='upper left', ncol=3); ax.set_title('v5.2 / Dorado 1.3.0, GM12878', fontsize=7)
    save(fig, out, 'Fig7')


# ================================================================ numbers for the placeholders
def numbers_new_sets(B):
    qc = read(f'{B}/ortho/{HEK[0]}_nofilter/tables/orthogonal_qc.tsv')
    if qc is not None:
        numbers.append(('truth sets after the reference check', qc.drop(columns=[c for c in qc.columns if c.startswith('match_') and c != 'match_offset_0'])
                        .round(4).to_string(index=False)))
    new = ['m6A-SAC-seq', 'miCLIP-CIMS', 'miCLIP-CITS', 'DART-seq', 'm6A_2+methods', 'PRAISE_PUS-dependent', 'UBS-seq_NSUN2-dependent']
    rows = []
    for arm in ARMS:
        for s in HEK:
            sm = read(f'{B}/ortho/{s}_{arm}/tables/summary.tsv')
            ci = read(f'{B}/ortho/{s}_{arm}/tables/ci.tsv')
            if sm is not None:
                rows.append(sm[sm['set'].isin(new)][['set', 'run', 'sites_covered_all_runs', 'recall_at_20', 'precision_at_20', 'auprc']].assign(arm=arm, sample=SLAB[s]))
            if ci is not None:
                d = ci[ci['set'].isin(new) & (ci['comparison'] == 'A_d212_v6 - C_d212_v52')].copy()
                d['call'] = np.select([d['ci_low'] > 0, d['ci_high'] < 0], ['v6 higher', 'v6 lower'], 'no clear difference')
                numbers.append((f'v6 - v5.2 on new sets, {SLAB[s]}, {ARMS[arm]}', d[['set', 'metric', 'estimate', 'ci_low', 'ci_high', 'call']].round(4).to_string(index=False)))
    if rows:
        numbers.append(('accuracy on the new sets', pd.concat(rows).round(4).to_string(index=False)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--base', default='/scratch/stein.an/basecall_cmp')
    ap.add_argument('--out', required=True)
    ap.add_argument('--style', help='matplotlib style sheet applied before the article defaults')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.style:
        plt.style.use(a.style)
    plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                         'legend.fontsize': 6, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
    jobs = [('Fig1', lambda: fig1(a.base, a.out)), ('Fig2', lambda: fig2(a.base, a.out)),
            ('Fig3', lambda: fig3(a.base, a.out)), ('FigS1', lambda: fig3(a.base, a.out, tol=1, name='FigS1')),
            ('Fig4', lambda: fig4(a.base, a.out)), ('Fig5', lambda: fig5(a.base, a.out)), ('Fig6', lambda: fig6(a.base, a.out)),
            ('Fig7', lambda: fig7(a.base, a.out)), ('numbers', lambda: numbers_new_sets(a.base))]
    failed = []
    for name, f in jobs:
        try:
            f()
        except Exception:
            failed.append(name); print(f'{name} FAILED:\n{traceback.format_exc()}', file=sys.stderr)
    with open(f'{a.out}/article_numbers.md', 'w') as fh:
        fh.write('# Numbers for the article placeholders\n\n')
        for title, body in numbers:
            fh.write(f'## {title}\n\n```\n{body}\n```\n\n')
        fh.write('## Missing inputs\n\n' + ('\n'.join(f'- {m}' for m in sorted(set(missing))) or 'none') + '\n')
    print(f'wrote {a.out}; failed: {", ".join(failed) or "none"}; missing inputs: {len(set(missing))}', file=sys.stderr)


if __name__ == '__main__':
    main()
