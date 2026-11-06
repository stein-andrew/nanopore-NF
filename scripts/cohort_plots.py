#!/usr/bin/env python3
"""Cohort figures for RNA_analysis.nf, in the style of Ghohabi Esfahani, Stein, Akeson et al.
(Genome Biology 2026, 27:190): multi-panel figures with bold panel letters, the lab style
sheet, the paper's colours, log-scaled 2D histograms with n / r / RMSE, and dashed medians.
One PDF per figure, editable in Illustrator (TrueType text, vector shapes only).

    cohort_plots.py --results cohort-v1/results --samplesheet analysis_samplesheet.csv --out plots/
    cohort_plots.py ... --compare-batch ver004          (figures 3 and 4C; default: the last batch)
    cohort_plots.py ... --group condition [--levels control case]      (adds figure 6)

  Fig1_cohort_qc.pdf          A read composition  B non-globin reads  C genes detected (>= 10 reads)
                              D modification sites vs non-globin depth  E variant calls by class
                              F A>G/T>C share of novel calls (RNA editing among novel calls)
  Fig2_modification_calls.pdf A, B sites per occupancy bin (main mods; 2'-O-Me)  C sites per sample
                              D samples calling each site  E modified sites per gene
                              F gene-level co-occurrence, P(observed | given)
  Fig3_concordance.pdf        A-H one batch vs the others, mean occupancy per site, per mod
                              I each sample vs the rest of the cohort, against depth
  Fig4_polya.pdf              A nuclear (no globin), globin and mitochondrial tails  B medians per
                              sample  C gene median poly(A), one batch vs the others
  Fig5_expression_isoforms.pdf A sample correlation matrix  B PCA  C isoform structural categories
  Fig6_group_comparison.pdf   (with --group) volcanoes for expression, modification sites and gene
                              poly(A), heatmap of the top genes; a TSV of every test beside it

Modification sites are the cohort's final-stage sites (>= 20 reads, >= 20% after gIVT
correction); a sample calls a site where its corrected percent there is >= 20. Empty cells in
the matrices mean not measured, so comparisons use only what both sides measured. Welch's
t-test with Benjamini-Hochberg FDR in figure 6 is exploratory; use DESeq2 or edgeR on the
counts for inference on expression. polars reads the tables (as text); numpy and scipy compute.
"""
import argparse
import gzip
import logging
import math
import os
import sys
import warnings

import numpy as np
import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib import font_manager                              # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm   # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter         # noqa: E402
from scipy import stats                                          # noqa: E402

LAB_STYLE = '/projects/Genometechlab/matplotlib_stylesheets/styles/genometechlab_main.mplstyle'
FONT_FALLBACKS = ['Helvetica', 'Arial', 'Liberation Sans', 'Nimbus Sans', 'TeX Gyre Heros']
GLOBIN = ['HBA1', 'HBA2', 'HBB', 'HBD', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ', 'HBE1']

# the paper's colours (sampled from its figures); 2'-O-Me: darker shades of the base's colour
MODS = [('a', 'm$^6$A', '#398DC3'), ('m', 'm$^5$C', '#D593BD'), ('17802', 'Ψ', '#DB8151'),
        ('17596', 'Ino', '#4BB193'), ('69426', 'Am', '#1F5F8B'), ('19228', 'Cm', '#9C5B8E'),
        ('19227', 'Um', '#9E4F25'), ('19229', 'Gm', '#8A7A3B')]
MOD_NAME = {c: n for c, n, _ in MODS}
MOD_COLOUR = {c: k for c, _, k in MODS}
MAIN_MODS = ['a', 'm', '17802', '17596']
OME_MODS = ['69426', '19228', '19227', '19229']
BATCH_COLOURS = ['#5769AB', '#EDB157', '#935193', '#6FBD81', '#5D9FC9', '#E76969', '#AB7B4B', '#63CFE1']
COMPARTMENT = [('nuclear', 'nuclear (no globin)', '#5769AB'), ('globin', 'globin', '#BD7B7B'),
               ('mito', 'mitochondrial', '#EDB157')]
MATRIX_CMAP = LinearSegmentedColormap.from_list('paper_matrix', ['#3B0000', '#7B3333', '#C9A9A2', '#FBF3EE'])
PANEL = (3.25, 2.25)                                               # the style sheet's figure size
OCC_BINS = list(range(20, 100, 10))                                # 20-30 ... 90-100


# ---------------------------------------------------------------- tables
def read_table(path, sep='\t', columns=None):
    data = gzip.open(path).read() if path.endswith('.gz') else open(path, 'rb').read()
    df = pl.read_csv(data, separator=sep, infer_schema_length=0, columns=columns)
    if df.columns and df.columns[0].startswith('#'):               # '#chrom' headers of the BEDs
        df = df.rename({df.columns[0]: df.columns[0].lstrip('#')})
    return df


def numbers(df, cols):
    return df.select([pl.col(c).cast(pl.Float64, strict=False) for c in cols]).to_numpy().astype(float)


def text(df, col):
    return df[col].to_list()


def write_rows(path, header, rows):
    with open(path, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for r in rows:
            out.write('\t'.join('NA' if (isinstance(v, float) and math.isnan(v)) else
                                (f'{v:.6g}' if isinstance(v, float) else str(v)) for v in r) + '\n')


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        ps = p[ok]
        order = np.argsort(ps)
        ranked = np.minimum.accumulate((ps[order] * len(ps) / np.arange(1, len(ps) + 1))[::-1])[::-1]
        out = np.empty_like(ps)
        out[order] = np.clip(ranked, 0, 1)
        q[ok] = out
    return q


def welch(x, ga, gb, min_n=2):
    """Welch's t-test per row, B minus A, ignoring missing values; rows with fewer than min_n
    values on either side get NaN. Returns (difference of means, p)."""
    a, b = x[:, ga], x[:, gb]
    na, nb = (~np.isnan(a)).sum(1), (~np.isnan(b)).sum(1)
    with np.errstate(invalid='ignore', divide='ignore'):
        ma, mb = np.nanmean(a, 1), np.nanmean(b, 1)
        va, vb = np.nanvar(a, 1, ddof=1), np.nanvar(b, 1, ddof=1)
        se2 = va / na + vb / nb
        t = (mb - ma) / np.sqrt(se2)
        df = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
        p = 2 * stats.t.sf(np.abs(t), df)
    p[(na < min_n) | (nb < min_n) | ~(se2 > 0)] = np.nan
    return mb - ma, p


# ---------------------------------------------------------------- style
def apply_style(path):
    notes = []
    if path and os.path.exists(path):
        plt.style.use(path)
        notes.append(f'style: {path}')
    elif path:
        notes.append(f'style {path} not found: matplotlib defaults')
    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 100,
                         'figure.constrained_layout.use': True, 'axes.titleweight': 'normal',
                         'axes.titlesize': 13, 'legend.frameon': False})
    installed = {f.name for f in font_manager.fontManager.ttflist}
    wanted = (plt.rcParams['font.sans-serif'] or ['Helvetica'])[0]
    chosen = next((f for f in [wanted] + FONT_FALLBACKS if f in installed), 'DejaVu Sans')
    plt.rcParams['font.sans-serif'] = [chosen, 'DejaVu Sans']
    plt.rcParams['mathtext.fontset'] = 'custom'
    for k, v in (('mathtext.rm', chosen), ('mathtext.sf', chosen), ('mathtext.it', f'{chosen}:italic'),
                 ('mathtext.bf', f'{chosen}:bold')):
        plt.rcParams[k] = v
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    notes.append(f'font: {chosen}' + ('' if chosen == wanted else f' ({wanted} is not installed)'))
    return notes


def figure(rows, cols, widths=None, heights=None):
    w = sum(widths) if widths else cols * PANEL[0]
    h = sum(heights) if heights else rows * PANEL[1]
    fig = plt.figure(figsize=(w + 0.4 * cols, h + 0.5 * rows))
    gs = fig.add_gridspec(rows, cols, width_ratios=widths, height_ratios=heights)
    return fig, gs


def letter(ax, s):
    ax.text(-0.02, 1.02, s, transform=ax.transAxes, fontsize=18, fontweight='bold', ha='right', va='bottom')


def stat_box(ax, lines, loc='upper left'):
    x, ha = (0.04, 'left') if 'left' in loc else (0.96, 'right')
    ax.text(x, 0.96, '\n'.join(lines), transform=ax.transAxes, ha=ha, va='top', fontsize=10,
            bbox=dict(boxstyle='square,pad=0.3', fc='white', ec='black', lw=0.6))


def sample_axis(ax, samples, batch, colours, fontsize=9):
    ax.set_xticks(range(len(samples)), samples, rotation=90, fontsize=fontsize)
    for lab, s in zip(ax.get_xticklabels(), samples):
        lab.set_color(colours[batch[s]])


def batch_legend(ax, colours, fontsize=9, **kw):
    handles = [plt.Line2D([], [], marker='o', ls='', color=c, label=b) for b, c in colours.items()]
    ax.legend(handles=handles, title='batch', fontsize=fontsize, title_fontsize=fontsize + 1, **kw)


def fig_batch_legend(fig, colours):
    handles = [plt.Line2D([], [], marker='o', ls='', color=c, label=b) for b, c in colours.items()]
    fig.legend(handles=handles, title=None, loc='lower center', bbox_to_anchor=(0.5, 1.0),
               ncol=len(handles), fontsize=11, handletextpad=0.3, columnspacing=1.2)


def log_axis(axis):
    """Labels at powers of ten only, as in the paper; no minor-tick labels."""
    axis.set_major_locator(LogLocator(base=10))
    axis.set_minor_formatter(NullFormatter())


def hist2d(ax, x, y, label_x, label_y, bins=40, rng=(0, 100), cmap='viridis', vmax=None):
    """The paper's 2D histograms: log colour scale, y = x dashed, n / r / RMSE boxed."""
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if len(x) < 3:
        ax.text(0.5, 0.5, 'too few sites', transform=ax.transAxes, ha='center')
        return None
    h, xe, ye = np.histogram2d(x, y, bins=bins, range=[rng, rng])
    h = np.ma.masked_equal(h.T, 0)
    m = ax.pcolormesh(xe, ye, h, cmap=cmap, norm=LogNorm(vmin=1, vmax=vmax or max(10, h.max())),
                      edgecolors='face', linewidth=0.2)
    ax.plot(rng, rng, ls='--', color='black', lw=0.8)
    ax.set_xlim(rng)
    ax.set_ylim(rng)
    ax.set_aspect('equal')
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    r = np.corrcoef(x, y)[0, 1]
    rmse = float(np.sqrt(np.mean((x - y) ** 2)))
    stat_box(ax, [f'n = {len(x):,}', f'r = {r:.3f}', f'RMSE = {rmse:.1f}'])
    return m


def colourbar(fig, mappable, ax, label='sites'):
    cb = fig.colorbar(mappable, ax=ax, fraction=0.05, pad=0.02, label=label)
    if getattr(cb, 'solids', None) is not None:
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor('face')
    return cb


# ---------------------------------------------------------------- data
class Cohort:
    def __init__(self, results, ref, sheet_path, batch_col):
        self.R, self.ref = results, ref
        self.cohort = os.path.join(results, 'cohort', ref)
        sheet = read_table(sheet_path, ',')
        ids = text(sheet, 'sample_id')
        b = text(sheet, batch_col) if batch_col in sheet.columns else ['all'] * len(ids)
        self.sheet = sheet
        self.batch = dict(zip(ids, b))
        self.samples = sorted(ids, key=lambda s: (self.batch[s], s))
        names = sorted(set(self.batch.values()))
        self.colours = {n: BATCH_COLOURS[i % len(BATCH_COLOURS)] for i, n in enumerate(names)}

    def path(self, *parts):
        return os.path.join(self.R, *parts)

    def mods_matrix(self):
        pm = read_table(os.path.join(self.cohort, f'cohort.{self.ref}.mods.pct_fpcorr.tsv.gz'))
        ms = [s for s in self.samples if s in pm.columns]
        return pm, ms, numbers(pm, ms), np.array(text(pm, 'mod'))


# ---------------------------------------------------------------- figure 1
def fig1_qc(c, notes):
    qc = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.read_qc.tsv'))
    row = {s: i for i, s in enumerate(text(qc, 'sample'))}
    ss = [s for s in c.samples if s in row]
    cols = {k: numbers(qc, [k]).ravel() for k in ('mapped_reads', 'globin_reads', 'mito_reads', 'non_globin_reads')}
    get = {k: np.array([v[row[s]] for s in ss]) for k, v in cols.items()}
    other = get['non_globin_reads'] - get['mito_reads']
    fig, gs = figure(2, 3)
    x = np.arange(len(ss))

    ax = fig.add_subplot(gs[0, 0])                                  # A composition, %
    tot = get['mapped_reads']
    parts = [('other', other, '#5769AB'), ('mitochondrial', get['mito_reads'], '#EDB157'), ('globin', get['globin_reads'], '0.78')]
    bottom = np.zeros(len(ss))
    for name, v, col in parts:
        ax.bar(x, 100 * v / tot, bottom=bottom, color=col, label=name, width=0.8)
        bottom += 100 * v / tot
    for i, t in enumerate(tot):
        ax.text(i, 101, f'{t / 1e6:.0f}', ha='center', va='bottom', fontsize=7)
    ax.set_ylim(0, 112)
    ax.set_ylabel('mapped reads (%)')
    sample_axis(ax, ss, c.batch, c.colours)
    ax.legend(fontsize=8, loc='lower left', bbox_to_anchor=(0, 1.08), ncol=3, handlelength=1, columnspacing=0.8)
    ax.text(1.0, 1.0, 'M reads', transform=ax.transAxes, ha='right', va='bottom', fontsize=7)
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B non-globin reads
    ax.bar(x, other / 1e6, color=[c.colours[c.batch[s]] for s in ss], width=0.8)
    med = np.median(other / 1e6)
    ax.axhline(med, ls='--', color='black', lw=0.8)
    ax.text(len(ss) - 0.5, med, f' {med:.2f}', va='bottom', ha='right', fontsize=9)
    ax.set_ylabel('non-globin reads (M)')
    sample_axis(ax, ss, c.batch, c.colours)
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C genes detected
    gc = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.gene_counts.tsv.gz'))
    gs_ = [s for s in ss if s in gc.columns]
    detected = (np.nan_to_num(numbers(gc, gs_)) >= 10).sum(axis=0)
    ax.scatter(np.array([depth for s, depth in zip(ss, other) if s in gs_]) / 1e6, detected,
               c=[c.colours[c.batch[s]] for s in gs_], s=30, zorder=3)
    for s, xx, yy in zip(gs_, [o for s2, o in zip(ss, other) if s2 in gs_], detected):
        ax.annotate(s[-2:], (xx / 1e6, yy), fontsize=7, xytext=(3, 2), textcoords='offset points')
    ax.set_xlabel('non-globin reads (M)')
    ax.set_ylabel('genes with ≥ 10 reads')
    letter(ax, 'C')

    ax = fig.add_subplot(gs[1, 0])                                  # D sites vs depth
    ms = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.mods.summary.tsv'))
    sites = dict(zip(text(ms, 'sample'), numbers(ms, ['final_sites']).ravel()))
    xs, ys = other / 1e6, np.array([sites.get(s, np.nan) for s in ss]) / 1e3
    ok = ~np.isnan(ys)
    ax.scatter(xs[ok], ys[ok], c=[c.colours[c.batch[s]] for s, k in zip(ss, ok) if k], s=30, zorder=3)
    for s, a, b in zip(ss, xs, ys):
        if not np.isnan(b):
            ax.annotate(s[-2:], (a, b), fontsize=7, xytext=(3, 2), textcoords='offset points')
    if ok.sum() >= 3:
        k, b0 = np.polyfit(xs[ok], ys[ok], 1)
        xx = np.linspace(0, xs[ok].max() * 1.05, 10)
        ax.plot(xx, k * xx + b0, ls='--', color='0.5', lw=0.8)
        stat_box(ax, [f'n = {int(ok.sum())}', f'r = {np.corrcoef(xs[ok], ys[ok])[0, 1]:.3f}'])
    ax.set_xlabel('non-globin reads (M)')
    ax.set_ylabel('modification sites (k)')
    letter(ax, 'D')

    summ = {}
    for s in ss:
        p = c.path('variants', c.ref, s, f'{s}.variants.summary.tsv')
        if os.path.exists(p):
            t = read_table(p)
            summ[s] = {k: (v1, v2) for k, v1, v2 in zip(text(t, 'class'), numbers(t, ['calls']).ravel(),
                                                       numbers(t, ['ag_tc_share']).ravel())}
    ax = fig.add_subplot(gs[1, 1])                                  # E variant classes
    if summ:
        bottom = np.zeros(len(ss))
        for cls, col in (('known', '#5769AB'), ('novel', '#E76969'), ('editing', '#4BB193')):
            v = np.array([summ.get(s, {}).get(cls, (0, 0))[0] for s in ss]) / 1e3
            ax.bar(x, v, bottom=bottom, color=col, label=cls, width=0.8)
            bottom += v
        masked = np.array([sum(n for k, (n, _) in summ.get(s, {}).items() if k.startswith('masked:')) for s in ss]) / 1e3
        ax.bar(x, masked, bottom=bottom, color='0.78', label='IG/TR/HLA', width=0.8)
        ax.set_ylabel('PASS variant calls (k)')
        sample_axis(ax, ss, c.batch, c.colours)
        ax.legend(fontsize=8, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=4, handlelength=1, columnspacing=0.6)
    letter(ax, 'E')

    ax = fig.add_subplot(gs[1, 2])                                  # F A>G share of novel calls
    if summ:
        for cls, col, mk in (('known', '#5769AB', 'o'), ('novel', '#E76969', 's')):
            v = np.array([summ.get(s, {}).get(cls, (0, np.nan))[1] for s in ss])
            ax.scatter(x, 100 * v, color=col, marker=mk, s=28, label=cls, zorder=3)
        ax.set_ylabel('A>G / T>C of SNVs (%)')
        ax.set_ylim(bottom=0)
        sample_axis(ax, ss, c.batch, c.colours)
        ax.legend(fontsize=8, loc='upper right')
    letter(ax, 'F')
    fig_batch_legend(fig, c.colours)
    return fig


# ---------------------------------------------------------------- figure 2
def per_sample_calls(pct, mods, samples):
    """called[mod] -> boolean matrix (sites of that mod x samples): corrected percent >= 20."""
    return {code: (pct[mods == code] >= 20) for code in MOD_NAME if (mods == code).any()}


def fig2_mods(c, notes):
    pm, ms, pct, mods = c.mods_matrix()
    called = per_sample_calls(pct, mods, ms)
    fig, gs = figure(2, 3, widths=[PANEL[0] * 1.15, PANEL[0] * 1.15, PANEL[0] * 1.3])
    bin_labels = [f'{b}' for b in OCC_BINS]

    for k, (group, title) in enumerate(((MAIN_MODS, 'A'), (OME_MODS, 'B'))):
        ax = fig.add_subplot(gs[0, k])
        present = [m for m in group if m in called]
        w = 0.8 / max(1, len(present))
        for j, code in enumerate(present):
            v = pct[mods == code]
            counts = np.array([[np.sum((v[:, si] >= b) & (v[:, si] < b + 10 if b < 90 else v[:, si] <= 100))
                                for si in range(len(ms))] for b in OCC_BINS])        # bins x samples
            med, lo, hi = np.median(counts, 1), counts.min(1), counts.max(1)
            xpos = np.arange(len(OCC_BINS)) - 0.4 + w * (j + 0.5)
            ax.bar(xpos, np.maximum(med, 0.8), width=w, color=MOD_COLOUR[code], label=MOD_NAME[code])
            ax.vlines(xpos, np.maximum(lo, 0.8), np.maximum(hi, 0.8), color='black', lw=0.6)
        ax.set_yscale('log')
        log_axis(ax.yaxis)
        ax.set_xticks(range(len(OCC_BINS)), bin_labels)
        ax.set_xlabel('modification occupancy (%, bin start)')
        ax.set_ylabel('sites per sample')
        ax.legend(fontsize=9, ncol=2, loc='upper right', handlelength=1)
        letter(ax, title)

    ax = fig.add_subplot(gs[0, 2])                                  # C sites per sample
    for code in [m for m in MAIN_MODS + OME_MODS if m in called]:
        ax.plot(range(len(ms)), called[code].sum(0), marker='o', ms=4, lw=1, color=MOD_COLOUR[code],
                ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    ax.set_ylabel('sites (≥ 20 reads, ≥ 20%)')
    sample_axis(ax, ms, c.batch, c.colours)
    ax.set_title("colours as in A, B; 2'-O-Me dashed", fontsize=9)
    letter(ax, 'C')

    ax = fig.add_subplot(gs[1, 0])                                  # D samples calling each site
    k = np.arange(1, len(ms) + 1)
    for code in [m for m in MAIN_MODS + OME_MODS if m in called]:
        n = called[code].sum(1)
        n = n[n > 0]
        frac = np.array([(n == i).mean() for i in k]) if len(n) else np.zeros(len(k))
        ax.plot(k, 100 * frac, marker='o', ms=3.5, lw=1, color=MOD_COLOUR[code],
                ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_xlabel('samples calling the site')
    ax.set_ylabel('sites (%)')
    ax.set_xticks(k[::2] if len(k) > 8 else k)
    letter(ax, 'D')

    genes = gene_sites(c)
    ax = fig.add_subplot(gs[1, 1])                                  # E sites per gene
    cap = 50
    present = [m for m in MAIN_MODS if genes.get(m)]
    if present:
        width = (cap + 2) + 10                                     # a gap between the blocks
        for j, code in enumerate(present):
            per_gene = np.array([len(v) for v in genes[code].values()])
            capped = np.minimum(per_gene, cap + 1)
            h = np.bincount(capped, minlength=cap + 2)[1:]
            xs = j * width + np.arange(1, cap + 2)
            ax.bar(xs, np.maximum(h, 0.8), width=1.0, color=MOD_COLOUR[code], label=MOD_NAME[code])
            ax.axvline(j * width + 0.5, color='black', lw=0.8)
        ax.set_yscale('log')
        log_axis(ax.yaxis)
        ticks = [j * width + t for j in range(len(present)) for t in (1, 25, cap + 1)]
        ax.set_xticks(ticks, [str(t) if t <= cap else f'>{cap}' for _ in present for t in (1, 25, cap + 1)], fontsize=8)
        ax.set_xlabel('modified sites per gene')
        ax.set_ylabel('genes')
        ax.legend(fontsize=8, loc='upper right', handlelength=1)
    else:
        ax.text(0.5, 0.5, 'no gene_mods.bed found', transform=ax.transAxes, ha='center')
    letter(ax, 'E')

    ax = fig.add_subplot(gs[1, 2])                                  # F co-occurrence
    present = [m for m in MAIN_MODS if genes.get(m)]
    if len(present) >= 2:
        gsets = {m: set(genes[m]) for m in present}
        w = 0.8 / (len(present) - 1)
        for i, obs in enumerate(present):
            givens = [g for g in present if g != obs]
            for j, giv in enumerate(givens):
                both = len(gsets[obs] & gsets[giv])
                p = both / len(gsets[giv]) if gsets[giv] else np.nan
                xpos = i - 0.4 + w * (j + 0.5)
                ax.bar(xpos, p, width=w, color=MOD_COLOUR[giv])
                ax.text(xpos, 0.03, f'{both:,} genes', rotation=90, ha='center', va='bottom', fontsize=7, color='white')
        handles = [plt.Rectangle((0, 0), 1, 1, color=MOD_COLOUR[m]) for m in present]
        ax.legend(handles, [MOD_NAME[m] for m in present], title='given', fontsize=9, title_fontsize=10,
                  loc='upper left', bbox_to_anchor=(1.0, 1.0), handlelength=1)
        ax.set_xticks(range(len(present)), [MOD_NAME[m] for m in present])
        ax.set_xlabel('observed modification')
        ax.set_ylabel('P(observed | given)')
        ax.set_ylim(0, 1.0)
    letter(ax, 'F')
    fig_batch_legend(fig, c.colours)
    notes.append(f'figure 2: {len(ms)} samples, {pct.shape[0]:,} cohort sites, '
                 f'{sum(len(v) for v in genes.values()):,} gene-modification pairs')
    return fig


def gene_sites(c):
    """mod code -> {gene_id: set of sites}, the union over samples of gene_mods.bed rows."""
    out = {}
    for s in c.samples:
        p = c.path('mods', c.ref, s, f'{s}.gene_mods.bed')
        if not os.path.exists(p):
            continue
        t = read_table(p, columns=None)
        if t.height == 0:
            continue
        for code, g, ch, st, sd in zip(text(t, 'mod_code'), text(t, 'gene_id'), text(t, 'chrom'),
                                       text(t, 'start'), text(t, 'strand')):
            out.setdefault(code, {}).setdefault(g, set()).add((ch, st, sd))
    return out


# ---------------------------------------------------------------- figure 3
def compare_groups(c, compare):
    names = sorted(set(c.batch.values()))
    if len(names) < 2:
        return None
    a = compare or names[-1]
    if a not in names:
        sys.exit(f'cohort_plots.py: --compare-batch {a} is not a batch ({", ".join(names)})')
    return a


def fig3_concordance(c, compare, notes):
    pm, ms, pct, mods = c.mods_matrix()
    a = compare_groups(c, compare)
    fig, gs = figure(3, 4)
    if a:
        ia = np.array([i for i, s in enumerate(ms) if c.batch[s] == a])
        ib = np.array([i for i, s in enumerate(ms) if c.batch[s] != a])
        with np.errstate(invalid='ignore'):
            meas_a = (~np.isnan(pct[:, ia])).sum(1)
            meas_b = (~np.isnan(pct[:, ib])).sum(1)
            mean_a = np.nanmean(np.where(np.isnan(pct[:, ia]), np.nan, pct[:, ia]), 1)
            mean_b = np.nanmean(np.where(np.isnan(pct[:, ib]), np.nan, pct[:, ib]), 1)
        keep = (meas_a >= min(2, len(ia))) & (meas_b >= min(2, len(ib)))
        vmax = 10                                                    # one colour scale for all eight
        for code in MAIN_MODS + OME_MODS:
            sel = keep & (mods == code) & ~np.isnan(mean_a) & ~np.isnan(mean_b)
            if sel.any():
                vmax = max(vmax, np.histogram2d(mean_b[sel], mean_a[sel], bins=40, range=[(0, 100), (0, 100)])[0].max())
        last, axes, rvals = None, [], {}
        for k, code in enumerate(MAIN_MODS + OME_MODS):
            ax = fig.add_subplot(gs[k // 4, k % 4])
            axes.append(ax)
            sel = keep & (mods == code)
            last = hist2d(ax, mean_b[sel], mean_a[sel], 'other batches, mean %', f'{a}, mean %', vmax=vmax) or last
            ok = sel & ~np.isnan(mean_a) & ~np.isnan(mean_b)
            rvals[code] = (np.corrcoef(mean_b[ok], mean_a[ok])[0, 1], int(ok.sum())) if ok.sum() >= 3 else (np.nan, 0)
            ax.set_title(MOD_NAME[code], color=MOD_COLOUR[code], fontweight='bold')
            letter(ax, 'ABCDEFGH'[k])
        if last is not None:
            colourbar(fig, last, axes, label='sites')
        ax = fig.add_subplot(gs[2, 2:])                                 # J the r of each modification
        codes = [m for m in MAIN_MODS + OME_MODS if not np.isnan(rvals.get(m, (np.nan, 0))[0])]
        ax.bar(range(len(codes)), [rvals[m][0] for m in codes], color=[MOD_COLOUR[m] for m in codes], width=0.7)
        for i, m in enumerate(codes):
            ax.text(i, rvals[m][0] + 0.01, f'{rvals[m][0]:.2f}', ha='center', va='bottom', fontsize=9)
            ax.text(i, 0.02, f'{rvals[m][1]:,}', ha='center', va='bottom', fontsize=7, color='white', rotation=90)
        ax.set_xticks(range(len(codes)), [MOD_NAME[m] for m in codes])
        ax.set_ylim(0, 1.08)
        ax.set_ylabel(f'r, {a} vs others')
        ax.set_title('per modification (sites in bars)', fontsize=11)
        letter(ax, 'J')
        notes.append(f'figure 3: {a} ({len(ia)} samples) vs the other batches ({len(ib)} samples); sites measured '
                     f'in >= 2 samples of each side')
    else:
        notes.append('figure 3: one batch only, so panels A-H are left empty')

    ax = fig.add_subplot(gs[2, :2])                                 # I each sample vs the rest
    qc = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.read_qc.tsv'))
    depth = dict(zip(text(qc, 'sample'), (numbers(qc, ['non_globin_reads']) - numbers(qc, ['mito_reads'])).ravel()))
    rs, ds, cs, ns = [], [], [], []
    for i, s in enumerate(ms):
        others = np.delete(pct, i, axis=1)
        with np.errstate(invalid='ignore'):
            ref = np.nanmedian(others, 1)
        ok = ~np.isnan(pct[:, i]) & ~np.isnan(ref)
        if ok.sum() >= 10:
            rs.append(np.corrcoef(pct[ok, i], ref[ok])[0, 1])
            ds.append(depth.get(s, np.nan) / 1e6)
            cs.append(c.colours[c.batch[s]])
            ns.append(s)
    ax.scatter(ds, rs, c=cs, s=36, zorder=3)
    for s, x_, y_ in zip(ns, ds, rs):
        ax.annotate(s, (x_, y_), fontsize=8, xytext=(4, 2), textcoords='offset points')
    ax.set_xlabel('non-globin reads (M)')
    ax.set_ylabel('r with the rest of the cohort')
    ax.set_title('each sample vs the median of the others, all modifications', fontsize=11)
    batch_legend(ax, c.colours, loc='lower right', fontsize=8, ncol=2)
    letter(ax, 'I')
    return fig


# ---------------------------------------------------------------- figure 4
def polya_hists(c, notes):
    """class -> {sample: (pt values, counts)} for nuclear without globin, globin, mito."""
    out = {k: {} for k, _, _ in COMPARTMENT}
    for s in c.samples:
        rp = c.path('polya', c.ref, s, f'{s}.polya_reads.tsv.gz')
        fp = c.path('polya', c.ref, s, f'{s}.polya_features.tsv.gz')
        if not os.path.exists(rp):
            continue
        globin_ids = set()
        if os.path.exists(fp):
            f = read_table(fp, columns=['level', 'feature_id', 'ref_gene_name'])
            globin_ids = {fid for lv, fid, nm in zip(text(f, 'level'), text(f, 'feature_id'), text(f, 'ref_gene_name'))
                          if lv == 'gene' and nm in GLOBIN}
        r = read_table(rp, columns=['compartment', 'pt', 'status', 'gene'])
        # classified column-wise in polars (millions of reads): only booleans and numbers reach numpy
        used = (r['status'] == 'used').fill_null(False).to_numpy().astype(bool)
        mito = (r['compartment'] == 'mito').fill_null(False).to_numpy().astype(bool)
        glob = r['gene'].is_in(sorted(globin_ids)).fill_null(False).to_numpy().astype(bool) & ~mito
        pt = numbers(r, ['pt']).ravel()
        used &= ~np.isnan(pt)
        masks = {'mito': mito, 'globin': glob, 'nuclear': ~mito & ~glob}
        for k, _, _ in COMPARTMENT:
            v = np.round(pt[used & masks[k]]).astype(int)
            if len(v):
                vals, cnt = np.unique(v, return_counts=True)
                out[k][s] = (vals, cnt)
        if not globin_ids:
            notes.append(f'{s}: no globin gene names in polya_features, so no globin class')
    return out


def hist_median(vals, cnt):
    cum = np.cumsum(cnt)
    return float(vals[np.searchsorted(cum, cum[-1] / 2)])


def fig4_polya(c, compare, notes):
    hists = polya_hists(c, notes)
    fig, gs = figure(1, 3, widths=[PANEL[0] * 1.25, PANEL[0], PANEL[0]])
    ax = fig.add_subplot(gs[0, 0])                                  # A pooled distributions
    top, meds = 0, []
    for k, label, col in COMPARTMENT:
        if not hists[k]:
            continue
        allv = np.concatenate([v for v, _ in hists[k].values()])
        allc = np.concatenate([n for _, n in hists[k].values()])
        vals, inv = np.unique(allv, return_inverse=True)
        cnt = np.bincount(inv, weights=allc)
        grid = np.arange(0, 301)
        dens = np.zeros(len(grid))
        inside = vals <= 300
        dens[vals[inside]] = cnt[inside]
        dens = np.convolve(dens / cnt.sum(), np.ones(5) / 5, mode='same')       # 5-nt smoothing
        ax.fill_between(grid, dens, color=col, alpha=0.3, lw=0)
        ax.plot(grid, dens, color=col, lw=1.2, label=f'{label}')
        med = hist_median(vals, cnt)
        ax.axvline(med, color=col, ls='--', lw=1)
        top = max(top, dens.max())
        meds.append((med, col))
    ax.set_xlim(0, 300)
    ax.set_ylim(0, top * 1.25 if top else 1)
    for j, (med, col) in enumerate(sorted(meds)):                   # labels above the curves, staggered
        ax.text(med + 3, top * (1.18 - 0.1 * j), f'{med:.0f}', color=col, fontsize=11, fontweight='bold', va='top')
    ax.set_xlabel('poly(A) length (nt)')
    ax.set_ylabel('density')
    ax.legend(fontsize=9, loc='upper right')
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B per-sample medians
    for j, (k, label, col) in enumerate(COMPARTMENT):
        xs = [i for i, s in enumerate(c.samples) if s in hists[k]]
        ys = [hist_median(*hists[k][c.samples[i]]) for i in xs]
        ax.plot(xs, ys, marker='o', ms=4, lw=0.8, color=col, label=label)
    sample_axis(ax, c.samples, c.batch, c.colours)
    ax.set_ylabel('median poly(A) (nt)')
    ax.set_ylim(bottom=0)
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C gene medians, batch vs rest
    a = compare_groups(c, compare)
    mp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.median_pt.tsv.gz')
    np_ = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.n_tails.tsv.gz')
    if a and os.path.exists(mp) and os.path.exists(np_):
        med, nt = read_table(mp), read_table(np_)
        ss = [s for s in c.samples if s in med.columns and s in nt.columns]
        mv, nv = numbers(med, ss), np.nan_to_num(numbers(nt, ss))
        mv[nv < 10] = np.nan                                         # the paper's >= 10 reads per gene
        ia = [i for i, s in enumerate(ss) if c.batch[s] == a]
        ib = [i for i, s in enumerate(ss) if c.batch[s] != a]
        with np.errstate(invalid='ignore'):
            ok = ((~np.isnan(mv[:, ia])).sum(1) >= min(2, len(ia))) & ((~np.isnan(mv[:, ib])).sum(1) >= min(2, len(ib)))
            ya, xb = np.nanmean(mv[ok][:, ia], 1), np.nanmean(mv[ok][:, ib], 1)
        hi = float(np.nanpercentile(np.concatenate([xb, ya]), 99.5)) if ok.any() else 200
        m = hist2d(ax, xb, ya, 'other batches (nt)', f'{a} (nt)', bins=40, rng=(0, math.ceil(hi / 10) * 10), cmap='magma')
        if m is not None:
            colourbar(fig, m, ax, label='genes')
        ax.set_title('gene medians (≥ 10 tails)', fontsize=11)
    else:
        ax.text(0.5, 0.5, 'needs two batches and the\ngene poly(A) matrices', transform=ax.transAxes, ha='center')
    letter(ax, 'C')
    fig_batch_legend(fig, c.colours)
    return fig


# ---------------------------------------------------------------- figure 5
def expression(c):
    cpm = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.gene_cpm.tsv.gz'))
    es = [s for s in c.samples if s in cpm.columns]
    m = np.nan_to_num(numbers(cpm, es))
    keep = np.median(m, axis=1) > 1
    names = [n if n not in ('.', '', None) else g for g, n in zip(text(cpm, 'gene_id'), text(cpm, 'gene_name'))]
    return es, np.log2(m[keep] + 1), [n for n, k in zip(names, keep) if k]


def fig5_expression(c, notes):
    es, lx, _ = expression(c)
    fig, gs = figure(1, 3, widths=[PANEL[0] * 1.6, PANEL[0] * 1.1, PANEL[0] * 1.2], heights=[PANEL[1] * 1.6])
    ax = fig.add_subplot(gs[0, 0])                                  # A correlation matrix
    r = np.corrcoef(lx, rowvar=False)
    off = r[~np.eye(len(r), dtype=bool)]
    shown = r.copy()
    np.fill_diagonal(shown, np.nan)
    cmap = MATRIX_CMAP.copy()
    cmap.set_bad('#FBF3EE')
    n = len(es)
    im = ax.pcolormesh(np.arange(n + 1) - 0.5, np.arange(n + 1) - 0.5, np.ma.masked_invalid(shown), cmap=cmap,
                       vmin=np.nanmin(off), vmax=np.nanmax(off), edgecolors='face', linewidth=0.3)
    for i in range(n):
        for j in range(n):
            v = r[i, j]
            ax.text(j, i, '1' if i == j else f'{v:.3f}'.lstrip('0'), ha='center', va='center', fontsize=6,
                    color='black' if (i == j or v > (np.nanmin(off) + np.nanmax(off)) / 2) else 'white')
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect('equal')
    sample_axis(ax, es, c.batch, c.colours, fontsize=8)
    ax.set_yticks(range(n), es, fontsize=8)
    for lab, s in zip(ax.get_yticklabels(), es):
        lab.set_color(c.colours[c.batch[s]])
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    colourbar(fig, im, ax, label='Pearson r, log2(CPM + 1)')
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B PCA
    keep = np.argsort(lx.var(axis=1))[::-1][:min(2000, len(lx))]
    xc = lx[keep] - lx[keep].mean(axis=1, keepdims=True)
    _, sv, vt = np.linalg.svd(xc, full_matrices=False)
    pcs, frac = vt.T * sv, sv ** 2 / (sv ** 2).sum()
    for i, s in enumerate(es):
        ax.scatter(pcs[i, 0], pcs[i, 1], color=c.colours[c.batch[s]], s=40, zorder=3)
        ax.annotate(s[-2:], (pcs[i, 0], pcs[i, 1]), fontsize=8, xytext=(4, 3), textcoords='offset points')
    ax.set_xlabel(f'PC1 ({100 * frac[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({100 * frac[1]:.1f}%)')
    ax.set_box_aspect(1)
    batch_legend(ax, c.colours, loc='best', fontsize=8)
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C isoform categories
    cp = c.path('isoforms', c.ref, 'high_confidence', f'cohort.{c.ref}.high_confidence.classification.tsv')
    if os.path.exists(cp):
        cats = text(read_table(cp), 'structural_category')
        names, counts = np.unique(cats, return_counts=True)
        order = np.argsort(counts)
        known = {'full-splice_match', 'incomplete-splice_match'}
        novel = {'novel_in_catalog', 'novel_not_in_catalog'}
        cols = ['#5769AB' if names[i] in known else '#EDB157' if names[i] in novel else '0.7' for i in order]
        ax.barh(range(len(names)), counts[order], color=cols)
        ax.set_yticks(range(len(names)), [names[i].replace('_', ' ') for i in order], fontsize=9)
        for i, v in enumerate(counts[order]):
            ax.text(v, i, f' {v:,}', va='center', fontsize=8)
        ax.set_xscale('log')
        ax.set_xlabel(f'high-confidence isoforms ({len(cats):,})')
        ax.set_xlim(right=counts.max() * 8)
    letter(ax, 'C')
    notes.append(f'figure 5: {lx.shape[0]:,} genes with median CPM > 1')
    return fig


# ---------------------------------------------------------------- figure 6
def volcano(ax, effect, p, q, labels, xlabel, colour, q_max, top=8):
    lp = -np.log10(np.clip(np.nan_to_num(p, nan=1.0), 1e-300, 1))
    sig = np.nan_to_num(q, nan=1.0) < q_max
    ax.scatter(effect[~sig], lp[~sig], s=5, color='0.75', lw=0)
    ax.scatter(effect[sig], lp[sig], s=7, color=colour, lw=0)
    for i in np.argsort(np.nan_to_num(p, nan=1.0))[:top]:
        ax.annotate(labels[i], (effect[i], lp[i]), fontsize=7, xytext=(3, 2), textcoords='offset points')
    ax.axvline(0, color='0.5', lw=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('-log$_{10}$ p (Welch)')
    stat_box(ax, [f'{int(sig.sum()):,} of {int(np.sum(~np.isnan(p))):,}', f'FDR < {q_max:g}'], loc='upper left')


def empty(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center', fontsize=10)
    ax.set_axis_off()


def fig6_groups(c, group, levels, fdr, out, notes):
    if group not in c.sheet.columns:
        sys.exit(f'cohort_plots.py: no column {group!r} in the samplesheet')
    grp = dict(zip(text(c.sheet, 'sample_id'), text(c.sheet, group)))
    levels = levels or sorted({v for v in grp.values() if v not in ('', '.', None)})
    if len(levels) != 2:
        sys.exit(f'cohort_plots.py: --group {group} has levels {levels}; pick two with --levels A B')
    a, b = levels
    label, tag = f'{b} vs {a}', f'{group}.{b}_vs_{a}'

    def idx(names):
        return (np.array([i for i, s in enumerate(names) if grp.get(s) == a]),
                np.array([i for i, s in enumerate(names) if grp.get(s) == b]))

    fig, gs = figure(2, 2, widths=[PANEL[0] * 1.2, PANEL[0] * 1.6], heights=[PANEL[1] * 1.3, PANEL[1] * 1.3])
    es, lx, names = expression(c)
    ga, gb = idx(es)
    if min(len(ga), len(gb)) < 2:
        sys.exit(f'cohort_plots.py: need at least 2 samples per level; {a}: {len(ga)}, {b}: {len(gb)}')
    lfc, p = welch(lx, ga, gb)
    q = bh(p)
    ax = fig.add_subplot(gs[0, 0])
    volcano(ax, lfc, p, q, names, f'log$_2$ fold change, {label}', '#E76969', fdr)
    letter(ax, 'A')
    write_rows(os.path.join(out, f'expression.{tag}.tsv'), ['gene', 'log2fc', 'p', 'q'],
               sorted(zip(names, lfc, p, q), key=lambda t: (math.isnan(t[2]), t[2])))

    ax = fig.add_subplot(gs[0, 1])
    top = list(np.argsort(np.nan_to_num(p, nan=1.0))[:30])
    cols = list(ga) + list(gb)
    z = lx[np.ix_(top, cols)]
    z = (z - z.mean(1, keepdims=True)) / (z.std(1, keepdims=True) + 1e-9)
    im = ax.pcolormesh(np.arange(len(cols) + 1) - 0.5, np.arange(len(top) + 1) - 0.5, z, cmap='RdBu_r',
                       vmin=-2.5, vmax=2.5, edgecolors='face', linewidth=0.3)
    ax.set_ylim(len(top) - 0.5, -0.5)
    ax.set_xticks(range(len(cols)), [es[i] for i in cols], rotation=90, fontsize=8)
    ax.set_yticks(range(len(top)), [names[i] for i in top], fontsize=7)
    ax.axvline(len(ga) - 0.5, color='black', lw=1)
    ax.text((len(ga) - 1) / 2, -1, a, ha='center', fontsize=9)
    ax.text(len(ga) + (len(gb) - 1) / 2, -1, b, ha='center', fontsize=9)
    ax.tick_params(length=0)
    colourbar(fig, im, ax, label='z-score')
    letter(ax, 'B')

    pm, ms, pct, mods = c.mods_matrix()
    ma, mb = idx(ms)
    min_n = min(3, len(ma), len(mb))
    d_all, p_all = welch(pct, ma, mb, min_n)
    full = ~np.isnan(p_all)                                         # >= 3 values per group
    ax = fig.add_subplot(gs[1, 0])
    if full.sum() >= 10:
        pf, mf = pct[full], mods[full]
        d, p = d_all[full], p_all[full]
        q = bh(p)
        lab = [f'{ch}:{int(float(st)) + 1} {MOD_NAME.get(mo, mo)}' for ch, st, mo, f in
               zip(text(pm, 'chrom'), text(pm, 'start'), text(pm, 'mod'), full) if f]
        sig = np.nan_to_num(q, nan=1.0) < fdr
        lp = -np.log10(np.clip(np.nan_to_num(p, nan=1.0), 1e-300, 1))
        ax.scatter(d[~sig], lp[~sig], s=5, color='0.75', lw=0)
        for code in MOD_NAME:
            s_ = sig & (mf == code)
            if s_.any():
                ax.scatter(d[s_], lp[s_], s=8, color=MOD_COLOUR[code], lw=0, label=MOD_NAME[code])
        for i in np.argsort(np.nan_to_num(p, nan=1.0))[:6]:
            ax.annotate(lab[i], (d[i], lp[i]), fontsize=6, xytext=(3, 2), textcoords='offset points')
        ax.axvline(0, color='0.5', lw=0.6)
        ax.set_xlabel(f'Δ occupancy (%), {label}')
        ax.set_ylabel('-log$_{10}$ p (Welch)')
        stat_box(ax, [f'{int(sig.sum()):,} of {len(p):,}', f'FDR < {fdr:g}', f'≥ {min_n} samples per group'])
        if sig.any():
            ax.legend(fontsize=7, loc='upper right', ncol=2, handlelength=1)
        write_rows(os.path.join(out, f'modifications.{tag}.tsv'), ['site', 'delta_pct', 'p', 'q'],
                   sorted(zip(lab, d, p, q), key=lambda t: (math.isnan(t[2]), t[2])))
    else:
        empty(ax, 'fewer than 10 sites measured in\n≥ 3 samples of each group')
    letter(ax, 'C')

    ax = fig.add_subplot(gs[1, 1])
    mp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.median_pt.tsv.gz')
    if os.path.exists(mp):
        med = read_table(mp)
        ss = [s for s in c.samples if s in med.columns]
        pv = numbers(med, ss)
        ntp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.n_tails.tsv.gz')
        if os.path.exists(ntp):
            pv[np.nan_to_num(numbers(read_table(ntp), ss)) < 10] = np.nan     # >= 10 tails, as in the paper
        pa, pb = idx(ss)
        d_all, p_all = welch(pv, pa, pb, min(3, len(pa), len(pb)))
        ok = ~np.isnan(p_all)
        if ok.sum() >= 10:
            gnames = [n if n not in ('.', '', None) else g for g, n, f in
                      zip(text(med, 'gene_id'), text(med, 'ref_gene_name'), ok) if f]
            d, p = d_all[ok], p_all[ok]
            q = bh(p)
            volcano(ax, d, p, q, gnames, f'Δ median poly(A) (nt), {label}', '#5769AB', fdr)
            write_rows(os.path.join(out, f'polya.{tag}.tsv'), ['gene', 'delta_median_pt', 'p', 'q'],
                       sorted(zip(gnames, d, p, q), key=lambda t: (math.isnan(t[2]), t[2])))
        else:
            empty(ax, 'fewer than 10 genes with ≥ 10 tails\nin ≥ 3 samples of each group')
    else:
        empty(ax, 'no gene poly(A) matrix')
    letter(ax, 'D')
    notes.append(f'figure 6: {label}, {len(ga)} vs {len(gb)} samples; Welch + BH, exploratory')
    return fig


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True)
    ap.add_argument('--samplesheet', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ref', default='GRCh38')
    ap.add_argument('--batch-column', default='batch')
    ap.add_argument('--compare-batch', help='the batch compared with all others in figures 3 and 4C (default: the last)')
    ap.add_argument('--group', help='samplesheet column to compare (two levels): adds figure 6')
    ap.add_argument('--levels', nargs=2, metavar=('A', 'B'))
    ap.add_argument('--fdr', type=float, default=0.1)
    ap.add_argument('--style', default=LAB_STYLE, help=f'matplotlib style sheet (default: {LAB_STYLE}); "" for none')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='(Mean of empty slice|All-NaN slice|Degrees of freedom)')
    notes = apply_style(args.style)
    c = Cohort(args.results, args.ref, args.samplesheet, args.batch_column)
    made = []
    figures = [('Fig1_cohort_qc', lambda: fig1_qc(c, notes)),
               ('Fig2_modification_calls', lambda: fig2_mods(c, notes)),
               ('Fig3_concordance', lambda: fig3_concordance(c, args.compare_batch, notes)),
               ('Fig4_polya', lambda: fig4_polya(c, args.compare_batch, notes)),
               ('Fig5_expression_isoforms', lambda: fig5_expression(c, notes))]
    if args.group:
        figures.append(('Fig6_group_comparison', lambda: fig6_groups(c, args.group, args.levels, args.fdr, args.out, notes)))
    for name, make in figures:
        fig = make()
        fig.savefig(os.path.join(args.out, f'{name}.pdf'), format='pdf', bbox_inches='tight', transparent=True)
        plt.close(fig)
        made.append(name)
    for n in notes:
        print(f'note: {n}', file=sys.stderr)
    print(f'cohort_plots.py: {len(made)} figures in {args.out}: {" ".join(f"{m}.pdf" for m in made)}')


if __name__ == '__main__':
    main()
