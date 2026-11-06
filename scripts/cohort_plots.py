#!/usr/bin/env python3
"""Plots of an RNA_analysis.nf cohort: one PDF per figure plus all of them in cohort_plots.pdf,
drawn in the Genometech lab style and editable in Illustrator (text stays text, as TrueType;
heatmap cells, points and bars are vector shapes, nothing rasterised).

    cohort_plots.py --results cohort-v1/results --samplesheet analysis_samplesheet.csv --out plots/
    cohort_plots.py ... --group condition [--levels control case]
    cohort_plots.py ... --style other.mplstyle   (default: the lab's genometechlab_main.mplstyle)

The style asks for Helvetica, which Linux clusters rarely have; then the closest installed
sans-serif (Arial, then Liberation Sans, which has Helvetica's metrics) is used and named, and
the font can be swapped in Illustrator.

Always (is the cohort comparable?):
  01 read composition per sample (globin, mitochondrial, other)
  02 non-globin reads against final modification sites
  03 sample correlation of expression, log2(CPM + 1), clustered, with batch
  04 PCA of expression (most variable genes), coloured by batch, sized by depth
  05 sample correlation of modification levels at sites measured in every sample
  06 modification levels per sample, one panel per modification type (same sites)
  07 poly(A) length distribution per sample, one panel per compartment
  08 structural categories of the high-confidence isoforms
  09 variant calls per sample by class
With --group (a samplesheet column with two levels, or --levels to pick two):
  10 expression volcano, 11 heatmap of the top genes, 12 modification-site volcano,
  13 poly(A) volcano (gene median tail length); each with a TSV of every test.
  Welch's t-test with Benjamini-Hochberg FDR: exploratory, for looking. For inference on
  expression use DESeq2 or edgeR on the counts matrix.

Empty cells in the modification and poly(A) matrices mean not measured, so those comparisons
use only features measured in every sample. polars reads the tables (as text, so no type
depends on the first rows); numpy and scipy do the numbers; matplotlib draws.
"""
import argparse
import gzip
import math
import os
import sys
import textwrap

import numpy as np
import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages            # noqa: E402
from matplotlib import font_manager                              # noqa: E402
import logging                                                   # noqa: E402
from scipy import stats                                          # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage        # noqa: E402
from scipy.spatial.distance import squareform                    # noqa: E402

LAB_STYLE = '/projects/Genometechlab/matplotlib_stylesheets/styles/genometechlab_main.mplstyle'
FONT_FALLBACKS = ['Helvetica', 'Arial', 'Liberation Sans', 'Nimbus Sans', 'TeX Gyre Heros']

MOD_NAMES = {'a': 'm6A', '17596': 'inosine', '69426': "2'-O-Me-A", 'm': 'm5C', '19228': "2'-O-Me-C",
             '19229': "2'-O-Me-G", '17802': 'pseudouridine', '19227': "2'-O-Me-U"}


# ---------------------------------------------------------------- tables
def read_table(path, sep='\t'):
    data = gzip.open(path).read() if path.endswith('.gz') else open(path, 'rb').read()
    return pl.read_csv(data, separator=sep, infer_schema_length=0)


def numbers(df, cols):
    """Columns as a float matrix; empty or non-numeric cells become NaN."""
    return df.select([pl.col(c).cast(pl.Float64, strict=False) for c in cols]).to_numpy().astype(float)


def write_rows(path, header, rows):
    with open(path, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for r in rows:
            out.write('\t'.join('NA' if (isinstance(v, float) and math.isnan(v)) else
                                (f'{v:.6g}' if isinstance(v, float) else str(v)) for v in r) + '\n')


def bh(p):
    """Benjamini-Hochberg q-values; NaN stays NaN."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return q
    ps = p[ok]
    order = np.argsort(ps)
    ranked = ps[order] * len(ps) / np.arange(1, len(ps) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ps)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


# ---------------------------------------------------------------- style
def apply_style(path):
    """The lab style, then the first of FONT_FALLBACKS that is installed. Returns a note."""
    notes = []
    if path and os.path.exists(path):
        plt.style.use(path)
        notes.append(f'style: {path}')
    elif path:
        notes.append(f'style {path} not found: matplotlib defaults')
    plt.rcParams['pdf.fonttype'] = 42                    # text stays editable text in Illustrator
    plt.rcParams['ps.fonttype'] = 42
    installed = {f.name for f in font_manager.fontManager.ttflist}
    wanted = plt.rcParams['font.sans-serif'][0] if plt.rcParams['font.sans-serif'] else 'Helvetica'
    chosen = next((f for f in [wanted] + FONT_FALLBACKS if f in installed), 'DejaVu Sans')
    plt.rcParams['font.sans-serif'] = [chosen, 'DejaVu Sans']
    for k in ('mathtext.rm', 'mathtext.sf'):
        plt.rcParams[k] = chosen
    plt.rcParams['mathtext.it'] = f'{chosen}:italic'
    plt.rcParams['mathtext.bf'] = f'{chosen}:bold'
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)   # no per-glyph fallback warnings
    notes.append(f'font: {chosen}' + ('' if chosen == wanted else f' ({wanted} is not installed)'))
    return notes


def tidy(fig):
    """The style's constrained layout when it is on, else tight_layout."""
    if fig.get_layout_engine() is None:
        tidy(fig)


def wrap(text, width=50):
    """Long titles on several lines: the lab style's 16 pt bold runs past a plot otherwise."""
    return textwrap.fill(text, width)


def palette():
    return plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'])


# ---------------------------------------------------------------- drawing helpers
class Figures:
    def __init__(self, out):
        self.out = out
        self.pdf = PdfPages(os.path.join(out, 'cohort_plots.pdf'))
        self.made = []

    def save(self, fig, name):
        fig.savefig(os.path.join(self.out, f'{name}.pdf'), format='pdf', bbox_inches='tight')
        self.pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        self.made.append(name)

    def close(self):
        self.pdf.close()


def batch_colours(batches):
    names = sorted(set(batches))
    pal = palette()
    return {b: pal[i % len(pal)] for i, b in enumerate(names)}


def legend_for(ax, colours, title, **kw):
    handles = [plt.Line2D([], [], marker='o', ls='', color=c, label=b) for b, c in colours.items()]
    ax.legend(handles=handles, title=title, frameon=False, **kw)


def clustered_heatmap(r, samples, batch, colours, title, label):
    order = leaves_list(linkage(squareform(np.clip(1 - r, 0, None), checks=False), 'average')) if len(samples) > 2 \
        else np.arange(len(samples))
    rr = r[np.ix_(order, order)]
    names = [samples[i] for i in order]
    fig, ax = plt.subplots(figsize=(0.6 * len(names) + 4, 0.6 * len(names) + 2))
    off = r[~np.eye(len(r), dtype=bool)]
    lo, hi = np.nanmin(off), np.nanmax(off)          # the scale spans the samples' pairs, not the 1.0 diagonal
    shown = rr.copy()
    np.fill_diagonal(shown, np.nan)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('0.85')
    n = len(names)
    im = ax.pcolormesh(np.arange(n + 1) - 0.5, np.arange(n + 1) - 0.5, np.ma.masked_invalid(shown),
                       cmap=cmap, vmin=lo, vmax=hi, edgecolors='face', linewidth=0.3)   # vector cells; 'face' closes seams
    ax.set_aspect('equal')
    ax.set_xticks(range(n), names, rotation=90)
    ax.set_yticks(range(n), names)
    for i, s in enumerate(names):                                # batch strips along both axes
        ax.add_patch(plt.Rectangle((i - 0.5, -1.4), 1, 0.7, color=colours[batch[s]], clip_on=False))
        ax.add_patch(plt.Rectangle((-1.4, i - 0.5), 0.7, 1, color=colours[batch[s]], clip_on=False))
    if len(names) <= 20:
        for i in range(len(names)):
            for j in range(len(names)):
                if i != j:
                    ax.text(j, i, f'{rr[i, j]:.3f}', ha='center', va='center', fontsize=5.5,
                            color='white' if rr[i, j] < (lo + hi) / 2 else 'black')
    ax.set_xlim(-1.5, n - 0.5)
    ax.set_ylim(n - 0.5, -1.5)
    for side in ('top', 'right', 'bottom', 'left'):
        ax.spines[side].set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03, label=label)
    cb.solids.set_rasterized(False)        # matplotlib rasterises gradients of 50+ steps; keep it vector
    cb.solids.set_edgecolor('face')        # and close the seams between its steps
    legend_for(ax, colours, 'batch', loc='upper left', bbox_to_anchor=(1.22, 1))
    ax.set_title(wrap(title))
    return fig


def volcano(ax, effect, p, q, labels, xlabel, colour=None, q_max=0.1, top=10):
    lp = -np.log10(np.clip(p, 1e-300, 1))
    sig = q < q_max
    ax.scatter(effect[~sig], lp[~sig], s=6, color='0.7', lw=0)                   # vector points, editable
    ax.scatter(effect[sig], lp[sig], s=8, color=colour if colour is not None else palette()[5], lw=0)
    for i in np.argsort(p)[:top]:
        ax.annotate(labels[i], (effect[i], lp[i]), fontsize=6, xytext=(3, 2), textcoords='offset points')
    ax.axvline(0, color='0.5', lw=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('-log10 p (Welch)')
    ax.set_title(f'{int(sig.sum())} of {len(p)} at FDR < {q_max:g}', fontsize=9)


def welch(x, ga, gb):
    """Welch's t-test per row of x between sample index sets ga and gb."""
    with np.errstate(invalid='ignore', divide='ignore'):
        t, p = stats.ttest_ind(x[:, gb], x[:, ga], axis=1, equal_var=False)
    return np.asarray(p, float)


# ---------------------------------------------------------------- the figures
def fig_read_composition(globin, mito, non_globin, samples, batch, colours):
    g, mt = np.array(globin) / 1e6, np.array(mito) / 1e6
    other = (np.array(non_globin) - np.array(mito)) / 1e6        # non_globin counts mitochondrial reads too
    fig, ax = plt.subplots(figsize=(0.5 * len(samples) + 2, 4))
    x = np.arange(len(samples))
    pal = palette()
    ax.bar(x, other, color=pal[0], label='other')
    ax.bar(x, mt, bottom=other, color=pal[4], label='mitochondrial')
    ax.bar(x, g, bottom=other + mt, color='0.8', label='globin')
    ax.set_xticks(x, samples, rotation=90)
    for lab, s in zip(ax.get_xticklabels(), samples):
        lab.set_color(colours[batch[s]])
    for i in range(len(samples)):
        ax.text(i, other[i] / 2, f'{other[i]:.1f}', ha='center', va='center', fontsize=6, color='white', rotation=90)
    ax.set_ylabel('mapped reads (millions)')
    ax.set_title(wrap('Read composition (sample names coloured by batch; white: other reads, M)'))
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1))
    return fig


def fig_depth_vs_sites(qc_samples, other_reads, msum, batch, colours):
    sites = dict(zip(msum['sample'].to_list(), numbers(msum, ['final_sites']).ravel()))
    xs = np.array(other_reads) / 1e6
    ys = np.array([sites.get(s, np.nan) for s in qc_samples]) / 1e3
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for s, x, y in zip(qc_samples, xs, ys):
        ax.scatter(x, y, color=colours[batch[s]], s=40)
        ax.annotate(s, (x, y), fontsize=7, xytext=(4, 2), textcoords='offset points')
    ok = ~np.isnan(ys)
    if ok.sum() >= 3:
        k, c = np.polyfit(xs[ok], ys[ok], 1)
        xx = np.linspace(0, xs.max() * 1.05, 10)
        ax.plot(xx, k * xx + c, color='0.6', lw=0.8, ls='--', label='least squares')
    ax.set_xlabel('non-globin reads (millions)')
    ax.set_ylabel('final modification sites (thousands)')
    ax.set_title('Depth against modification sites')
    legend_for(ax, colours, 'batch', loc='lower right')
    return fig


def fig_pca(x, samples, batch, colours, depth, top=2000):
    var = x.var(axis=1)
    keep = np.argsort(var)[::-1][:min(top, len(var))]
    xc = x[keep] - x[keep].mean(axis=1, keepdims=True)
    _, s, vt = np.linalg.svd(xc, full_matrices=False)
    pcs = vt.T * s
    frac = s ** 2 / (s ** 2).sum()
    fig, ax = plt.subplots(figsize=(6, 5))
    size = 30 + 170 * (np.array(depth) - min(depth)) / (max(depth) - min(depth) + 1e-9)
    for i, sm in enumerate(samples):
        ax.scatter(pcs[i, 0], pcs[i, 1], color=colours[batch[sm]], s=size[i], alpha=0.85)
        ax.annotate(sm, (pcs[i, 0], pcs[i, 1]), fontsize=7, xytext=(4, 3), textcoords='offset points')
    ax.set_xlabel(f'PC1 ({100 * frac[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({100 * frac[1]:.1f}%)' if len(frac) > 1 else 'PC2')
    ax.set_title(wrap(f'PCA of expression, {len(keep)} most variable genes (point size: non-globin depth)'))
    legend_for(ax, colours, 'batch', loc='best')
    return fig


def fig_mod_distributions(pct, mods, samples, batch, colours):
    codes = [c for c in MOD_NAMES if c in set(mods)] + sorted(set(mods) - set(MOD_NAMES))
    n = len(codes)
    cols = min(4, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.2 * rows), squeeze=False)
    for k, code in enumerate(codes):
        ax = axes[k // cols][k % cols]
        sel = np.array([m == code for m in mods])
        data = [pct[sel, j] for j in range(len(samples))]
        bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6, medianprops={'color': 'black'})
        for patch, s in zip(bp['boxes'], samples):
            patch.set_facecolor(colours[batch[s]])
            patch.set_alpha(0.8)
        ax.set_xticks(range(1, len(samples) + 1), samples, rotation=90, fontsize=6)
        ax.set_title(f'{MOD_NAMES.get(code, code)} ({int(sel.sum())} sites)', fontsize=9)
        ax.set_ylabel('% modified (gIVT-corrected)', fontsize=8)
    for k in range(n, rows * cols):
        axes[k // cols][k % cols].axis('off')
    fig.suptitle(wrap('Modification levels at sites measured in every sample (boxes coloured by batch)'))
    tidy(fig)
    return fig


def fig_polya(hists, samples, batch, colours):
    comps = sorted({c for h in hists.values() for c in h})
    fig, axes = plt.subplots(1, len(comps), figsize=(5.5 * len(comps), 4), squeeze=False)
    for k, comp in enumerate(comps):
        ax = axes[0][k]
        pooled = []
        for s in samples:
            if comp not in hists.get(s, {}):
                continue
            pt, reads = hists[s][comp]
            order = np.argsort(pt)
            pt, reads = pt[order], reads[order]
            total = reads.sum()
            if total == 0:
                continue
            ax.plot(pt, reads / total, color=colours[batch[s]], lw=0.9, alpha=0.85)
            pooled.append(np.repeat(pt, np.minimum(reads, 10000).astype(int)))
            med = pt[np.searchsorted(np.cumsum(reads), total / 2)]
            ax.axvline(med, color=colours[batch[s]], lw=0.4, alpha=0.6)
        if pooled:
            ax.set_xlim(0, np.percentile(np.concatenate(pooled), 99))
        ax.set_xlabel('poly(A) length (nt)')
        ax.set_ylabel('fraction of reads')
        ax.set_title(f'{comp} (thin lines: each sample\'s median)', fontsize=9)
    legend_for(axes[0][-1], colours, 'batch', loc='upper right')
    fig.suptitle('Poly(A) tail length per sample')
    tidy(fig)
    return fig


def fig_isoform_categories(cls):
    cats = cls['structural_category'].to_list()
    names, counts = np.unique(cats, return_counts=True)
    order = np.argsort(counts)
    fig, ax = plt.subplots(figsize=(6.5, 0.35 * len(names) + 1.5))
    ax.barh(range(len(names)), counts[order], color=palette()[0])
    ax.set_yticks(range(len(names)), names[order])
    for i, c in enumerate(counts[order]):
        ax.text(c, i, f' {c}', va='center', fontsize=7)
    ax.set_xlabel('high-confidence isoforms')
    ax.set_title(f'Structural categories ({len(cats)} isoforms)')
    return fig


def fig_variants(summaries, samples, batch, colours):
    parts = ['known', 'novel', 'editing', 'masked']
    pal = palette()
    colour = {'known': pal[0], 'novel': pal[5], 'editing': pal[1], 'masked': '0.75'}
    vals = {p: [] for p in parts}
    for s in samples:
        d = summaries.get(s, {})
        for p in parts[:3]:
            vals[p].append(d.get(p, 0))
        vals['masked'].append(sum(v for k, v in d.items() if k.startswith('masked:')))
    fig, ax = plt.subplots(figsize=(0.5 * len(samples) + 2, 4))
    x = np.arange(len(samples))
    bottom = np.zeros(len(samples))
    for p in parts:
        ax.bar(x, vals[p], bottom=bottom, color=colour[p], label=p)
        bottom += np.array(vals[p], float)
    ax.set_xticks(x, samples, rotation=90)
    for lab, s in zip(ax.get_xticklabels(), samples):
        lab.set_color(colours[batch[s]])
    ax.set_ylabel('variant calls (PASS classes)')
    ax.set_title(wrap('Variant calls by class (masked: IG/TR/HLA regions; names coloured by batch)'))
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1))
    return fig


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True, help='the run outdir (contains cohort/, mods/, ...)')
    ap.add_argument('--samplesheet', required=True)
    ap.add_argument('--out', required=True, help='folder for the PNGs, PDF and comparison tables')
    ap.add_argument('--ref', default='GRCh38')
    ap.add_argument('--batch-column', default='batch')
    ap.add_argument('--group', help='samplesheet column to compare (two levels)')
    ap.add_argument('--levels', nargs=2, metavar=('A', 'B'), help='the two levels of --group, A as reference')
    ap.add_argument('--fdr', type=float, default=0.1)
    ap.add_argument('--style', default=LAB_STYLE, help=f'matplotlib style sheet (default: {LAB_STYLE}); "" for none')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    style_notes = apply_style(args.style)
    R, ref = args.results, args.ref
    cohort = os.path.join(R, 'cohort', ref)

    sheet = read_table(args.samplesheet, ',')
    ids = sheet['sample_id'].to_list()
    bcol = sheet[args.batch_column].to_list() if args.batch_column in sheet.columns else ['all'] * len(ids)
    batch = dict(zip(ids, bcol))
    samples = sorted(ids, key=lambda s: (batch[s], s))            # grouped by batch on every axis
    colours = batch_colours(batch.values())
    figs = Figures(args.out)
    notes = []

    # 01-02 read QC and depth
    qc = read_table(os.path.join(cohort, f'cohort.{ref}.read_qc.tsv'))
    qrow = {s: i for i, s in enumerate(qc['sample'].to_list())}
    qc_samples = [s for s in samples if s in qrow]
    qnum = {c: numbers(qc, [c]).ravel() for c in ('globin_reads', 'mito_reads', 'non_globin_reads')}
    col = {c: [qnum[c][qrow[s]] for s in qc_samples] for c in qnum}
    figs.save(fig_read_composition(col['globin_reads'], col['mito_reads'], col['non_globin_reads'],
                                   qc_samples, batch, colours), '01_read_composition')
    other = [qnum['non_globin_reads'][qrow[s]] - qnum['mito_reads'][qrow[s]] for s in qc_samples]
    depth = dict(zip(qc_samples, other))
    msum = read_table(os.path.join(cohort, f'cohort.{ref}.mods.summary.tsv'))
    figs.save(fig_depth_vs_sites(qc_samples, other, msum, batch, colours), '02_depth_vs_sites')

    # 03-04 expression
    cpm = read_table(os.path.join(cohort, f'cohort.{ref}.gene_cpm.tsv.gz'))
    es = [s for s in samples if s in cpm.columns]
    m = numbers(cpm, es)
    m = np.nan_to_num(m, nan=0.0)
    keep = np.median(m, axis=1) > 1
    lx = np.log2(m[keep] + 1)
    gene_names = [n if n not in ('.', '', None) else g for g, n in zip(cpm['gene_id'].to_list(), cpm['gene_name'].to_list())]
    gene_names = [n for n, k in zip(gene_names, keep) if k]
    r = np.corrcoef(lx, rowvar=False)
    figs.save(clustered_heatmap(r, es, batch, colours,
                                f'Expression: sample correlation, log2(CPM + 1), {int(keep.sum())} genes with median CPM > 1',
                                'Pearson r'), '03_expression_correlation')
    figs.save(fig_pca(lx, es, batch, colours, [depth.get(s, 0) for s in es]), '04_expression_pca')

    # 05-06 modifications at sites measured in every sample
    pmat = read_table(os.path.join(cohort, f'cohort.{ref}.mods.pct_fpcorr.tsv.gz'))
    ms = [s for s in samples if s in pmat.columns]
    pct = numbers(pmat, ms)
    full = ~np.isnan(pct).any(axis=1)
    mods_all = pmat['mod'].to_list()
    mods = [mo for mo, f in zip(mods_all, full) if f]
    site_ids = [f'{c}:{int(float(st)) + 1}{sd} {MOD_NAMES.get(mo, mo)}' for c, st, sd, mo, f in
                zip(pmat['chrom'].to_list(), pmat['start'].to_list(), pmat['strand'].to_list(), mods_all, full) if f]
    pct_full = pct[full]
    notes.append(f'modification sites measured in all {len(ms)} samples: {int(full.sum())} of {len(full)}')
    if full.sum() >= 10:
        rm = np.corrcoef(pct_full, rowvar=False)
        figs.save(clustered_heatmap(rm, ms, batch, colours,
                                    f'Modification levels: sample correlation at {int(full.sum())} sites measured in every sample',
                                    'Pearson r'), '05_modification_correlation')
        figs.save(fig_mod_distributions(pct_full, mods, ms, batch, colours), '06_modification_levels')
    else:
        notes.append('fewer than 10 sites measured in every sample: figures 05-06 skipped')

    # 07 poly(A)
    hists = {}
    for s in samples:
        p = os.path.join(R, 'polya', ref, s, f'{s}.polya_hist.tsv')
        if not os.path.exists(p):
            continue
        h = read_table(p)
        comp = h['compartment'].to_list()
        v = numbers(h, ['pt', 'reads'])
        hists[s] = {c: (v[[i for i, x in enumerate(comp) if x == c], 0], v[[i for i, x in enumerate(comp) if x == c], 1])
                    for c in set(comp)}
    if hists:
        figs.save(fig_polya(hists, samples, batch, colours), '07_polya_lengths')

    # 08 isoform categories
    cp = os.path.join(R, 'isoforms', ref, 'high_confidence', f'cohort.{ref}.high_confidence.classification.tsv')
    if os.path.exists(cp):
        figs.save(fig_isoform_categories(read_table(cp)), '08_isoform_categories')

    # 09 variants
    summaries = {}
    for s in samples:
        p = os.path.join(R, 'variants', ref, s, f'{s}.variants.summary.tsv')
        if os.path.exists(p):
            t = read_table(p)
            summaries[s] = dict(zip(t['class'].to_list(), numbers(t, ['calls']).ravel()))
    if summaries:
        figs.save(fig_variants(summaries, samples, batch, colours), '09_variant_classes')

    # 10-13 two-group comparisons
    if args.group:
        if args.group not in sheet.columns:
            sys.exit(f'cohort_plots.py: no column {args.group!r} in {args.samplesheet}')
        grp = dict(zip(ids, sheet[args.group].to_list()))
        levels = args.levels or sorted({v for v in grp.values() if v not in ('', '.', None)})
        if len(levels) != 2:
            sys.exit(f'cohort_plots.py: --group {args.group} has levels {levels}; pick two with --levels A B')
        a, b = levels
        label = f'{b} vs {a}'

        def idx(names):
            return (np.array([i for i, s in enumerate(names) if grp.get(s) == a]),
                    np.array([i for i, s in enumerate(names) if grp.get(s) == b]))

        ga, gb = idx(es)
        if min(len(ga), len(gb)) < 2:
            sys.exit(f'cohort_plots.py: need at least 2 samples per level; {a}: {len(ga)}, {b}: {len(gb)}')
        if min(len(ga), len(gb)) < 3:
            notes.append(f'only {min(len(ga), len(gb))} samples in the smaller group: p-values are very rough')
        tag = f'{args.group}.{b}_vs_{a}'

        # 10-11 expression
        lfc = lx[:, gb].mean(axis=1) - lx[:, ga].mean(axis=1)
        p = welch(lx, ga, gb)
        q = bh(p)
        fig, ax = plt.subplots(figsize=(6, 5))
        volcano(ax, lfc, np.nan_to_num(p, nan=1.0), np.nan_to_num(q, nan=1.0), gene_names,
                f'log2 fold change, {label} (on log2 CPM + 1)', q_max=args.fdr)
        fig.suptitle(wrap(f'Expression: {label} (exploratory; use DESeq2/edgeR for inference)'), fontsize=10)
        figs.save(fig, '10_expression_volcano')
        write_rows(os.path.join(args.out, f'expression.{tag}.tsv'), ['gene', 'log2fc', 'p', 'q'],
                   sorted(zip(gene_names, lfc, p, q), key=lambda t: (math.isnan(t[2]), t[2])))
        top = [i for i in np.argsort(np.nan_to_num(p, nan=1.0))[:40]]
        cols = list(ga) + list(gb)
        z = lx[np.ix_(top, cols)]
        z = (z - z.mean(axis=1, keepdims=True)) / (z.std(axis=1, keepdims=True) + 1e-9)
        fig, ax = plt.subplots(figsize=(0.45 * len(cols) + 3, 0.22 * len(top) + 2))
        im = ax.pcolormesh(np.arange(z.shape[1] + 1) - 0.5, np.arange(z.shape[0] + 1) - 0.5, z,
                           cmap='RdBu_r', vmin=-2.5, vmax=2.5, edgecolors='face', linewidth=0.3)   # vector cells, no seams
        ax.set_ylim(len(top) - 0.5, -0.5)
        ax.set_xticks(range(len(cols)), [f'{es[i]} ({grp[es[i]]})' for i in cols], rotation=90, fontsize=7)
        ax.set_yticks(range(len(top)), [gene_names[i] for i in top], fontsize=6)
        ax.axvline(len(ga) - 0.5, color='k', lw=1)
        cb = fig.colorbar(im, ax=ax, fraction=0.03, label='z-score of log2(CPM + 1)')
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor('face')
        ax.set_title(wrap(f'Top {len(top)} genes by p, {label}'))
        figs.save(fig, '11_expression_top_genes')

        # 12 modification sites
        if full.sum() >= 10:
            ma, mb = idx(ms)
            d = pct_full[:, mb].mean(axis=1) - pct_full[:, ma].mean(axis=1)
            p = welch(pct_full, ma, mb)
            q = bh(p)
            fig, ax = plt.subplots(figsize=(6, 5))
            volcano(ax, d, np.nan_to_num(p, nan=1.0), np.nan_to_num(q, nan=1.0), site_ids,
                    f'difference in % modified, {label}', q_max=args.fdr)
            fig.suptitle(wrap(f'Modification sites measured in every sample: {label}'), fontsize=10)
            figs.save(fig, '12_modification_volcano')
            write_rows(os.path.join(args.out, f'modifications.{tag}.tsv'), ['site', 'delta_pct', 'p', 'q'],
                       sorted(zip(site_ids, d, p, q), key=lambda t: (math.isnan(t[2]), t[2])))

        # 13 poly(A) per gene
        pp = os.path.join(cohort, f'cohort.{ref}.polya.gene.median_pt.tsv.gz')
        if os.path.exists(pp):
            pg = read_table(pp)
            ps = [s for s in samples if s in pg.columns]
            pv = numbers(pg, ps)
            pfull = ~np.isnan(pv).any(axis=1)
            if pfull.sum() >= 10:
                pa, pb = idx(ps)
                pvf = pv[pfull]
                names = [n if n not in ('.', '', None) else g for g, n, f in
                         zip(pg['gene_id'].to_list(), pg['ref_gene_name'].to_list(), pfull) if f]
                d = pvf[:, pb].mean(axis=1) - pvf[:, pa].mean(axis=1)
                p = welch(pvf, pa, pb)
                q = bh(p)
                fig, ax = plt.subplots(figsize=(6, 5))
                volcano(ax, d, np.nan_to_num(p, nan=1.0), np.nan_to_num(q, nan=1.0), names,
                        f'difference in median poly(A) length (nt), {label}', q_max=args.fdr)
                fig.suptitle(wrap(f'Gene poly(A) length, {int(pfull.sum())} genes measured in every sample: {label}'), fontsize=10)
                figs.save(fig, '13_polya_volcano')
                write_rows(os.path.join(args.out, f'polya.{tag}.tsv'), ['gene', 'delta_median_pt', 'p', 'q'],
                           sorted(zip(names, d, p, q), key=lambda t: (math.isnan(t[2]), t[2])))

    figs.close()
    for n in style_notes + notes:
        print(f'note: {n}', file=sys.stderr)
    print(f'cohort_plots.py: {len(figs.made)} figures in {args.out} (one PDF each + cohort_plots.pdf): {" ".join(figs.made)}')


if __name__ == '__main__':
    main()
