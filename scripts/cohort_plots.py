#!/usr/bin/env python3
"""Cohort figures for RNA_analysis.nf results, laid out like Ghohabi Esfahani, Stein, Akeson et al.
(Genome Biology 2026, 27:190): bold panel letters, the lab style sheet, the paper's colours,
log-scaled 2D histograms, dashed medians. One vector PDF per figure, editable in Illustrator.

    cohort_plots.py --results <results> --samplesheet <analysis_samplesheet.csv> --out <dir>
                    [--mask-bed regions.bed] [--mask-contigs REGEX] [--read-stats <dir>]
                    [--group <column> [--levels A B] [--fdr 0.1]] [--compare-batch <batch>]
                    [--only Fig3 Fig5 ...] [--threads 8] [--ref GRCh38] [--style <mplstyle>]

  Fig1_read_qc                per-read identity (dashed medians), identity by aligned length, indels
                              per aligned base, substitution matrix per batch (paper Fig. 1). Needs
                              scripts/read_stats.sh outputs in --read-stats (default <results>/read_stats)
  Fig2_cohort_overview        read composition, non-globin depth, modification sites against depth,
                              expression PCA and correlation matrix
  Fig3_ivt_correction         Dorado occupancy against gIVT-corrected occupancy per modification, with
                              the 20% lines; sites >= 20% per sample before and after (paper Fig. 3,
                              Table 4); per-sample counts in ivt_correction_sites.tsv
  Fig4_modification_calls     sites per occupancy bin, samples calling each site, sites per gene,
                              co-occurrence of modifications in genes (paper Fig. 5A, B, D)
  Fig5_replicate_concordance  pairwise ICC(2,1) of occupancy per modification, ICC over all samples, ICC
                              within against between batches, and a depth-or-batch check (the paper's
                              HEK293 replicate analysis, over every pair of samples)
  Fig6_isoform_modifications  isoform-specific occupancy differences (paper Fig. 5C)
  Fig7_polya                  per-read and per-gene poly(A), per-sample medians, isoform-specific
                              tails, long- and short-tailed genes (paper Fig. S2)
  Fig8_group_comparison       with --group: expression, modification and poly(A) volcanoes, heatmap
  Fig9_isoform_usage          isoform usage between samples: variation beyond read sampling, PCA of
                              isoform fractions, per-sample isoform switching, example genes
  FigS1_batch_comparison      with --compare-batch: one batch's mean occupancy against the others'
  FigS2_jaccard               pairwise Jaccard of called sites, within against between batches, depth check

Next to the PDFs: figure_legends.md (a legend per figure, with its counts and whether each panel is pooled
over samples or per sample; also one file per figure in legends/) and tables/ (Table1_cohort_summary.tsv,
one value per sample, and the numbers behind each figure, named after it).

Modification figures (2C, 3-6, 8, S1, S2) leave out sites on contigs matching --mask-contigs (default:
unplaced, random, alt, fix, decoy, EBV and HLA contigs) and inside --mask-bed intervals; the notes
report what each rule removed. Sites are the cohort's final-stage sites (>= 20 reads, >= 20% after
gIVT correction); a sample calls a site where its corrected percent is >= 20, and an empty cell
means not measured there (below 20 reads, or no correction possible). Figure 4E counts a gene's
site only if no site of another modification lies within 5 nt. Figure 7's per-read poly(A) uses
reads assigned to a gene (mitochondrial reads by chromosome). Welch's t-test with Benjamini-
Hochberg FDR in figure 8 is exploratory; use DESeq2 or edgeR on the counts for expression.

Figure 3 streams each <sample>.cov20_fpcorr.bed through awk (tens of millions of rows each) and
caches the binned counts in <out>/.cache, so only the first run reads them. polars reads the
other tables (as text); numpy and scipy compute.
"""
import argparse
import gzip
import re
import logging
import math
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib import font_manager                              # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm   # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter         # noqa: E402
from scipy import stats                                          # noqa: E402

log = logging.getLogger('cohort_plots')
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


# ---------------------------------------------------------------- v4: data access, masks, statistics
SITE_FILTER_DEFAULT = r'(_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-)'


class SiteMask:
    """Sites kept for modification figures: not on contigs matching the regex (unplaced, random,
    alt...), and outside --mask-bed intervals. Counts what each rule removed from the cohort
    matrix, per mod."""

    def __init__(self, contig_regex, bed):
        self.pattern = contig_regex or ''
        self.rx = re.compile(contig_regex) if contig_regex else None
        self.iv, raw = {}, {}
        if bed:
            with (gzip.open(bed, 'rt') if bed.endswith('.gz') else open(bed)) as fh:
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    if len(f) < 3 or line.startswith(('#', 'track', 'browser')):
                        continue
                    try:
                        raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
                    except ValueError:
                        continue                                    # a header line
        for c_, v in raw.items():                                   # merged, sorted
            v.sort()
            merged = [list(v[0])]
            for s_, e_ in v[1:]:
                if s_ <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e_)
                else:
                    merged.append([s_, e_])
            a = np.array(merged, dtype=float)
            self.iv[c_] = (a[:, 0], a[:, 1])
        self.removed = {}

    def intervals(self):
        return [(c_, int(a), int(b)) for c_, (s_, e_) in sorted(self.iv.items()) for a, b in zip(s_, e_)]

    def keep(self, chroms, starts, mods=None):
        chroms = np.asarray(chroms, dtype=object)
        starts = np.asarray(starts, dtype=float)
        ok = np.ones(len(chroms), bool)
        if self.rx is not None:
            uniq = {c_: bool(self.rx.search(c_)) for c_ in set(chroms.tolist())}
            by_contig = np.array([uniq[c_] for c_ in chroms.tolist()], bool) if len(chroms) else np.zeros(0, bool)
            self._count('contig', by_contig, mods)
            ok &= ~by_contig
        if self.iv:
            in_bed = np.zeros(len(chroms), bool)
            for c_, (s_, e_) in self.iv.items():
                sel = np.where(chroms == c_)[0]
                if len(sel):
                    j = np.searchsorted(s_, starts[sel], side='right') - 1
                    hit = (j >= 0) & (starts[sel] < e_[np.clip(j, 0, None)])
                    in_bed[sel[hit]] = True
            self._count('mask-bed', in_bed & ok, mods)
            ok &= ~in_bed
        return ok

    def _count(self, rule, mask, mods):
        if mods is None:
            return
        for code in set(np.asarray(mods)[mask].tolist()):
            self.removed[(rule, code)] = self.removed.get((rule, code), 0) + int(np.sum(mask & (np.asarray(mods) == code)))


def icc21(x):
    """ICC(2,1), two-way random effects, absolute agreement, single measure (Shrout & Fleiss)."""
    n, k = x.shape
    if n < 3 or k < 2:
        return np.nan
    g = x.mean()
    ssr = k * ((x.mean(1) - g) ** 2).sum()
    ssc = n * ((x.mean(0) - g) ** 2).sum()
    sse = ((x - g) ** 2).sum() - ssr - ssc
    msr, msc, mse = ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)


def jaccard_matrix(pct, nv, min_cov=20, min_pct=20):
    """Pairwise Jaccard of called sites among sites both samples measure (>= min_cov reads and a
    corrected percent). Returns (Jaccard, sites called in either)."""
    k = pct.shape[1]
    cov = (np.nan_to_num(nv) >= min_cov) & ~np.isnan(pct)          # empty percent: gIVT could not correct
    called = np.nan_to_num(pct) >= min_pct
    jac, union = np.full((k, k), np.nan), np.zeros((k, k), int)
    for i in range(k):
        for j in range(k):
            both = cov[:, i] & cov[:, j]
            a, b = called[:, i] & both, called[:, j] & both
            u = int((a | b).sum())
            union[i, j] = u
            jac[i, j] = (a & b).sum() / u if u else np.nan
    return jac, union


def isolated(chroms, starts, strands, mods, min_gap=5):
    """True for sites with no site of a different modification within min_gap nt on the same
    chrom and strand (adjacent calls by two models can echo one another)."""
    starts = np.asarray(starts, float)
    mods = np.asarray(mods, dtype=object)
    ok = np.ones(len(starts), bool)
    groups = {}
    for i, k_ in enumerate(zip(chroms, strands)):
        groups.setdefault(k_, []).append(i)
    for idx in groups.values():
        idx = np.array(idx)
        order = idx[np.argsort(starts[idx], kind='stable')]
        pos, md = starts[order], mods[order]
        for j in range(len(order)):
            lo = j - 1
            while lo >= 0 and pos[j] - pos[lo] <= min_gap:
                if md[lo] != md[j]:
                    ok[order[j]] = False
                    break
                lo -= 1
            if not ok[order[j]]:
                continue
            hi = j + 1
            while hi < len(order) and pos[hi] - pos[j] <= min_gap:
                if md[hi] != md[j]:
                    ok[order[j]] = False
                    break
                hi += 1
    return ok


class Cohort:
    def __init__(self, results, ref, sheet_path, batch_col, mask, out):
        self.R, self.ref, self.mask, self.out = results, ref, mask, out
        self.tables = os.path.join(out, 'tables')
        os.makedirs(self.tables, exist_ok=True)
        self.legends = {}                                           # figure name -> legend text
        self.summary = {}                                           # sample -> {column: value}, Table 1
        self._conc = None                                           # figures 5 and S2 share one computation
        self.cohort = os.path.join(results, 'cohort', ref)
        sheet = read_table(sheet_path, ',')
        ids = text(sheet, 'sample_id')
        b = text(sheet, batch_col) if batch_col in sheet.columns else ['all'] * len(ids)
        self.sheet = sheet
        self.batch = dict(zip(ids, b))
        self.samples = sorted(ids, key=lambda s: (self.batch[s], s))
        names = sorted(set(self.batch.values()))
        self.colours = {n: BATCH_COLOURS[i % len(BATCH_COLOURS)] for i, n in enumerate(names)}
        self._mods = None

    def path(self, *parts):
        return os.path.join(self.R, *parts)

    def mods_matrix(self):
        """Masked cohort sites: (table, samples, percent, mods, n_valid or None)."""
        if self._mods is None:
            pm = read_table(os.path.join(self.cohort, f'cohort.{self.ref}.mods.pct_fpcorr.tsv.gz'))
            ms = [s for s in self.samples if s in pm.columns]
            mods = np.array(text(pm, 'mod'))
            keep = self.mask.keep(text(pm, 'chrom'), numbers(pm, ['start']).ravel(), mods)
            pct = numbers(pm, ms)[keep]
            nvp = os.path.join(self.cohort, f'cohort.{self.ref}.mods.n_valid.tsv.gz')
            nv = None
            if os.path.exists(nvp):
                nt = read_table(nvp)
                if all(text(nt, k) == text(pm, k) for k in ('chrom', 'start', 'mod', 'strand')):
                    nv = numbers(nt, ms)[keep]
                else:
                    log.warning('n_valid and pct_fpcorr matrices list different sites; figure 5 skipped')
            rows = {k: np.array(text(pm, k), dtype=object)[keep] for k in ('chrom', 'start', 'strand')}
            self._mods = (rows, ms, pct, mods[keep], nv)
        return self._mods

    def table(self, name, header, rows):
        write_rows(os.path.join(self.tables, name), header, rows)
        return name

    def note_sample(self, s, key, value):
        self.summary.setdefault(s, {})[key] = value

    def depth(self):
        qc = read_table(os.path.join(self.cohort, f'cohort.{self.ref}.read_qc.tsv'))
        other = (numbers(qc, ['non_globin_reads']) - numbers(qc, ['mito_reads'])).ravel()
        return dict(zip(text(qc, 'sample'), other))


# ---------------------------------------------------------------- figure 1: read accuracy
def fig1_reads(c, rs_dir, notes):
    batches = sorted(set(c.batch.values()))
    ident, subs = {b: [] for b in batches}, {b: np.zeros((4, 4)) for b in batches}
    found = 0
    for s in c.samples:
        ip = os.path.join(rs_dir, f'{s}.read_identity.tsv.gz')
        sp = os.path.join(rs_dir, f'{s}.substitutions.tsv')
        if os.path.exists(ip):
            t = read_table(ip)
            ident[c.batch[s]].append(numbers(t, ['identity', 'aligned_len', 'ins_rate', 'del_rate']))
            found += 1
        if os.path.exists(sp):
            t = read_table(sp)
            idx = {b: i for i, b in enumerate('ACGU')}
            for r_, b_, n_ in zip(text(t, 'ref_base'), text(t, 'read_base'), numbers(t, ['count']).ravel()):
                subs[c.batch[s]][idx[r_], idx[b_]] += n_
    if not found:
        notes.append(f'figure 1 skipped: no read_stats outputs in {rs_dir} (run scripts/read_stats.sh per sample)')
        return None
    ident = {b: np.vstack(v) for b, v in ident.items() if v}
    nb = len(ident)
    fig, gs = figure(2, max(3, nb))
    ax = fig.add_subplot(gs[0, 0])                                  # A identity
    bins = np.linspace(0.80, 1.0, 81)
    ymax = 0
    for b, v in ident.items():
        h, _, _ = ax.hist(v[:, 0], bins=bins, density=True, color=c.colours[b], alpha=0.55, label=b)
        ymax = max(ymax, h.max())
    for j, (b, v) in enumerate(ident.items()):
        med = np.median(v[:, 0])
        ax.axvline(med, ls='--', color=c.colours[b], lw=1)
        ax.text(med, ymax * (1.05 - 0.09 * j), f'{med:.4f}', color=c.colours[b], fontsize=9, ha='right', va='top')
    ax.set_xlim(0.80, 1.0)
    ax.set_xlabel('identity')
    ax.set_ylabel('density')
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B identity by aligned length
    edges = [0, 500, 1000, 1500, 2000, np.inf]
    labels = ['0-500', '500-1000', '1000-1500', '1500-2000', '2000+']
    w = 0.8 / nb
    for j, (b, v) in enumerate(ident.items()):
        data = [v[(v[:, 1] >= lo) & (v[:, 1] < hi), 0] for lo, hi in zip(edges[:-1], edges[1:])]
        pos = np.arange(len(labels)) - 0.4 + w * (j + 0.5)
        bp = ax.boxplot([d if len(d) else [np.nan] for d in data], positions=pos, widths=w * 0.9, showfliers=False,
                        patch_artist=True, medianprops={'color': 'black', 'lw': 0.8},
                        whiskerprops={'lw': 0.6}, capprops={'lw': 0.6}, boxprops={'lw': 0.6})
        for patch in bp['boxes']:
            patch.set_facecolor(c.colours[b])
    ax.set_xticks(range(len(labels)), labels, rotation=30, fontsize=9)
    ax.set_xlabel('aligned length (nt)')
    ax.set_ylabel('identity')
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C indels
    grid = np.linspace(0, 0.10, 201)
    top = 0
    for b, v in ident.items():
        indel = v[:, 2] + v[:, 3]
        h, e = np.histogram(indel, bins=grid, density=True)
        mid = (e[:-1] + e[1:]) / 2
        ax.fill_between(mid, h, color=c.colours[b], alpha=0.3, lw=0)
        ax.plot(mid, h, color=c.colours[b], lw=1)
        top = max(top, h.max())
    for j, (b, v) in enumerate(ident.items()):
        med = np.median(v[:, 2] + v[:, 3])
        ax.axvline(med, ls='--', color=c.colours[b], lw=1)
        ax.text(med + 0.001, top * (1.0 - 0.09 * j), f'{med:.4f}', color=c.colours[b], fontsize=9, va='top')
    ax.set_xlim(0, 0.10)
    ax.set_xlabel('indels per aligned base')
    ax.set_ylabel('density')
    letter(ax, 'C')

    for j, b in enumerate(batches[:max(3, nb)]):                     # D... substitution matrices
        if b not in ident:
            continue
        ax = fig.add_subplot(gs[1, j])
        m = subs[b]
        rate = 100 * m / np.maximum(m.sum(1, keepdims=True), 1)
        show = np.log10(np.clip(rate, 0.01, 100))                    # colour on a log scale, as the matrix spans 0.1-99%
        im = ax.pcolormesh(np.arange(5) - 0.5, np.arange(5) - 0.5, show, cmap=MATRIX_CMAP, vmin=-1.3, vmax=2,
                           edgecolors='face', linewidth=0.3)
        for r_ in range(4):
            for b_ in range(4):
                ax.text(b_, r_, f'{rate[r_, b_]:.2f}', ha='center', va='center', fontsize=9,
                        color='black' if r_ == b_ else 'white')
        ax.set_xlim(-0.5, 3.5)
        ax.set_ylim(3.5, -0.5)
        ax.set_aspect('equal')
        ax.set_xticks(range(4), list('ACGU'))
        ax.set_yticks(range(4), list('ACGU'))
        ax.set_xlabel(f'{b} read base')
        ax.set_ylabel('reference base')
        ax.tick_params(length=0)
        for side in ('left', 'bottom'):
            ax.spines[side].set_visible(False)
        letter(ax, 'DEFGHIJK'[j])
    fig_batch_legend(fig, {b: c.colours[b] for b in ident})
    nreads = sum(len(v) for v in ident.values())
    notes.append(f'figure 1: {nreads:,} sampled reads (>= 200 aligned bases), {found} samples')
    c.legends['Fig1_read_qc'] = f"""Read-level accuracy per batch, as the paper's Fig. 1 ({found} samples, {nreads:,} reads).
Reads are a random 2% of each sample's primary alignments with MAPQ >= 40 and >= 200 aligned bases
(scripts/read_stats.sh), pooled within each batch: every read counts once, so deeper samples weigh more.
A, identity, matches / (matches + mismatches + insertions + deletions); dashed lines and numbers: medians.
B, identity by aligned length (reference bases), boxes per batch (median, quartiles; whiskers 1.5 IQR).
C, insertions plus deletions per aligned base; dashed: medians.
D onwards, substitution matrix per batch: % of reference bases (rows, RNA sense) read as each base
(columns), pooled over positions of the read_stats.sh region (default chr22) with >= 20 reads on that strand."""
    return fig


# ---------------------------------------------------------------- figure 2: cohort overview
def fig2_overview(c, notes):
    qc = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.read_qc.tsv'))
    row = {s: i for i, s in enumerate(text(qc, 'sample'))}
    ss = [s for s in c.samples if s in row]
    v = {k: numbers(qc, [k]).ravel() for k in ('mapped_reads', 'globin_reads', 'mito_reads', 'non_globin_reads')}
    g = {k: np.array([v[k][row[s]] for s in ss]) for k in v}
    other = g['non_globin_reads'] - g['mito_reads']
    fig, gs = figure(2, 3, widths=[PANEL[0], PANEL[0], PANEL[0] * 1.6], heights=[PANEL[1], PANEL[1] * 1.6])
    x = np.arange(len(ss))
    ax = fig.add_subplot(gs[0, 0])                                  # A composition
    bottom = np.zeros(len(ss))
    for name, val, col in (('other', other, '#5769AB'), ('mitochondrial', g['mito_reads'], '#EDB157'),
                           ('globin', g['globin_reads'], '0.78')):
        ax.bar(x, 100 * val / g['mapped_reads'], bottom=bottom, color=col, label=name, width=0.8)
        bottom += 100 * val / g['mapped_reads']
    for i, t in enumerate(g['mapped_reads']):
        ax.text(i, 101, f'{t / 1e6:.0f}', ha='center', va='bottom', fontsize=7)
    ax.set_ylim(0, 112)
    ax.set_ylabel('mapped reads (%)')
    sample_axis(ax, ss, c.batch, c.colours)
    ax.legend(fontsize=8, loc='lower left', bbox_to_anchor=(0, 1.08), ncol=3, handlelength=1, columnspacing=0.8)
    ax.text(1.0, 1.0, 'M reads', transform=ax.transAxes, ha='right', va='bottom', fontsize=7)
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B non-globin depth
    ax.bar(x, other / 1e6, color=[c.colours[c.batch[s]] for s in ss], width=0.8)
    med = np.median(other / 1e6)
    ax.axhline(med, ls='--', color='black', lw=0.8)
    ax.text(len(ss) - 0.5, med, f'{med:.2f}', va='bottom', ha='right', fontsize=9)
    ax.set_ylabel('non-globin reads (M)')
    sample_axis(ax, ss, c.batch, c.colours)
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C sites and genes vs depth
    _, mss, mpct, _, _ = c.mods_matrix()
    sites = dict(zip(mss, (np.nan_to_num(mpct) >= 20).sum(0).astype(float)))   # masked, as every mod figure
    ys = np.array([sites.get(s, np.nan) for s in ss]) / 1e3
    ok = ~np.isnan(ys)
    ax.scatter(other[ok] / 1e6, ys[ok], c=[c.colours[c.batch[s]] for s, k in zip(ss, ok) if k], s=30, zorder=3)
    for s, a, b in zip(ss, other / 1e6, ys):
        if not np.isnan(b):
            ax.annotate(s[-2:], (a, b), fontsize=7, xytext=(3, 2), textcoords='offset points')
    if ok.sum() >= 3:
        k, b0 = np.polyfit(other[ok] / 1e6, ys[ok], 1)
        xx = np.linspace(0, (other[ok] / 1e6).max() * 1.05, 10)
        ax.plot(xx, k * xx + b0, ls='--', color='0.5', lw=0.8)
        stat_box(ax, [f'n = {int(ok.sum())}', f'r = {np.corrcoef(other[ok], ys[ok])[0, 1]:.3f}'])
    ax.set_xlabel('non-globin reads (M)')
    ax.set_ylabel('modification sites (k)')
    letter(ax, 'C')

    es, lx, _ = expression(c)
    ax = fig.add_subplot(gs[1, 0:2])                                # D PCA
    keep = np.argsort(lx.var(axis=1))[::-1][:min(2000, len(lx))]
    xc = lx[keep] - lx[keep].mean(axis=1, keepdims=True)
    _, sv, vt = np.linalg.svd(xc, full_matrices=False)
    pcs, frac = vt.T * sv, sv ** 2 / (sv ** 2).sum()
    for i, s in enumerate(es):
        ax.scatter(pcs[i, 0], pcs[i, 1], color=c.colours[c.batch[s]], s=45, zorder=3)
        ax.annotate(s, (pcs[i, 0], pcs[i, 1]), fontsize=8, xytext=(4, 3), textcoords='offset points')
    ax.set_xlabel(f'PC1 ({100 * frac[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({100 * frac[1]:.1f}%)')
    letter(ax, 'D')

    ax = fig.add_subplot(gs[1, 2])                                  # E correlation matrix
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
            if i != j:
                ax.text(j, i, f'{r[i, j]:.3f}'.lstrip('0'), ha='center', va='center', fontsize=5.5,
                        color='black' if r[i, j] > (np.nanmin(off) + np.nanmax(off)) / 2 else 'white')
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect('equal')
    sample_axis(ax, es, c.batch, c.colours, fontsize=7)
    ax.set_yticks(range(n), es, fontsize=7)
    for lab, s in zip(ax.get_yticklabels(), es):
        lab.set_color(c.colours[c.batch[s]])
    ax.tick_params(length=0)
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(False)
    colourbar(fig, im, ax, label='Pearson r, log$_2$(CPM + 1)')
    letter(ax, 'E')
    fig_batch_legend(fig, c.colours)
    notes.append(f'figure 2: {lx.shape[0]:,} genes with median CPM > 1')
    for i, s in enumerate(ss):
        c.note_sample(s, 'mapped_reads', int(g['mapped_reads'][i]))
        c.note_sample(s, 'globin_pct', round(100 * g['globin_reads'][i] / g['mapped_reads'][i], 1))
        c.note_sample(s, 'mito_pct', round(100 * g['mito_reads'][i] / g['mapped_reads'][i], 1))
        c.note_sample(s, 'non_globin_non_mito_reads', int(other[i]))
    gcp = os.path.join(c.cohort, f'cohort.{c.ref}.gene_counts.tsv.gz')
    if os.path.exists(gcp):
        gct = read_table(gcp)
        gs_ = [s for s in ss if s in gct.columns]
        for s, d in zip(gs_, (np.nan_to_num(numbers(gct, gs_)) >= 10).sum(0)):
            c.note_sample(s, 'genes_ge10_reads', int(d))
    c.legends['Fig2_cohort_overview'] = f"""Cohort overview, one value per sample throughout ({len(ss)} samples, coloured by batch).
A, mapped reads in globin genes, on chrM, and everything else, as % of the sample's mapped reads; the number
above each bar is the sample's mapped reads (millions). B, mapped reads outside globin genes and chrM: the
depth that the rest of the analysis works with; dashed: cohort median. C, modification sites per sample (all
eight modifications; >= 20 reads and >= 20% after gIVT correction in that sample; masked regions left out)
against B; dashed: least-squares line, r = Pearson. D, PCA of log2(CPM + 1) over the 2,000 most variable of
the {lx.shape[0]:,} genes with median CPM > 1. E, Pearson r between samples over the same {lx.shape[0]:,} genes.
Per-sample numbers: tables/Table1_cohort_summary.tsv."""
    return fig


# ---------------------------------------------------------------- figure 3: gIVT correction
# Each <sample>.cov20_fpcorr.bed has a row for every position x mod code with >= 20 reads: tens of
# millions of rows, too many to load. awk streams each one into 2%-wide bins of Dorado percent
# (column 11) x corrected percent (column 19; -1 where gIVT could not correct), applying the same
# masks as the other modification figures. Bins are cached in <out>/.cache, keyed by the file's
# size and time and by the masks, so a rerun only reads files that changed.
IVT_AWK = r"""
BEGIN { FS = "\t"; OFS = "\t" }
FNR == NR { if ($1 != "__none__") { k = ++nm[$1]; ms[$1, k] = $2 + 0; me[$1, k] = $3 + 0 }; next }
NF < 19 { short++; next }
{
    if (rx != "" && $1 ~ rx) { mc[$4]++; next }
    if ($1 in nm) { p = $2 + 0; for (i = 1; i <= nm[$1]; i++) if (p >= ms[$1, i] && p < me[$1, i]) { mb[$4]++; next } }
    bx = int($11 / 2); if (bx > 49) bx = 49; if (bx < 0) bx = 0
    if ($19 == "") by = -1; else { by = int($19 / 2); if (by > 49) by = 49; if (by < 0) by = 0 }
    h[$4 OFS bx OFS by]++
}
END {
    for (k in h) print "bin", k, h[k]
    for (m in mc) print "contig", m, "", "", mc[m]
    for (m in mb) print "bed", m, "", "", mb[m]
    if (short) print "short", "", "", "", short
}
"""


def ivt_bins(c, out, threads, notes):
    """sample -> list of (kind, mod, bx, by, n) rows from each sample's cov20_fpcorr.bed."""
    awk = shutil.which('mawk') or shutil.which('gawk') or shutil.which('awk')
    if not awk:
        sys.exit('cohort_plots.py: figure 3 needs awk on PATH')
    cache = os.path.join(out, '.cache')
    os.makedirs(cache, exist_ok=True)
    mask_lines = ['__none__\t0\t0'] + [f'{c_}\t{a}\t{b}' for c_, a, b in c.mask.intervals()]
    mask_text = '\n'.join(mask_lines) + '\n'
    mask_path = os.path.join(cache, 'mask_intervals.tsv')
    with open(mask_path, 'w') as fh:
        fh.write(mask_text)
    sig = hashlib.md5((c.mask.pattern + '\n' + mask_text + IVT_AWK).encode()).hexdigest()
    todo, paths = [], {}
    for s in c.samples:
        p = c.path('mods', c.ref, s, f'{s}.cov20_fpcorr.bed')
        if not os.path.exists(p):
            continue
        st = os.stat(p)
        key = f'#{os.path.abspath(p)}\t{st.st_size}\t{int(st.st_mtime)}\t{sig}'
        cp = os.path.join(cache, f'{s}.ivt_bins.tsv')
        paths[s] = cp
        if os.path.exists(cp):
            with open(cp) as fh:
                if fh.readline().rstrip('\n') == key:
                    continue
        todo.append((s, p, cp, key))

    def run(job):
        s, p, cp, key = job
        fd, tmp = tempfile.mkstemp(dir=cache, prefix=f'{s}.', suffix='.partial')
        with os.fdopen(fd, 'w') as fh:
            fh.write(key + '\n')
            fh.flush()
            r = subprocess.run([awk, '-v', f'rx={c.mask.pattern}', IVT_AWK, mask_path, p], stdout=fh,
                               stderr=subprocess.PIPE, text=True)
        if r.returncode:
            os.unlink(tmp)
            return s, r.stderr.strip() or f'awk exit {r.returncode}'
        os.chmod(tmp, 0o644)
        os.replace(tmp, cp)
        return s, None

    if todo:
        log.info('figure 3: binning %d cov20_fpcorr.bed file(s) with %s, %d at a time', len(todo), os.path.basename(awk), threads)
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            for s, err in pool.map(run, todo):
                if err:
                    sys.exit(f'cohort_plots.py: binning {s}.cov20_fpcorr.bed failed: {err}')
    bins = {}
    for s, cp in paths.items():
        rows = []
        with open(cp) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip('\n').split('\t')
                rows.append((f[0], f[1], int(f[2]) if f[2] else 0, int(f[3]) if f[3] else 0, int(f[4])))
        bins[s] = rows
    return bins


def fig3_ivt(c, out, threads, notes):
    bins = ivt_bins(c, out, threads, notes)
    if not bins:
        notes.append('figure 3 skipped: no <sample>.cov20_fpcorr.bed in results/mods')
        return None
    hists = {code: np.zeros((50, 50)) for code in MOD_NAME}
    table, masked, short = [], {}, 0
    for s in c.samples:
        if s not in bins:
            continue
        per = {}
        for kind, code, bx, by, n in bins[s]:
            if kind == 'bin':
                per.setdefault(code, []).append((bx, by, n))
                if by >= 0 and code in hists:
                    hists[code][bx, by] += n
            elif kind in ('contig', 'bed'):
                masked[(kind, code)] = masked.get((kind, code), 0) + n
            elif kind == 'short':
                short += n
        for code in MAIN_MODS + OME_MODS:
            if code not in per:
                continue
            a = np.array(per[code])
            bx, by, n = a[:, 0], a[:, 1], a[:, 2]
            r20, c20, unc = bx >= 10, by >= 10, by < 0
            table.append((s, plain(code), int(n.sum()), int(n[r20].sum()), int(n[c20].sum()), int(n[r20 & c20].sum()),
                          int(n[r20 & ~c20 & ~unc].sum()), int(n[~r20 & c20].sum()), int(n[r20 & unc].sum())))
    tname = c.table('Fig3_ivt_correction_sites.tsv',
                    ['sample', 'mod', 'rows_cov20', 'dorado_ge20', 'corrected_ge20', 'kept', 'removed', 'raised',
                     'not_correctable_ge20'], table)

    fig, gs = figure(3, 4, heights=[PANEL[1] * 1.15, PANEL[1] * 1.15, PANEL[1]])
    edges = np.linspace(0, 100, 51)
    vmax = max(10, max(h.max() for h in hists.values()))
    axes, last = [], None
    box = dict(boxstyle='square,pad=0.2', fc='white', ec='none', alpha=0.85)
    for k, code in enumerate(MAIN_MODS + OME_MODS):
        ax = fig.add_subplot(gs[k // 4, k % 4])
        axes.append(ax)
        h = np.ma.masked_equal(hists[code].T, 0)
        if h.count():
            last = ax.pcolormesh(edges, edges, h, cmap='viridis', norm=LogNorm(vmin=1, vmax=vmax), edgecolors='face', linewidth=0.2)
        ax.axvline(20, color='#E00000', ls='--', lw=0.9)
        ax.axhline(20, color='#E00000', ls='--', lw=0.9)
        hh = hists[code]
        kept, removed, raised = hh[10:, 10:].sum(), hh[10:, :10].sum(), hh[:10, 10:].sum()
        ax.text(23, 97, f'kept {kept:,.0f}', ha='left', va='top', fontsize=8, bbox=box)
        ax.text(98, 2, f'removed {removed:,.0f}', ha='right', va='bottom', fontsize=8, bbox=box)
        if raised:
            ax.text(2, 98, f'raised {raised:,.0f}', ha='left', va='top', fontsize=8, bbox=box)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_aspect('equal')
        ax.set_xlabel('Dorado occupancy (%)')
        ax.set_ylabel('after gIVT (%)')
        ax.set_title(MOD_NAME[code], color=MOD_COLOUR[code], fontweight='bold')
        letter(ax, 'ABCDEFGH'[k])
    if last is not None:
        colourbar(fig, last, axes, label='site rows (all samples)')

    ax = fig.add_subplot(gs[2, :])                                  # I sites before and after, per sample
    codes = [m for m in MAIN_MODS + OME_MODS if any(r[1] == plain(m) for r in table)]
    w = 0.38
    for i, code in enumerate(codes):
        before = np.array([r[3] for r in table if r[1] == plain(code)])
        after = np.array([r[4] for r in table if r[1] == plain(code)])
        for off, vals, alpha in ((-w / 2, before, 0.45), (w / 2, after, 1.0)):
            ax.bar(i + off, max(np.median(vals), 0.8), width=w, color=MOD_COLOUR[code], alpha=alpha)
            ax.vlines(i + off, max(vals.min(), 0.8), max(vals.max(), 0.8), color='black', lw=0.6)
        drop = 100 * (1 - np.median(after) / max(np.median(before), 1))
        ax.text(i, max(np.median(before), 1) * 1.6, f'−{drop:.0f}%', ha='center', fontsize=9)
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    ax.set_xticks(range(len(codes)), [MOD_NAME[m] for m in codes])
    ax.set_ylabel('sites per sample (≥ 20%)')
    ax.set_title('before (light) and after (dark) gIVT subtraction: median and range over samples', fontsize=11)
    letter(ax, 'I')
    unc = sum(r[8] for r in table)
    notes.append(f'figure 3: {len(bins)} samples; {unc:,} site rows >= 20% that gIVT could not correct are left out of A-H '
                 f'(and out of every final site set); per-sample counts in tables/{tname}')
    c.legends['Fig3_ivt_correction'] = f"""Effect of gIVT false-positive subtraction, as the paper's Fig. 3 and Table 4 ({len(bins)} samples).
A-H, POOLED over samples: every position x modification with >= 20 reads in a sample is one entry, so a
position covered in all {len(bins)} samples appears {len(bins)} times and deeper samples contribute more. x: Dorado
occupancy (modkit percent modified: modified calls / valid calls, after the fixed per-modification confidence
thresholds); y: the same after subtracting the gIVT error rate of its 9-mer. 2%-wide bins, log colour scale shared by all panels; red dashed lines at 20%. kept: >= 20% before
and after; removed: >= 20% before, < 20% after. {unc:,} entries >= 20% that gIVT could not correct (no 9-mer,
base mismatch, 9-mer not in the table) are not drawn; they are not in any final site set either.
I, PER SAMPLE: sites >= 20% before (light) and after (dark) correction; bars are the median over samples, lines
the range; the number is the change in the median. Per-sample counts for every modification: tables/{tname}."""
    for kind in ('contig', 'bed'):
        parts = [f'{plain(m)} {masked[(kind, m)]:,}' for m in MAIN_MODS + OME_MODS if (kind, m) in masked]
        if parts:
            notes.append(f'figure 3, masked by {"contig" if kind == "contig" else "--mask-bed"}: site rows ' + ', '.join(parts))
    if short:
        notes.append(f'figure 3: {short:,} rows with fewer than 19 columns skipped')
    return fig


def plain(code):
    """ASCII modification name for tables and notes (m6A, m5C, Psi, Ino, Am, Cm, Um, Gm)."""
    return MOD_NAME.get(code, code).replace('$', '').replace('^', '').replace('Ψ', 'Psi')


# ---------------------------------------------------------------- figure 4: modification calls
def gene_sites(c):
    """mod code -> {gene_id: set of (chrom, start, strand)}, union over samples, masked."""
    out = {}
    for s in c.samples:
        p = c.path('mods', c.ref, s, f'{s}.gene_mods.bed')
        if not os.path.exists(p):
            continue
        t = read_table(p)
        if t.height == 0:
            continue
        ch, st, codes = text(t, 'chrom'), numbers(t, ['start']).ravel(), np.array(text(t, 'mod_code'))
        keep = c.mask.keep(ch, st, None)
        for code, g, c_, s_, sd, k in zip(codes, text(t, 'gene_id'), ch, st, text(t, 'strand'), keep):
            if k:
                out.setdefault(code, {}).setdefault(g, set()).add((c_, int(s_), sd))
    return out


def fig4_mods(c, notes):
    rows, ms, pct, mods, _ = c.mods_matrix()
    fig, gs = figure(2, 3, widths=[PANEL[0] * 1.15, PANEL[0] * 1.15, PANEL[0] * 1.3])
    bin_rows, co_rows = [], []
    for k, (group, lab) in enumerate(((MAIN_MODS, 'A'), (OME_MODS, 'B'))):
        ax = fig.add_subplot(gs[0, k])
        present = [m for m in group if (mods == m).any()]
        w = 0.8 / max(1, len(present))
        for j, code in enumerate(present):
            v = pct[mods == code]
            counts = np.array([[np.sum((v[:, si] >= b) & ((v[:, si] < b + 10) if b < 90 else (v[:, si] <= 100)))
                                for si in range(len(ms))] for b in OCC_BINS])
            bin_rows.extend((s_, plain(code), f'{b}-{b + 10}', int(counts[bi, si]))
                            for bi, b in enumerate(OCC_BINS) for si, s_ in enumerate(ms))
            xpos = np.arange(len(OCC_BINS)) - 0.4 + w * (j + 0.5)
            ax.bar(xpos, np.maximum(np.median(counts, 1), 0.8), width=w, color=MOD_COLOUR[code], label=MOD_NAME[code])
            ax.vlines(xpos, np.maximum(counts.min(1), 0.8), np.maximum(counts.max(1), 0.8), color='black', lw=0.6)
        ax.set_yscale('log')
        log_axis(ax.yaxis)
        ax.set_xticks(range(len(OCC_BINS)), [str(b) for b in OCC_BINS])
        ax.set_xlabel('modification occupancy (%, bin start)')
        ax.set_ylabel('sites per sample')
        ax.legend(fontsize=9, ncol=2, loc='upper right', handlelength=1)
        letter(ax, lab)

    ax = fig.add_subplot(gs[0, 2])                                  # C samples calling each site
    kk = np.arange(1, len(ms) + 1)
    for code in [m for m in MAIN_MODS + OME_MODS if (mods == m).any()]:
        n = (pct[mods == code] >= 20).sum(1)
        n = n[n > 0]
        ax.plot(kk, [100 * (n == i).mean() if len(n) else 0 for i in kk], marker='o', ms=3.5, lw=1,
                color=MOD_COLOUR[code], ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_xlabel('samples calling the site')
    ax.set_ylabel('sites (%)')
    ax.set_xticks(kk[::2] if len(kk) > 8 else kk)
    ax.set_title("colours as in A, B; 2'-O-Me dashed", fontsize=9)
    letter(ax, 'C')

    genes = gene_sites(c)
    ax = fig.add_subplot(gs[1, 0:2])                                # D sites per gene
    cap = 50
    present = [m for m in MAIN_MODS if genes.get(m)]
    width = cap + 12
    for j, code in enumerate(present):
        per_gene = np.array([len(v) for v in genes[code].values()])
        h = np.bincount(np.minimum(per_gene, cap + 1), minlength=cap + 2)[1:]
        ax.bar(j * width + np.arange(1, cap + 2), np.maximum(h, 0.8), width=1.0, color=MOD_COLOUR[code], label=MOD_NAME[code])
        ax.axvline(j * width + 0.5, color='black', lw=0.8)
    if present:
        ax.set_yscale('log')
        log_axis(ax.yaxis)
        ax.set_xticks([j * width + t for j in range(len(present)) for t in (1, 25, cap + 1)],
                      [str(t) if t <= cap else f'>{cap}' for _ in present for t in (1, 25, cap + 1)], fontsize=9)
        ax.legend(fontsize=9, loc='upper right', handlelength=1)
    ax.set_xlabel('unique modified sites per gene (≥ 20 reads, ≥ 20%)')
    ax.set_ylabel('genes')
    letter(ax, 'D')

    ax = fig.add_subplot(gs[1, 2])                                  # E co-occurrence, 5-nt rule
    allsites = [(m, g, s_) for m in MAIN_MODS if m in genes for g, ss in genes[m].items() for s_ in ss]
    if allsites:
        iso = isolated([s_[0] for _, _, s_ in allsites], [s_[1] for _, _, s_ in allsites], [s_[2] for _, _, s_ in allsites],
                       [m for m, _, _ in allsites])
        gsets = {}
        for (m, g, _), ok in zip(allsites, iso):
            if ok:
                gsets.setdefault(m, set()).add(g)
        present = [m for m in MAIN_MODS if gsets.get(m)]
        if len(present) >= 2:
            w = 0.8 / (len(present) - 1)
            for i, obs in enumerate(present):
                for j, giv in enumerate([g for g in present if g != obs]):
                    both = len(gsets[obs] & gsets[giv])
                    co_rows.append((plain(obs), plain(giv), len(gsets[giv]), both, both / len(gsets[giv])))
                    xpos = i - 0.4 + w * (j + 0.5)
                    ax.bar(xpos, both / len(gsets[giv]), width=w, color=MOD_COLOUR[giv])
                    ax.text(xpos, 0.03, f'{both:,}', rotation=90, ha='center', va='bottom', fontsize=7, color='white')
            ax.legend([plt.Rectangle((0, 0), 1, 1, color=MOD_COLOUR[m]) for m in present], [MOD_NAME[m] for m in present],
                      title='given', fontsize=9, title_fontsize=10, loc='upper left', bbox_to_anchor=(1.0, 1.0), handlelength=1)
            ax.set_xticks(range(len(present)), [MOD_NAME[m] for m in present])
            ax.set_ylim(0, 1.0)
        notes.append(f'figure 4E: {int(iso.sum()):,} of {len(iso):,} gene sites have no other modification within 5 nt')
    ax.set_xlabel('observed modification')
    ax.set_ylabel('P(observed | given)')
    ax.set_title('no other modification within 5 nt (genes in bars)', fontsize=9)
    letter(ax, 'E')
    fig_batch_legend(fig, c.colours)
    t1 = c.table('Fig4_sites_per_occupancy_bin.tsv', ['sample', 'mod', 'occupancy_bin', 'sites'], bin_rows)
    t2 = c.table('Fig4_cooccurrence.tsv', ['observed', 'given', 'genes_with_given', 'genes_with_both', 'p_observed_given'], co_rows)
    ngenes = len({g for m in genes for g in genes[m]})
    c.legends['Fig4_modification_calls'] = f"""Modification calls in the cohort, as the paper's Fig. 5A, B, D ({len(ms)} samples, {pct.shape[0]:,} cohort sites).
A, B, PER SAMPLE: sites with >= 20 reads and >= 20% gIVT-corrected occupancy in that sample, in 10% occupancy
bins (bin start shown; 90 is 90-100%); bars: median over samples, lines: range (A, m6A, m5C, Ψ, inosine;
B, 2'-O-methylations). Per-sample counts: tables/{t1}.
C, COHORT: each cohort site (called in at least one sample) by how many samples call it; % of that
modification's sites. A site measured in few samples can only be called in few, so this mixes depth and
reproducibility; figure 5 separates them.
D, POOLED (union over samples): gene-level sites from each sample's gene_mods.bed, each site counted once per
gene however many samples call it; genes by number of sites (> 50 in the last bar), {ngenes:,} genes.
E, POOLED: of the genes carrying a site of the given modification (bar colour), the fraction that also carry the
observed modification (x); a site counts only if no site of another modification lies within 5 nt, so adjacent
calls by two models cannot make co-occurrence. Genes with both in the bars; tables/{t2}."""
    return fig


# ---------------------------------------------------------------- figure 5 (ICC) and S2 (Jaccard): concordance
def pairwise_icc(pct, nv, min_cov=20, min_pct=20, min_sites=10):
    """ICC(2,1) for every pair of samples, over the sites both measure (>= min_cov reads and a corrected
    value in each) that are >= min_pct in either. Returns (ICC, sites used); NaN below min_sites."""
    k = pct.shape[1]
    cov = (np.nan_to_num(nv) >= min_cov) & ~np.isnan(pct)
    called = np.nan_to_num(pct) >= min_pct
    out, used = np.full((k, k), np.nan), np.zeros((k, k), int)
    for i in range(k):
        for j in range(i + 1, k):
            sel = cov[:, i] & cov[:, j] & (called[:, i] | called[:, j])
            used[i, j] = used[j, i] = int(sel.sum())
            if sel.sum() >= min_sites:
                out[i, j] = out[j, i] = icc21(pct[sel][:, [i, j]])
    return out, used


def concordance(c):
    """Pairwise ICC and Jaccard per modification, cohort ICC, and each pair's smaller depth; computed once
    for figures 5 and S2."""
    if getattr(c, '_conc', None) is not None:
        return c._conc
    rows, ms, pct, mods, nv = c.mods_matrix()
    if nv is None:
        c._conc = {}
        return c._conc
    n = len(ms)
    groups = np.array([c.batch[s] for s in ms])
    depth = c.depth()
    d = np.array([depth.get(s, np.nan) for s in ms]) / 1e6
    out = {'ms': ms, 'n': n, 'groups': groups, 'same': groups[:, None] == groups[None, :], 'offd': ~np.eye(n, dtype=bool),
           'mindepth': np.minimum(d[:, None], d[None, :]), 'icc': {}, 'icc_n': {}, 'jac': {}, 'jac_n': {}, 'icc_all': {}}
    for code in MAIN_MODS + OME_MODS:
        sel = mods == code
        if sel.sum() < 3:
            continue
        out['icc'][code], out['icc_n'][code] = pairwise_icc(pct[sel], nv[sel])
        out['jac'][code], out['jac_n'][code] = jaccard_matrix(pct[sel], nv[sel])
        full = ((np.nan_to_num(nv[sel]) >= 20) & ~np.isnan(pct[sel])).all(1) & (np.nan_to_num(pct[sel]) >= 20).any(1)
        out['icc_all'][code] = (icc21(pct[sel][full]) if full.sum() >= 3 else np.nan, int(full.sum()))
    summary, pairs = [], []
    same, offd = out['same'], out['offd']
    for code in out['icc']:
        im, jm = out['icc'][code], out['jac'][code]
        w = lambda m, msk: float(np.nanmean(m[msk])) if np.isfinite(m[msk]).any() else np.nan
        summary.append((plain(code), out['icc_all'][code][1], out['icc_all'][code][0], w(im, offd), w(im, same & offd),
                        w(im, ~same), w(jm, offd), w(jm, same & offd), w(jm, ~same)))
        for a in range(n):
            for b in range(a + 1, n):
                pairs.append((plain(code), ms[a], ms[b], int(same[a, b]), out['mindepth'][a, b], im[a, b],
                              int(out['icc_n'][code][a, b]), jm[a, b], int(out['jac_n'][code][a, b])))
    out['t1'] = c.table('Fig5_replicate_concordance.tsv',
                        ['mod', 'icc_sites_all_samples', 'icc21_all_samples', 'mean_pairwise_icc21', 'icc21_within_batch',
                         'icc21_between_batches', 'mean_jaccard', 'jaccard_within_batch', 'jaccard_between_batches'], summary)
    out['t2'] = c.table('Fig5_pairwise.tsv', ['mod', 'sample_a', 'sample_b', 'same_batch', 'smaller_depth_M', 'icc21',
                                              'icc21_sites', 'jaccard', 'jaccard_sites_called_in_either'], pairs)
    c._conc = out
    return out


def pair_heatmaps(fig, gs, c, conc, mats, title_fmt):
    """Eight sample x sample heatmaps (one per modification), each on its own colour scale, batches outlined."""
    ms, n, groups = conc['ms'], conc['n'], conc['groups']
    cuts = [i + 0.5 for i in range(n - 1) if groups[i] != groups[i + 1]]
    cmap = MATRIX_CMAP.copy()
    cmap.set_bad('0.82')
    for k, code in enumerate(MAIN_MODS + OME_MODS):
        ax = fig.add_subplot(gs[k // 4, k % 4])
        if code in mats:
            m = mats[code].copy()
            np.fill_diagonal(m, np.nan)
            off = m[np.isfinite(m)]
            lo, hi = (float(np.percentile(off, 2)), float(np.percentile(off, 98))) if len(off) else (0.0, 1.0)
            mesh = ax.pcolormesh(np.arange(n + 1) - 0.5, np.arange(n + 1) - 0.5, np.ma.masked_invalid(m), cmap=cmap,
                                 vmin=lo, vmax=max(hi, lo + 1e-3), edgecolors='face', linewidth=0.2)
            for e in cuts:
                ax.axhline(e, color='black', lw=0.8)
                ax.axvline(e, color='black', lw=0.8)
            cb = colourbar(fig, mesh, ax, label='')
            cb.ax.tick_params(labelsize=7)
            cb.locator = plt.MaxNLocator(4)
            cb.update_ticks()
            ax.set_title(title_fmt.format(name=MOD_NAME[code], mean=np.nanmean(m)), color=MOD_COLOUR[code],
                         fontsize=12, fontweight='bold')
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_aspect('equal')
        ax.set_xticks(range(n), [s[-2:] for s in ms], fontsize=6, rotation=90)
        ax.set_yticks(range(n), [s[-2:] for s in ms], fontsize=6)
        for labs in (ax.get_xticklabels(), ax.get_yticklabels()):
            for lab, s in zip(labs, ms):
                lab.set_color(c.colours[c.batch[s]])
        ax.tick_params(length=0)
        for side in ('left', 'bottom'):
            ax.spines[side].set_visible(False)
        letter(ax, 'ABCDEFGH'[k])


def within_between(ax, conc, mats, ylabel):
    same, offd = conc['same'], conc['offd']
    codes = [m for m in MAIN_MODS + OME_MODS if m in mats]
    for i, code in enumerate(codes):
        m = mats[code]
        wv = np.nanmean(m[same & offd]) if np.isfinite(m[same & offd]).any() else np.nan
        bv = np.nanmean(m[~same]) if np.isfinite(m[~same]).any() else np.nan
        ax.plot([i - 0.15, i + 0.15], [wv, bv], color='0.6', lw=0.8, zorder=1)
        ax.scatter([i - 0.15], [wv], color=MOD_COLOUR[code], marker='o', s=36, zorder=3)
        ax.scatter([i + 0.15], [bv], facecolor='white', edgecolor=MOD_COLOUR[code], marker='o', s=36, zorder=3, lw=1.2)
    ax.set_xticks(range(len(codes)), [MOD_NAME[m] for m in codes], fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_title('same batch (filled), different batches (open)', fontsize=9)


def depth_check(ax, c, conc, mats, ylabel):
    """One point per sample pair: its agreement averaged over the modifications, against the pair's smaller
    depth. Same-batch pairs filled in the batch colour, different-batch pairs open."""
    ms, n, same = conc['ms'], conc['n'], conc['same']
    stack = np.array([mats[m] for m in mats])
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(stack, 0)
    xs, ys, sm = [], [], []
    for a in range(n):
        for b in range(a + 1, n):
            if np.isfinite(mean[a, b]) and np.isfinite(conc['mindepth'][a, b]):
                xs.append(conc['mindepth'][a, b])
                ys.append(mean[a, b])
                sm.append(same[a, b])
                if same[a, b]:
                    ax.scatter(xs[-1], ys[-1], color=c.colours[c.batch[ms[a]]], s=30, zorder=3)
                else:
                    ax.scatter(xs[-1], ys[-1], facecolor='white', edgecolor='0.45', s=22, lw=0.8, zorder=2)
    if len(xs) >= 3:
        rho = stats.spearmanr(xs, ys)[0]
        stat_box(ax, [f'{len(xs)} pairs', f'Spearman ρ = {rho:.2f}'], loc='upper left' if rho > 0 else 'upper right')
    ax.set_xlabel('smaller depth of the pair (M non-globin reads)')
    ax.set_ylabel(ylabel)
    ax.set_title('each pair of samples; filled: same batch (its colour)', fontsize=9)


def fig5_concordance(c, notes):
    conc = concordance(c)
    if not conc:
        notes.append('figure 5 skipped: needs cohort.<ref>.mods.n_valid.tsv.gz (coverage per site)')
        return None
    n = conc['n']
    fig, gs = figure(3, 4, widths=[PANEL[0] * 1.08] * 4, heights=[PANEL[1] * 1.25, PANEL[1] * 1.25, PANEL[1] * 1.1])
    pair_heatmaps(fig, gs, c, conc, conc['icc'], '{name}   mean {mean:.2f}')

    ax = fig.add_subplot(gs[2, 0])                                  # I cohort ICC, all samples
    ia = conc['icc_all']
    codes = [m for m in MAIN_MODS + OME_MODS if m in ia and not np.isnan(ia[m][0])]
    ax.bar(range(len(codes)), [ia[m][0] for m in codes], color=[MOD_COLOUR[m] for m in codes], width=0.7)
    for i, m in enumerate(codes):
        ax.text(i, ia[m][0] + 0.01, f'{ia[m][0]:.2f}', ha='center', va='bottom', fontsize=7)
        ax.text(i, 0.02, f'{ia[m][1]:,}', rotation=90, ha='center', va='bottom', fontsize=6, color='white')
    ax.set_xticks(range(len(codes)), [MOD_NAME[m] for m in codes], fontsize=8)
    if not codes:
        ax.text(0.5, 0.5, f'fewer than 3 sites measured\nin all {n} samples', transform=ax.transAxes, ha='center', va='center')
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('ICC(2,1), all samples')
    ax.set_title(f'sites measured in all {n} (in bars)', fontsize=9)
    letter(ax, 'I')

    ax = fig.add_subplot(gs[2, 1])                                  # J within against between batches
    within_between(ax, conc, conc['icc'], 'mean pairwise ICC(2,1)')
    letter(ax, 'J')
    ax = fig.add_subplot(gs[2, 2:])                                 # K depth or batch
    depth_check(ax, c, conc, conc['icc'], 'pairwise ICC(2,1), mean over mods')
    letter(ax, 'K')
    fig_batch_legend(fig, c.colours)
    c.legends['Fig5_replicate_concordance'] = f"""Agreement of modification occupancy between samples ({n} samples), as the paper's HEK293 replicates.
A-H, ICC(2,1) for every PAIR of samples (two-way random effects, absolute agreement, single measurement: a
sample that reads systematically higher counts as disagreement, so a batch offset lowers it). Each pair uses
the cohort sites both samples measure (>= 20 reads and a gIVT-corrected value in each) and that are >= 20% in
either; pairs with fewer than 10 such sites are grey, as is the diagonal. Each panel has its own colour scale
(2nd to 98th percentile of its pairs); compare modifications by the mean in the title, remembering that ICC also
rises with how widely occupancy varies between sites. Black lines separate batches: a batch effect would show as
lighter blocks on the diagonal.
I, ICC(2,1) over all {n} samples at once, at the sites measured in every sample and >= 20% in at least one (the
paper's statistic); the number of sites in the bars. The shallowest sample decides which sites qualify.
J, mean pairwise ICC over pairs from the same batch (filled) and from different batches (open).
K, depth or batch: one point per pair, its ICC averaged over the eight modifications, against the smaller of the
two samples' depths (reads outside globin genes and chrM). Filled, in the batch colour: pairs within a batch;
open: pairs across batches. If within-batch pairs sit above between-batch pairs of the same depth, that is a
batch effect; if all pairs follow depth together, it is depth.
Call-level agreement (Jaccard of the sites called >= 20%) is figure S2. Tables: tables/{conc['t1']} (per
modification), tables/{conc['t2']} (every pair: depth, ICC, Jaccard and the sites each used)."""
    notes.append(f'figure 5: pairwise ICC(2,1); tables {conc["t1"]}, {conc["t2"]}')
    return fig


def figS2_jaccard(c, notes):
    conc = concordance(c)
    if not conc:
        notes.append('figure S2 skipped: needs cohort.<ref>.mods.n_valid.tsv.gz (coverage per site)')
        return None
    n = conc['n']
    fig, gs = figure(3, 4, widths=[PANEL[0] * 1.08] * 4, heights=[PANEL[1] * 1.25, PANEL[1] * 1.25, PANEL[1] * 1.1])
    pair_heatmaps(fig, gs, c, conc, conc['jac'], '{name}   mean {mean:.2f}')
    ax = fig.add_subplot(gs[2, 0:2])
    within_between(ax, conc, conc['jac'], 'mean Jaccard index')
    letter(ax, 'I')
    ax = fig.add_subplot(gs[2, 2:])
    depth_check(ax, c, conc, conc['jac'], 'Jaccard, mean over mods')
    letter(ax, 'J')
    fig_batch_legend(fig, c.colours)
    c.legends['FigS2_jaccard'] = f"""Agreement of modification CALLS between samples ({n} samples): do two samples call the same sites?
A-H, Jaccard index for every pair of samples: among the cohort sites both samples measure (>= 20 reads and a
gIVT-corrected value in each), the sites called (>= 20%) in both divided by the sites called in either. Each
panel has its own colour scale; black lines separate batches; grey: the sample itself. Jaccard is lower than
ICC (figure 5) mostly because of the cutoff: most sites lie at 20-30% (figure 4A), where read sampling alone moves
a site across 20% from one sample to the next. Modifications with more high-occupancy sites (m6A) score higher.
I, mean Jaccard over pairs from the same batch (filled) and different batches (open).
J, as figure 5K, with Jaccard: one point per pair, mean over the modifications, against the pair's smaller depth.
Numbers: tables/{conc['t2']}."""
    return fig


# ---------------------------------------------------------------- figure 6: isoform-specific modifications
def fig6_isoforms(c, notes, min_delta=20):
    per_sample, deltas, site_sum = {}, {}, {}
    gene_name = {}
    found = 0
    for s in c.samples:
        p = c.path('mods', c.ref, s, f'{s}.feature_mods.tsv.gz')
        if not os.path.exists(p):
            continue
        found += 1
        t = read_table(p, columns=['level', 'feature', 'gene_id', 'ref_gene_name', 'chrom', 'pos', 'strand', 'mod',
                                   'n_valid', 'percent_fpcorr'])
        lv, gid, gnm = text(t, 'level'), text(t, 'gene_id'), text(t, 'ref_gene_name')
        chs, sds, mcs = text(t, 'chrom'), text(t, 'strand'), text(t, 'mod')
        for l_, g_, n_ in zip(lv, gid, gnm):
            if l_ == 'gene' and n_ not in ('.', None):
                gene_name[g_] = n_
        v = numbers(t, ['pos', 'n_valid', 'percent_fpcorr'])
        iso = np.array([l_ == 'isoform' for l_ in lv]) & (np.nan_to_num(v[:, 1]) >= 20) & ~np.isnan(v[:, 2])
        keep = c.mask.keep(chs, v[:, 0], None)
        groups = {}
        for i in np.where(iso & keep)[0]:
            key = (gid[i], chs[i], int(v[i, 0]), sds[i], mcs[i])
            groups.setdefault(key, []).append(v[i, 2])
        n_hits = {}
        for key, occ in groups.items():
            if len(occ) < 2:
                continue
            d = max(occ) - min(occ)
            deltas.setdefault(key[4], []).append(d)
            site_sum.setdefault(key, []).append(d)
            if d >= min_delta:
                n_hits[key[4]] = n_hits.get(key[4], 0) + 1
        per_sample[s] = n_hits
    if not found:
        notes.append('figure 6 skipped: no <sample>.feature_mods.tsv.gz')
        return None
    fig, gs = figure(2, 2, widths=[PANEL[0] * 1.3, PANEL[0] * 1.6])
    ax = fig.add_subplot(gs[0, 0])                                  # A per sample
    ss = [s for s in c.samples if s in per_sample]
    for code in [m for m in MAIN_MODS + OME_MODS if m in deltas]:
        ax.plot(range(len(ss)), [per_sample[s].get(code, 0) for s in ss], marker='o', ms=4, lw=1, color=MOD_COLOUR[code],
                ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylabel(f'sites with isoform Δ ≥ {min_delta}%')
    sample_axis(ax, ss, c.batch, c.colours)
    ax.legend(fontsize=7, ncol=2, loc='upper right', handlelength=1.5)
    letter(ax, 'A')
    ax = fig.add_subplot(gs[0, 1])                                  # B distribution of deltas
    bins = np.linspace(0, 100, 41)
    for code in [m for m in MAIN_MODS if m in deltas]:
        h, e = np.histogram(deltas[code], bins=bins)
        ax.step(e[:-1], np.maximum(h, 0.8), where='post', color=MOD_COLOUR[code], lw=1.3, label=MOD_NAME[code])
    ax.axvline(min_delta, color='#E00000', ls='--', lw=0.9)
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    ax.set_xlabel('max occupancy difference between isoforms (%)')
    ax.set_ylabel('site-sample pairs')
    ax.legend(fontsize=9)
    letter(ax, 'B')
    # C, D: the two genes with the most recurrent isoform-specific sites
    rec = {}
    for key, ds in site_sum.items():
        if np.mean(ds) >= min_delta and len(ds) >= max(2, len(ss) // 3):
            rec[key[0]] = rec.get(key[0], 0) + 1
    top = sorted(rec, key=lambda g: -rec[g])[:2]
    for j, g in enumerate(top):
        ax = fig.add_subplot(gs[1, j])
        pts = [(k_[2], np.mean(d), k_[4]) for k_, d in site_sum.items() if k_[0] == g]
        pts.sort()
        xs_ = np.array([p_[0] for p_ in pts])
        ys_ = np.array([p_[1] for p_ in pts])
        cs_ = [MOD_COLOUR.get(p_[2], '0.5') for p_ in pts]
        ax.scatter(xs_, ys_, c=cs_, s=28, zorder=3)
        big = ys_ >= min_delta
        ax.scatter(xs_[big], ys_[big], facecolor='none', edgecolor='#E00000', s=80, lw=1, zorder=4)
        ax.axhline(min_delta, color='#E76969', ls='--', lw=0.8)
        ax.set_ylim(0, 105)
        ax.set_xlabel('genomic position')
        ax.set_ylabel('isoform Δ (%), mean over samples')
        ax.set_title(gene_name.get(g, g), fontstyle='italic', fontsize=12)
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.tick_params(axis='x', labelsize=9, rotation=20)
        letter(ax, 'CD'[j])
    if not top:
        ax = fig.add_subplot(gs[1, :])
        empty(ax, f'no gene with recurrent isoform-specific sites (Δ ≥ {min_delta}% in ≥ a third of samples)')
    trows = []
    for key, ds in site_sum.items():
        big = sum(d >= min_delta for d in ds)
        if big:
            trows.append((key[0], gene_name.get(key[0], '.'), key[1], key[2], key[3], plain(key[4]), len(ds), big,
                          float(np.mean(ds)), float(np.max(ds))))
    trows.sort(key=lambda r: (-r[7], -r[8]))
    tname = c.table('Fig6_isoform_specific_sites.tsv', ['gene_id', 'gene_name', 'chrom', 'pos', 'strand', 'mod',
                                                         'samples_measured', f'samples_delta_ge{min_delta}', 'mean_delta',
                                                         'max_delta'], trows)
    c.legends['Fig6_isoform_modifications'] = f"""Isoform-specific modification, as the paper's Fig. 5C ({found} samples).
Within each sample, for every site measured on two or more isoforms of a gene (>= 20 reads of each isoform at
the site, gIVT-corrected), Δ is the largest minus the smallest isoform occupancy. This compares isoforms WITHIN
a sample; it does not compare samples.
A, PER SAMPLE: sites with Δ >= {min_delta} points, per modification (2'-O-Me dashed). Deeper samples measure more
sites on more isoforms, so this tracks depth.
B, POOLED: every site x sample with a Δ (m6A, m5C, Ψ, inosine), by Δ; red line at {min_delta}.
C, D, the two genes with the most sites whose Δ, averaged over the samples measuring it, is >= {min_delta} in at
least a third of samples; each point is a site, coloured by modification; circled: mean Δ >= {min_delta}.
Every site with Δ >= {min_delta} in at least one sample: tables/{tname} ({len(trows):,} sites)."""
    notes.append(f'figure 6: isoforms with >= 20 reads at a site; Δ = max - min corrected occupancy over isoforms; '
                 f'top genes: {", ".join(gene_name.get(g, g) for g in top) or "none"}')
    return fig


# ---------------------------------------------------------------- figure 7: poly(A)
def polya_data(c, notes, isoforms=()):
    """Read-level histograms per class and sample, plus per-read tails of the chosen isoforms."""
    out = {k: {} for k, _, _ in COMPARTMENT}
    iso_pt = {i: [] for i in isoforms}
    dropped = [0, 0]                                                 # nuclear tails without a gene, of all
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
        cols = ['compartment', 'pt', 'status', 'gene'] + (['isoform'] if isoforms else [])
        r = read_table(rp, columns=cols)
        used = (r['status'] == 'used').fill_null(False).to_numpy().astype(bool)
        mito = (r['compartment'] == 'mito').fill_null(False).to_numpy().astype(bool)
        glob = r['gene'].is_in(sorted(globin_ids)).fill_null(False).to_numpy().astype(bool) & ~mito
        unassigned = (r['gene'] == '.').fill_null(True).to_numpy().astype(bool)
        pt = numbers(r, ['pt']).ravel()
        used &= ~np.isnan(pt)
        nuclear = ~mito & ~glob & ~unassigned
        dropped[0] += int(np.sum(used & ~mito & unassigned))
        dropped[1] += int(np.sum(used & ~mito))
        for k, m in (('mito', mito), ('globin', glob), ('nuclear', nuclear)):
            v = np.round(pt[used & m]).astype(int)
            if len(v):
                out[k][s] = np.unique(v, return_counts=True)
        if isoforms:
            hit = r['isoform'].is_in(list(isoforms)).fill_null(False).to_numpy().astype(bool) & used
            ids = np.array(r['isoform'].to_list(), dtype=object)[hit] if hit.any() else []
            for i_, p_ in zip(ids, pt[hit]):
                iso_pt[i_].append(p_)
    if dropped[1]:
        notes.append(f'figure 7: {dropped[0]:,} of {dropped[1]:,} non-mitochondrial tails ({100 * dropped[0] / dropped[1]:.1f}%) '
                     'are on reads assigned to no gene (MAPQ < 40, or outside the catalogue) and are left out')
    return out, iso_pt


def isoform_candidates(c, n=100, min_tails=10):
    """Genes whose two best-covered isoforms differ most in median tail, consistently over samples.
    Each: (|difference|, gene_id, gene name, isoform a, isoform b, label a, label b, median a, median b,
    samples a, samples b, tails a, tails b); labels are reference transcript IDs where known."""
    per, names, labels = {}, {}, {}
    for s in c.samples:
        p = c.path('polya', c.ref, s, f'{s}.isoform_polya.bed')
        if not os.path.exists(p):
            continue
        t = read_table(p)
        if t.height == 0:
            continue
        v = numbers(t, ['n_tails', 'median_pt'])
        tx = text(t, 'ref_tx_id') if 'ref_tx_id' in t.columns else [None] * t.height
        cat = text(t, 'structural_category') if 'structural_category' in t.columns else [None] * t.height
        for fid, gid, nm, rt, ct, (nt, md) in zip(text(t, 'feature_id'), text(t, 'gene_id'), text(t, 'ref_gene_name'), tx,
                                                  cat, v):
            if nt >= min_tails and not is_globin(nm):
                per.setdefault(gid, {}).setdefault(fid, []).append((nt, md))
                if nm not in ('.', '', None):
                    names[gid] = nm
                base = rt if rt not in ('.', '', None) else fid
                labels[fid] = f'{base} ({CATEGORY_SHORT.get(ct, ct)})' if ct not in ('.', '', None) else base
    best = []
    for gid, isos in per.items():
        if len(isos) < 2:
            continue
        a, b = sorted(isos, key=lambda f: -sum(x[0] for x in isos[f]))[:2]
        ma, mb = float(np.median([x[1] for x in isos[a]])), float(np.median([x[1] for x in isos[b]]))
        la, lb = (labels[a], labels[b]) if labels[a] != labels[b] else (a, b)   # two models of one reference transcript
        if min(len(isos[a]), len(isos[b])) >= max(2, len(c.samples) // 3):
            best.append((abs(ma - mb), gid, names.get(gid, gid), a, b, la, lb, ma, mb, len(isos[a]),
                         len(isos[b]), int(sum(x[0] for x in isos[a])), int(sum(x[0] for x in isos[b]))))
    best.sort(key=lambda x: -x[0])
    return best[:n]


def fig7_polya(c, notes):
    cands = isoform_candidates(c)
    isoforms = [x for b in cands[:2] for x in (b[3], b[4])]
    hists, iso_pt = polya_data(c, notes, isoforms)
    fig, gs = figure(2, 3)
    grid = np.arange(0, 301)

    def density(vals, cnt):
        d = np.zeros(len(grid))
        inside = vals <= 300
        d[vals[inside]] = cnt[inside]
        return np.convolve(d / max(cnt.sum(), 1), np.ones(5) / 5, mode='same')

    ax = fig.add_subplot(gs[0, 0])                                  # A read-weighted
    top, meds = 0, []
    for k, label, col in COMPARTMENT:
        if not hists[k]:
            continue
        vals, inv = np.unique(np.concatenate([v for v, _ in hists[k].values()]), return_inverse=True)
        cnt = np.bincount(inv, weights=np.concatenate([n for _, n in hists[k].values()]))
        d = density(vals, cnt)
        ax.fill_between(grid, d, color=col, alpha=0.3, lw=0)
        ax.plot(grid, d, color=col, lw=1.2, label=label)
        med = hist_median(vals, cnt)
        ax.axvline(med, color=col, ls='--', lw=1)
        meds.append((med, col))
        top = max(top, d.max())
    ax.set_xlim(0, 300)
    ax.set_ylim(0, top * 1.25 if top else 1)
    for j, (med, col) in enumerate(sorted(meds)):
        ax.text(med + 3, top * (1.18 - 0.1 * j), f'{med:.0f}', color=col, fontsize=11, fontweight='bold', va='top')
    ax.set_xlabel('poly(A) length (nt), every assigned read')
    ax.set_ylabel('density')
    ax.legend(fontsize=8, loc='upper right')
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B gene-level medians
    mp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.median_pt.tsv.gz')
    ntp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.n_tails.tsv.gz')
    if os.path.exists(mp) and os.path.exists(ntp):
        med, ntt = read_table(mp), read_table(ntp)
        ss = [s for s in c.samples if s in med.columns]
        mv, nv = numbers(med, ss), np.nan_to_num(numbers(ntt, ss))
        mv[nv < 10] = np.nan
        names = text(med, 'ref_gene_name')
        chroms = text(med, 'chrom')
        with np.errstate(invalid='ignore'):
            gmed = np.nanmedian(mv, 1)
        cls = np.array(['mito' if ch == 'chrM' else 'globin' if nm in GLOBIN else 'nuclear' for ch, nm in zip(chroms, names)])
        top2, meds2 = 0, []
        for k, label, col in COMPARTMENT:
            v = gmed[(cls == k) & ~np.isnan(gmed)]
            if len(v) < 3:
                continue
            h, e = np.histogram(v, bins=np.arange(0, 305, 5), density=True)
            ax.fill_between(e[:-1] + 2.5, h, color=col, alpha=0.3, lw=0, step=None)
            ax.plot(e[:-1] + 2.5, h, color=col, lw=1.2, label=f'{label} ({len(v):,} genes)')
            m_ = float(np.median(v))
            ax.axvline(m_, color=col, ls='--', lw=1)
            meds2.append((m_, col))
            top2 = max(top2, h.max())
        for j, (m_, col) in enumerate(sorted(meds2)):
            ax.text(m_ + 3, top2 * (1.18 - 0.1 * j), f'{m_:.0f}', color=col, fontsize=11, fontweight='bold', va='top')
        ax.set_ylim(0, top2 * 1.25 if top2 else 1)
        ax.set_xlim(0, 300)
        ax.legend(fontsize=7, loc='upper right')
    ax.set_xlabel('gene median poly(A) (nt), ≥ 10 tails')
    ax.set_ylabel('density')
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C per-sample medians
    for (k, label, col), mk in zip(COMPARTMENT, ('s', '^', 'D')):  # markers too: the colours overlap the batch colours
        xs = [i for i, s in enumerate(c.samples) if s in hists[k]]
        ax.plot(xs, [hist_median(*hists[k][c.samples[i]]) for i in xs], marker=mk, ms=4, lw=0.8, color=col, label=label)
    sample_axis(ax, c.samples, c.batch, c.colours)
    ax.set_ylabel('median poly(A), per read (nt)')
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, loc='best', handlelength=1.5)
    letter(ax, 'C')

    for j, cand in enumerate(cands[:2]):                             # D, E isoform-specific tails
        ax = fig.add_subplot(gs[1, j])
        _, gid, gname, a, b, la, lb = cand[:7]
        data = [np.array(iso_pt.get(a, [])), np.array(iso_pt.get(b, []))]
        if all(len(d) >= 5 for d in data):
            parts = ax.violinplot(data, positions=[0, 1], showextrema=False, widths=0.8)
            for body, col in zip(parts['bodies'], ('0.55', '0.82')):
                body.set_facecolor(col)
                body.set_alpha(0.6)
                body.set_edgecolor('black')
                body.set_linewidth(0.6)
            top = max(np.percentile(d, 99) for d in data)
            for i, d in enumerate(data):
                ax.hlines(np.median(d), i - 0.25, i + 0.25, color='black', lw=1.2)
                ax.text(i, np.percentile(d, 99), f'n = {len(d):,}', ha='center', va='bottom', fontsize=8)
            p = stats.mannwhitneyu(data[0], data[1]).pvalue
            stat_box(ax, [f'Mann-Whitney p = {p:.1e}' if p > 0 else 'Mann-Whitney p < 1e-300'], loc='upper left')
            ax.set_ylim(0, top * 1.45)
        else:
            ax.text(0.5, 0.5, 'fewer than 5 tails on one isoform', transform=ax.transAxes, ha='center', va='center', fontsize=9)
        ax.set_xticks([0, 1], [la, lb], fontsize=8, rotation=12)
        ax.set_ylabel('poly(A) length (nt)')
        ax.set_title(gname, fontsize=11, fontstyle='italic')
        letter(ax, 'DE'[j])
    if not cands:
        ax = fig.add_subplot(gs[1, 0:2])
        empty(ax, 'no gene with two isoforms of >= 10 tails in enough samples')

    ax = fig.add_subplot(gs[1, 2])                                  # F long / short tailed genes
    if os.path.exists(mp) and os.path.exists(ntp):
        want = [g for g in ('DDX17', 'RPS24') if g in names]
        if len(want) < 2:                                           # otherwise the longest / shortest well-covered
            ok = (nv >= 10).all(1) & ~np.isnan(gmed) & (cls == 'nuclear')
            if ok.any():
                idx = np.where(ok)[0]
                want = [names[idx[np.argmax(gmed[idx])]], names[idx[np.argmin(gmed[idx])]]]
        for i, gname in enumerate(want[:2]):
            r_ = names.index(gname)
            y_ = mv[r_]
            xs = np.arange(len(ss))
            col, mk = ('0.15', 'o') if i == 0 else ('0.6', 's')          # greys: the colours mean batches here
            ax.scatter(xs[~np.isnan(y_)], y_[~np.isnan(y_)], color=col, marker=mk, s=26, label=gname, zorder=3)
            ax.axhline(np.nanmedian(y_), color=col, ls='--', lw=0.8)
        sample_axis(ax, ss, c.batch, c.colours)
        ax.set_ylabel('gene median poly(A) (nt)')
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.0, 1.0), handletextpad=0.2)
    letter(ax, 'F')
    fig_batch_legend(fig, c.colours)
    tname = c.table('Fig7_isoform_polya_pairs.tsv', ['abs_difference', 'gene_id', 'gene_name', 'isoform_a', 'isoform_b',
                                                      'label_a', 'label_b', 'median_of_sample_medians_a',
                                                      'median_of_sample_medians_b', 'samples_a', 'samples_b', 'tails_a',
                                                      'tails_b'], cands)
    for s in c.samples:
        if s in hists['nuclear']:
            c.note_sample(s, 'polya_median_per_read_nuclear', hist_median(*hists['nuclear'][s]))
        if s in hists['mito']:
            c.note_sample(s, 'polya_median_per_read_mito', hist_median(*hists['mito'][s]))
    if os.path.exists(mp) and os.path.exists(ntp):
        for j, s in enumerate(ss):
            v = mv[(cls == 'nuclear') & ~np.isnan(mv[:, j]), j]
            if len(v):
                c.note_sample(s, 'polya_median_of_gene_medians_nuclear', float(np.median(v)))
    pair = '; '.join(f'{x[2]} ({x[5]} vs {x[6]})' for x in cands[:2]) or 'none'
    c.legends['Fig7_polya'] = f"""Poly(A) tails, as the paper's Fig. S2.
A, POOLED READS: every tail from every sample on a read assigned to a gene (nuclear without globin; globin) or
on chrM (mitochondrial), smoothed over 5 nt; dashed lines and numbers: medians. Pooling weights deep samples and
highly expressed genes: the nuclear curve is dominated by a few abundant short-tailed transcripts. Reads assigned
to no gene are left out: reads below MAPQ 40, which gene assignment requires (among them the duplicated HBA1 and
HBA2), and reads outside the catalogue.
B, ONE VALUE PER GENE: for each gene with >= 10 tails in a sample, its median tail there; then the median of those
over samples. This is the typical gene, not the typical read, which is why it sits far above A.
C, PER SAMPLE: the per-read median of each class in A.
D, E, POOLED READS of all samples for the two best-covered isoforms of the genes whose isoforms differ most in
median tail, consistently in >= a third of samples: {pair}. Lines: medians. The Mann-Whitney p treats reads as
independent, which they are not across samples: read it as a size, not as a test.
F, PER SAMPLE: gene median tail of DDX17 and RPS24 (the paper's long- and short-tailed examples; otherwise the
longest- and shortest-tailed well-covered genes); dashed: median over samples.
The 100 gene isoform pairs with the largest differences: tables/{tname}. Per-sample medians: Table 1."""
    notes.append(f'figure 7: isoform pairs: {pair}')
    return fig


# ---------------------------------------------------------------- figure 9: isoform usage between samples
CATEGORY_SHORT = {'full-splice_match': 'FSM', 'incomplete-splice_match': 'ISM', 'novel_in_catalog': 'NIC',
                  'novel_not_in_catalog': 'NNC'}


def is_globin(name):
    """True for a catalogue gene whose name (several, comma-joined, for merged genes) includes a globin gene."""
    return bool(set(str(name).split(',')) & set(GLOBIN))


def short_name(name):
    """First name of a merged catalogue gene, with how many more it holds."""
    parts = [p_ for p_ in str(name).split(',') if p_]
    return parts[0] + (f' +{len(parts) - 1}' if len(parts) > 1 else '') if parts else str(name)


def gene_names(c):
    if not hasattr(c, '_gene_names'):
        p = os.path.join(c.cohort, f'cohort.{c.ref}.gene_cpm.tsv.gz')
        c._gene_names = {}
        if os.path.exists(p):
            t = read_table(p, columns=['gene_id', 'gene_name'])
            c._gene_names = {g: n for g, n in zip(text(t, 'gene_id'), text(t, 'gene_name')) if n not in ('.', '', None)}
    return c._gene_names


def fig9_usage(c, group, levels, fdr, notes, min_reads=20):
    base = c.path('isoforms', c.ref, 'high_confidence')
    cp = os.path.join(base, f'cohort.{c.ref}.counts.tsv')
    clp = os.path.join(base, f'cohort.{c.ref}.high_confidence.classification.tsv')
    if not (os.path.exists(cp) and os.path.exists(clp)):
        notes.append(f'figure 9 skipped: needs {cp} and {clp}')
        return None
    cl = read_table(clp)
    col = (lambda k: text(cl, k) if k in cl.columns else ['.'] * cl.height)
    info = {i: (rt, cat) for i, rt, cat in zip(text(cl, 'isoform_id'), col('ref_tx_id'), col('structural_category'))}
    cnt = read_table(cp)
    ss = [s for s in c.samples if s in cnt.columns]
    n = len(ss)
    keep = np.array([t in info for t in text(cnt, 'transcript_id')])
    X = np.nan_to_num(numbers(cnt, ss))[keep]
    tx = [t for t, k in zip(text(cnt, 'transcript_id'), keep) if k]
    gid = [g for g, k in zip(text(cnt, 'gene_id'), keep) if k]
    names = gene_names(c)
    label = {t: (info[t][0] if info[t][0] not in ('.', '', None) else t) for t in tx}
    by_gene = {}
    for i, g in enumerate(gid):
        by_gene.setdefault(g, []).append(i)
    min_samples = max(3, math.ceil(0.75 * n))
    grp = None
    if group and group in c.sheet.columns:
        gv = dict(zip(text(c.sheet, 'sample_id'), text(c.sheet, group)))
        lv = levels or sorted({v for v in gv.values() if v not in ('', '.', None)})
        if len(lv) == 2:
            grp = (lv, np.array([gv.get(s) for s in ss], dtype=object))
    F_all = np.full(X.shape, np.nan)
    rows, sw, seen, examples = [], np.zeros(n), np.zeros(n), []
    globin_skipped = 0
    for g, idx in by_gene.items():
        if len(idx) < 2:
            continue
        if is_globin(names.get(g, '')):                             # near-identical copies: shares measure mapping, not usage
            globin_skipped += 1
            continue
        C = X[idx]
        tot = C.sum(0)
        ok = tot >= min_reads
        if ok.sum() < min_samples:
            continue
        F = np.full(C.shape, np.nan)
        F[:, ok] = C[:, ok] / tot[ok]
        pooled = C[:, ok].sum(1) / tot[ok].sum()
        major = int(np.argmax(pooled))
        exp_sd = float(np.sqrt(np.mean(pooled[major] * (1 - pooled[major]) / tot[ok])))
        if exp_sd == 0:
            continue                                                # one isoform takes every read
        F_all[idx] = F
        fm = F[major, ok]
        obs_sd = float(np.std(fm, ddof=1))
        sub = C[:, ok]
        sub = sub[sub.sum(1) > 0]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p = float(stats.chi2_contingency(sub)[1]) if sub.shape[0] >= 2 else np.nan
        dom = np.argmax(F[:, ok], axis=0)
        sw[ok] += dom != major
        seen[ok] += 1
        row = [g, names.get(g, g), max(1, g.count('ENSG')), len(idx), int(ok.sum()), float(np.median(tot[ok])), tx[idx[major]],
               label[tx[idx[major]]],
               float(pooled[major]), float(fm.min()), float(fm.max()), obs_sd, exp_sd, obs_sd / exp_sd, p,
               int((dom != major).sum())]
        if grp:
            (la, lb), gl = grp
            a_, b_ = fm[gl[ok] == la], fm[gl[ok] == lb]
            if min(len(a_), len(b_)) >= 2:
                row += [float(b_.mean() - a_.mean()), float(stats.ttest_ind(b_, a_, equal_var=False).pvalue)]
            else:
                row += [np.nan, np.nan]
        rows.append(row)
        if ok.all() and np.median(tot) >= 100 and fm.max() - fm.min() >= 0.2:
            examples.append((obs_sd / exp_sd, g, idx, pooled))
    if not rows:
        notes.append(f'figure 9 skipped: no gene with two high-confidence isoforms and >= {min_reads} reads in {min_samples} samples')
        return None
    q = bh([r[14] for r in rows])
    for r, q_ in zip(rows, q):
        r.insert(15, q_)
    header = ['gene_id', 'gene_name', 'catalogue_genes_merged', 'isoforms', 'samples_used', 'median_gene_reads',
              'major_isoform', 'major_label',
              'major_fraction_pooled', 'major_fraction_min', 'major_fraction_max', 'sd_observed', 'sd_expected_binomial',
              'sd_ratio', 'chi2_p', 'chi2_q', 'samples_major_not_dominant']
    if grp:
        (la, lb), _ = grp
        qg = bh([r[-1] for r in rows])
        for r, q_ in zip(rows, qg):
            r.append(q_)
        header += [f'delta_major_fraction_{lb}_vs_{la}', 'welch_p', 'welch_q']
    rows.sort(key=lambda r: -r[13])
    t1 = c.table('Fig9_isoform_usage_genes.tsv', header, rows)
    used = {r[0] for r in rows}
    frows = [[t, g, names.get(g, g), info[t][0], info[t][1]] + [F_all[i, j] for j in range(n)]
             for i, (t, g) in enumerate(zip(tx, gid)) if g in used]
    t2 = c.table('Fig9_isoform_fractions.tsv', ['transcript_id', 'gene_id', 'gene_name', 'ref_tx_id', 'structural_category'] + ss,
                 frows)
    cats = {}
    for t in info:
        cats[info[t][1]] = cats.get(info[t][1], 0) + 1
    t3 = c.table('Fig9_isoform_categories.tsv', ['structural_category', 'high_confidence_isoforms'],
                 sorted(cats.items(), key=lambda kv: -kv[1]))
    for j, s in enumerate(ss):
        if seen[j]:
            c.note_sample(s, 'pct_genes_major_isoform_not_dominant', round(100 * sw[j] / seen[j], 1))

    fig, gs = figure(2, 3, widths=[PANEL[0] * 1.1, PANEL[0] * 1.1, PANEL[0] * 1.2])
    ax = fig.add_subplot(gs[0, 0])                                  # A observed against expected variation
    obs, exp = np.array([r[11] for r in rows]), np.array([r[12] for r in rows])
    ratio, rng_ = obs / exp, np.array([r[10] - r[9] for r in rows])
    hit = (np.nan_to_num(np.array([r[15] for r in rows]), nan=1) < fdr) & (ratio >= 2) & (rng_ >= 0.2)
    lo_, hi_ = max(1e-3, min(exp.min(), obs[obs > 0].min() if (obs > 0).any() else 1e-3)), max(exp.max(), obs.max())
    ax.scatter(exp[~hit], np.maximum(obs[~hit], lo_), s=4, color='0.7', lw=0)
    ax.scatter(exp[hit], np.maximum(obs[hit], lo_), s=6, color='#9E2A2B', lw=0)
    ax.plot([lo_, hi_], [lo_, hi_], ls='--', color='black', lw=0.8)
    ax.plot([lo_, hi_], [2 * lo_, 2 * hi_], ls=':', color='black', lw=0.8)
    ax.set_xscale('log')
    ax.set_yscale('log')
    log_axis(ax.xaxis)
    log_axis(ax.yaxis)
    ax.set_xlabel('expected SD (read sampling)')
    ax.set_ylabel('observed SD between samples')
    stat_box(ax, [f'{len(rows):,} genes', f'{int(hit.sum()):,} ≥ 2x, Δ ≥ 0.2'], loc='upper left')
    ax.set_title('major isoform fraction', fontsize=10)
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B PCA of isoform fractions
    full = ~np.isnan(F_all).any(1) & np.isin(np.array(gid, dtype=object), list(used))
    Z = F_all[full]
    if len(Z) >= 10:
        top = np.argsort(Z.var(1))[::-1][:min(500, len(Z))]
        zc = Z[top] - Z[top].mean(1, keepdims=True)
        _, sv, vt = np.linalg.svd(zc, full_matrices=False)
        pcs, fr = vt.T * sv, sv ** 2 / (sv ** 2).sum()
        for i, s in enumerate(ss):
            ax.scatter(pcs[i, 0], pcs[i, 1], color=c.colours[c.batch[s]], s=40, zorder=3)
            ax.annotate(s[-2:], (pcs[i, 0], pcs[i, 1]), fontsize=8, xytext=(3, 2), textcoords='offset points')
        ax.set_xlabel(f'PC1 ({100 * fr[0]:.1f}%)')
        ax.set_ylabel(f'PC2 ({100 * fr[1]:.1f}%)')
        ax.set_title(f'isoform fractions, {len(top)} most variable isoforms', fontsize=10)
    else:
        empty(ax, 'too few genes measured in every sample')
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C per sample: major isoform not dominant
    pct_sw = np.where(seen > 0, 100 * sw / np.maximum(seen, 1), np.nan)
    ax.bar(range(n), pct_sw, color=[c.colours[c.batch[s]] for s in ss], width=0.8)
    ax.axhline(np.nanmedian(pct_sw), ls='--', color='black', lw=0.8)
    sample_axis(ax, ss, c.batch, c.colours)
    ax.set_ylabel('genes, other top isoform (%)')
    letter(ax, 'C')

    examples.sort(key=lambda e: -e[0])
    shades = ['#3B0000', '#B4766F', '#E8CFC4']
    for j, (rt, g, idx, pooled) in enumerate(examples[:3]):         # D-F examples
        ax = fig.add_subplot(gs[1, j])
        order = list(np.argsort(pooled)[::-1])
        C = X[idx]
        F = C / C.sum(0)
        bottom = np.zeros(n)
        for k, o in enumerate(order[:3]):
            t = tx[idx[o]]
            ax.bar(range(n), F[o], bottom=bottom, color=shades[k], width=0.8,
                   label=f'{label[t]} ({CATEGORY_SHORT.get(info[t][1], info[t][1])})')
            bottom += F[o]
        if len(order) > 3:
            ax.bar(range(n), 1 - bottom, bottom=bottom, color='white', edgecolor='0.6', lw=0.5, width=0.8,
                   label=f'{len(order) - 3} other')
        sample_axis(ax, ss, c.batch, c.colours)
        ax.set_ylim(0, 1)
        ax.set_ylabel('fraction of the gene’s reads')
        ax.set_title(f'{short_name(names.get(g, g))}  (median {np.median(C.sum(0)):,.0f} reads, SD {rt:.1f}x expected)',
                     fontsize=9, fontstyle='italic')
        ax.legend(fontsize=6.5, loc='upper left', bbox_to_anchor=(0, -0.32), ncol=2, handlelength=1, columnspacing=0.8)
        letter(ax, 'DEF'[j])
    if not examples:
        ax = fig.add_subplot(gs[1, :])
        empty(ax, 'no gene with >= 100 reads in every sample and a major-isoform range >= 0.2')
    fig_batch_legend(fig, c.colours)
    shown = '; '.join(short_name(names.get(e[1], e[1])) for e in examples[:3]) or 'none'
    c.legends['Fig9_isoform_usage'] = f"""Isoform usage compared BETWEEN samples ({n} samples, {len(rows):,} genes; {globin_skipped} globin genes left out,
since reads split between near-identical globin copies by mapping, not by isoform use).
For every gene with >= 2 high-confidence isoforms, the fraction of the gene's reads on each isoform, per sample
(IsoQuant requantification counts; the gene total is the sum over its high-confidence isoforms). A sample is
used for a gene where the gene has >= {min_reads} reads; genes need >= {min_samples} such samples. The major
isoform has the largest share of the gene's reads pooled over samples.
A, ONE POINT PER GENE: standard deviation of the major isoform's fraction across samples (y) against the SD that
read sampling alone would give at each sample's depth (binomial, x). Dashed: equal; dotted: twice. Points
far above the line vary between samples more than counting noise explains; dark red: >= 2x, a range of the
major fraction >= 0.2, and chi-square FDR < {fdr:g} (the chi-square test assumes binomial counts and so finds
almost everything: rank by the SD ratio).
B, PCA of samples on isoform fractions (the 500 most variable isoforms among genes measured in every sample),
the isoform counterpart of figure 2D. C, PER SAMPLE: % of measured genes whose major isoform is not that
sample's most used isoform. D-F, the three genes with >= 100 reads in every sample, a major-fraction range
>= 0.2 and the highest SD ratio ({shown}): stacked fractions of their three main isoforms (reference
transcript ID where known; FSM / ISM / NIC / NNC: SQANTI categories).
Exploratory. For tests between groups use DRIMSeq or satuRn on the counts in
isoforms/{c.ref}/high_confidence/cohort.{c.ref}.counts.tsv. Tables: {t1} (per gene{', with a Welch test of the major fraction between ' + ' and '.join(grp[0]) if grp else ''}), {t2} (per isoform and sample), {t3}."""
    notes.append(f'figure 9: {len(rows):,} genes ({globin_skipped} globin genes left out); {int(hit.sum()):,} with SD >= 2x '
                 f'expected and range >= 0.2; examples {shown}')
    return fig


def hist_median(vals, cnt):
    cum = np.cumsum(cnt)
    return float(vals[np.searchsorted(cum, cum[-1] / 2)])
def expression(c):
    cpm = read_table(os.path.join(c.cohort, f'cohort.{c.ref}.gene_cpm.tsv.gz'))
    es = [s for s in c.samples if s in cpm.columns]
    m = np.nan_to_num(numbers(cpm, es))
    keep = np.median(m, axis=1) > 1
    names = [n if n not in ('.', '', None) else g for g, n in zip(text(cpm, 'gene_id'), text(cpm, 'gene_name'))]
    return es, np.log2(m[keep] + 1), [n for n, k in zip(names, keep) if k]



def compare_groups(c, compare):
    names = sorted(set(c.batch.values()))
    if len(names) < 2:
        return None
    a = compare or names[-1]
    if a not in names:
        sys.exit(f'cohort_plots.py: --compare-batch {a} is not a batch ({", ".join(names)})')
    return a


def figS1_batch(c, compare, notes):
    _, ms, pct, mods, _ = c.mods_matrix()
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
        ax = fig.add_subplot(gs[2, :2])                                 # I the r of each modification
        codes = [m for m in MAIN_MODS + OME_MODS if not np.isnan(rvals.get(m, (np.nan, 0))[0])]
        ax.bar(range(len(codes)), [rvals[m][0] for m in codes], color=[MOD_COLOUR[m] for m in codes], width=0.7)
        for i, m in enumerate(codes):
            ax.text(i, rvals[m][0] + 0.01, f'{rvals[m][0]:.2f}', ha='center', va='bottom', fontsize=9)
            ax.text(i, 0.02, f'{rvals[m][1]:,}', ha='center', va='bottom', fontsize=7, color='white', rotation=90)
        ax.set_xticks(range(len(codes)), [MOD_NAME[m] for m in codes])
        ax.set_ylim(0, 1.08)
        ax.set_ylabel(f'r, {a} vs others')
        ax.set_title('per modification (sites in bars)', fontsize=11)
        letter(ax, 'I')
        notes.append(f'figure S1: {a} ({len(ia)} samples) vs the other batches ({len(ib)} samples); sites measured '
                     f'in >= 2 samples of each side')
        c.legends['FigS1_batch_comparison'] = f"""One batch against the rest: {a} ({len(ia)} samples) against the other batches ({len(ib)}).
A-H, per cohort site measured (>= 20 reads, a corrected value) in >= 2 samples on each side: the mean corrected
occupancy over the other batches' samples (x) and over {a}'s (y); log colour, one scale for all panels; dashed
y = x; n sites, Pearson r, RMSE. I, r per modification, sites in the bars. Figure 5 makes the same comparison
for every batch at once (within against between batches)."""
    else:
        notes.append('figure S1: one batch only, so panels A-H are left empty')

    return fig


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


def fig8_groups(c, group, levels, fdr, out, notes):
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

    rows, ms, pct, mods, _ = c.mods_matrix()
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
               zip(rows['chrom'], rows['start'], mods, full) if f]
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
    notes.append(f'figure 8: {label}, {len(ga)} vs {len(gb)} samples; Welch + BH, exploratory')
    c.legends['Fig8_group_comparison'] = f"""{group}: {b} ({len(gb)} samples) against {a} ({len(ga)} samples). Per feature, Welch's t-test on the samples'
values, Benjamini-Hochberg FDR; coloured: FDR < {fdr:g}. Exploratory: use DESeq2 or edgeR on the counts for
expression, and a beta-binomial model on n_mod / n_valid for modifications, before claiming differences.
A, expression, log2(CPM + 1), genes with median CPM > 1. B, z-scores of the 30 genes with the lowest p.
C, modification occupancy at cohort sites measured in >= 3 samples per group. D, gene median poly(A), genes with
>= 10 tails in >= 3 samples per group. Every test: tables/expression.{tag}.tsv, tables/modifications.{tag}.tsv,
tables/polya.{tag}.tsv."""
    return fig


# ---------------------------------------------------------------- main
FIGURES = ['Fig1', 'Fig2', 'Fig3', 'Fig4', 'Fig5', 'Fig6', 'Fig7', 'Fig8', 'Fig9', 'FigS1', 'FigS2']


def table1(c):
    """Per-sample summary (the paper's Table 1 role): depth, composition, sites per modification, tails."""
    sites = {}
    if os.path.exists(os.path.join(c.cohort, f'cohort.{c.ref}.mods.pct_fpcorr.tsv.gz')):
        _, ms, pct, mods, _ = c.mods_matrix()
        for code in MAIN_MODS + OME_MODS:
            called = (np.nan_to_num(pct[mods == code]) >= 20).sum(0)
            for s, n in zip(ms, called):
                sites.setdefault(s, {})[f'sites_{plain(code)}'] = int(n)
    cols = (['mapped_reads', 'globin_pct', 'mito_pct', 'non_globin_non_mito_reads', 'genes_ge10_reads']
            + [f'sites_{plain(m)}' for m in MAIN_MODS + OME_MODS] + ['sites_all']
            + ['polya_median_per_read_nuclear', 'polya_median_of_gene_medians_nuclear', 'polya_median_per_read_mito',
               'pct_genes_major_isoform_not_dominant'])
    rows = []
    for s in c.samples:
        d = dict(c.summary.get(s, {}))
        d.update(sites.get(s, {}))
        if s in sites:
            d['sites_all'] = sum(sites[s].values())
        rows.append([s, c.batch[s]] + [d.get(k, float('nan')) for k in cols])
    return c.table('Table1_cohort_summary.tsv', ['sample', 'batch'] + cols, rows)


def write_legends(c, made):
    """One legend file per figure made (legends/), and all of them, in figure order, in figure_legends.md."""
    d = os.path.join(c.out, 'legends')
    os.makedirs(d, exist_ok=True)
    for name in made:
        if name in c.legends:
            with open(os.path.join(d, f'{name}.md'), 'w') as fh:
                fh.write(f'## {name}\n\n{c.legends[name]}\n')
    order = {k: i for i, k in enumerate(FIGURES)}
    files = sorted((f for f in os.listdir(d) if f.endswith('.md')), key=lambda f: (order.get(f.split('_')[0], 99), f))
    with open(os.path.join(c.out, 'figure_legends.md'), 'w') as out:
        out.write('# Figure legends\n\nWritten by cohort_plots.py from the data each figure was drawn from. POOLED means '
                  'every sample\'s values are put together before counting; PER SAMPLE means one value per sample.\n\n')
        for f in files:
            out.write(open(os.path.join(d, f)).read() + '\n')



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True)
    ap.add_argument('--samplesheet', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ref', default='GRCh38')
    ap.add_argument('--batch-column', default='batch')
    ap.add_argument('--read-stats', help='directory of read_stats.sh outputs (default: <results>/read_stats)')
    ap.add_argument('--mask-contigs', default=SITE_FILTER_DEFAULT,
                    help=f'regex of contigs to drop, POSIX ERE (awk reads it too; default: {SITE_FILTER_DEFAULT}); "" for none')
    ap.add_argument('--mask-bed', help='BED of regions to drop from the modification figures')
    ap.add_argument('--group', help='samplesheet column to compare (two levels): adds figure 8 and a test in figure 9\'s table')
    ap.add_argument('--levels', nargs=2, metavar=('A', 'B'))
    ap.add_argument('--fdr', type=float, default=0.1)
    ap.add_argument('--compare-batch', help='adds figure S1: this batch against all others')
    ap.add_argument('--only', nargs='+', choices=FIGURES, help='make only these figures')
    ap.add_argument('--threads', type=int, default=min(8, os.cpu_count() or 1),
                    help='files binned at once for figure 3 (default: up to 8)')
    ap.add_argument('--style', default=LAB_STYLE, help=f'matplotlib style sheet (default: {LAB_STYLE}); "" for none')
    args = ap.parse_args()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s %(message)s', datefmt='%H:%M:%S'))
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    log.propagate = False
    os.makedirs(args.out, exist_ok=True)
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='(Mean of empty slice|All-NaN slice|Degrees of freedom|invalid value)')
    notes = apply_style(args.style)
    mask = SiteMask(args.mask_contigs, args.mask_bed)
    c = Cohort(args.results, args.ref, args.samplesheet, args.batch_column, mask, args.out)
    rs_dir = args.read_stats or os.path.join(args.results, 'read_stats')
    figures = [('Fig1', 'Fig1_read_qc', lambda: fig1_reads(c, rs_dir, notes)),
               ('Fig2', 'Fig2_cohort_overview', lambda: fig2_overview(c, notes)),
               ('Fig3', 'Fig3_ivt_correction', lambda: fig3_ivt(c, args.out, args.threads, notes)),
               ('Fig4', 'Fig4_modification_calls', lambda: fig4_mods(c, notes)),
               ('Fig5', 'Fig5_replicate_concordance', lambda: fig5_concordance(c, notes)),
               ('Fig6', 'Fig6_isoform_modifications', lambda: fig6_isoforms(c, notes)),
               ('Fig7', 'Fig7_polya', lambda: fig7_polya(c, notes))]
    if args.group:
        figures.append(('Fig8', 'Fig8_group_comparison', lambda: fig8_groups(c, args.group, args.levels, args.fdr, c.tables, notes)))
    figures.append(('Fig9', 'Fig9_isoform_usage', lambda: fig9_usage(c, args.group, args.levels, args.fdr, notes)))
    if args.compare_batch:
        figures.append(('FigS1', 'FigS1_batch_comparison', lambda: figS1_batch(c, args.compare_batch, notes)))
    figures.append(('FigS2', 'FigS2_jaccard', lambda: figS2_jaccard(c, notes)))
    made = []
    for key, name, make in figures:
        if args.only and key not in args.only:
            continue
        log.info('%s ...', name)
        fig = make()
        if fig is None:
            continue
        fig.savefig(os.path.join(args.out, f'{name}.pdf'), format='pdf', bbox_inches='tight', transparent=True)
        plt.close(fig)
        made.append(name)
    write_legends(c, made)
    if args.only:
        notes.append('Table 1 not rewritten: it needs every figure (run without --only)')
    else:
        notes.append(f'tables in {c.tables}: {table1(c)} and one or more per figure; legends in figure_legends.md')
    for rule, label in (('contig', 'contig'), ('mask-bed', '--mask-bed')):
        parts = [f'{plain(m)} {mask.removed[(rule, m)]:,}' for m in MAIN_MODS + OME_MODS if (rule, m) in mask.removed]
        if parts:
            notes.append(f'masked by {label}: cohort sites ' + ', '.join(parts))
    for n in notes:
        print(f'note: {n}', file=sys.stderr)
    print(f'cohort_plots.py: {len(made)} figures in {args.out}: {" ".join(f"{m}.pdf" for m in made)}')


if __name__ == '__main__':
    main()