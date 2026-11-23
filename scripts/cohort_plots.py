#!/usr/bin/env python3
"""Cohort figures for RNA_analysis.nf results (v5), laid out like Ghohabi Esfahani, Stein, Akeson et al.
(Genome Biology 2026, 27:190): bold panel letters, the lab style sheet and Helvetica, the paper's colours,
log-scaled 2D histograms, dashed medians, statistics as plain text (no boxes). One vector PDF per figure,
editable in Illustrator, a legend per figure and the numbers behind every panel.

    cohort_plots.py --results <results> --samplesheet <analysis_samplesheet.csv> --out <dir>
                    [--group <column> [--levels A B] [--fdr 0.1]] [--covariates COL ...]
                    [--reference GRCh38.fa] [--gtf gencode.gtf] [--errortable gIVT_9mer.tsv.gz]
                    [--composition cbc.tsv | --composition-columns COL ...] [--known-sites MOD=BED ...]
                    [--repeats rmsk.bed] [--sample-qc <dir>] [--read-stats <dir>]
                    [--mask-bed regions.bed] [--mask-contigs REGEX] [--compare-batch <batch>]
                    [--only Fig3 Fig12 ...] [--threads 8] [--ref GRCh38] [--style <mplstyle>] [--font-dir DIR]

Figures kept from v4, with the fixes from the cohort review:
  Fig1_read_qc                identity (dashed medians) and indels per batch; PER-SAMPLE medians beside them;
                              identity by aligned length; substitution matrix per batch (paper Fig. 1)
  Fig2_cohort_overview        read composition, depth, sites against depth, expression PCA and a clustered
                              correlation matrix, on CPM renormalised WITHOUT globin and chrM genes
  Fig3_ivt_correction         Dorado against gIVT-corrected occupancy (paper Fig. 3, Table 4); per-sample
                              sites before and after as paired points, not bars on a log axis
  Fig4_modification_calls     sites per occupancy bin (per-sample points), how consistently measuring samples
                              call a site, sites per gene, gene co-occurrence against an expression-matched null
  Fig5_between_donor_agreement  pairwise ICC(2,1); sites chosen from the OTHER samples (leave-pair-out), so the
                              selection does not condition on the pair's own values
  Fig6_isoform_modifications  isoform differences tested (CMH across samples) and set against the binomial null
  Fig7_polya                  pooled and per-sample tails (ECDFs), isoform tails as per-sample paired medians
  Fig8_group_comparison       design balance, composition by group, moderated-t expression and poly(A)
                              models with covariates, random-effects logit model for occupancy, forest plots,
                              site-level change split into within-isoform and isoform-usage parts, power
  Fig9_isoform_usage          v4 panels plus how much of each gene's usage variation tracks composition
  FigS1_batch_comparison      (--compare-batch) one batch against the others
  FigS2_jaccard               pairwise Jaccard, leave-pair-out site selection
New:
  Fig10_cohort_drivers        leukocyte composition, what drives each data layer (R2 of PCs and per-sample
                              summaries on every covariate), sex check, genotype concordance, design balance
  Fig11_integrity_signal      5'-completeness by transcript length, read lengths, read categories,
                              translocation speed and qscore, modkit thresholds and fail fraction per sample
  Fig12_ivt_calibration       funnel plots against the binomial gIVT null, randomized-PIT null p-value
                              histograms, lambda and pi0 per sample, the 20% rule against FDR across coverage
  Fig13_controls_plausibility mt-rRNA known sites, model specificity at them, cross-talk offset profiles,
                              DRACH fraction, metagene, known-site overlap and occupancy, donor SNPs in the 9-mer
  Fig14_measurement_vs_biology split-half against donor pairs at matched depth, observed against binomial SD,
                              thinning rarefaction, power
  Fig15_bayesian_occupancy    measurement model: n_mod ~ Bin(n, s*pi + fp*(1-pi)), empirical-Bayes prior per
                              modification, posterior calls against the 20% rule, shrinkage, forest plot
  Fig16_orthogonal_sites      (--known-sites) calls in each orthogonal set by occupancy, and occupancy against it

Inputs are RNA_analysis.nf's outputs under --results (cohort/, mods/, polya/, isoforms/, variants/, qc/).
Optional inputs switch panels on; a panel without its input says what it needs:
  --reference (uncompressed FASTA + .fai)   DRACH fraction; exact gIVT rates (with --errortable)
  --errortable + --reference                exact 9-mer FP rate for sites gIVT clipped to 0 (Fig 12 null)
  --gtf                                     metagene (Fig 13); 5' completeness needs cohort_extract.py
  --composition / --composition-columns     measured composition (CBC or deconvolution) instead of marker scores
  --known-sites MOD=BED                     orthogonal sites (Fig 16): BED6 from orthogonal_to_bed.py, occupancy in
                                            column 5 (0-1 or 0-100, '.' for none); several sets per modification
  --repeats                                 RepeatMasker BED, for the inosine-in-Alu fraction
  --sample-qc                               cohort_extract.py outputs (default <results>/sample_qc)
Variants come from variants/<ref>/<s>/<s>.variants.tsv.gz (known/novel PASS SNVs, editing sites excluded).

Fig 12-15 stream every <sample>.cov20_fpcorr.bed once (pandas, in chunks) and cache the summaries in
<out>/.cache keyed by file size, time and settings; Fig 3 keeps its awk pass and cache. Everything is
exploratory unless a legend says otherwise; the Welch and moderated tests treat samples, not reads, as units.
"""
import argparse
import csv
import glob
import gzip
import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib import font_manager                              # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm   # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter, MaxNLocator   # noqa: E402
from scipy import stats                                          # noqa: E402
from scipy.special import digamma, expit, logit, polygamma       # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage         # noqa: E402
from scipy.spatial.distance import squareform                    # noqa: E402

log = logging.getLogger('cohort_plots')
LAB_STYLE = '/projects/Genometechlab/matplotlib_stylesheets/styles/genometechlab_main.mplstyle'
FONT_DIRS = ['/projects/Genometechlab/matplotlib_stylesheets/fonts']
FONT_FALLBACKS = ['Helvetica', 'Arial', 'Liberation Sans', 'Nimbus Sans', 'TeX Gyre Heros']
CLUSTER_REFERENCE = '/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
CLUSTER_GTF = '/projects/Genometechlab/referenceFastaFiles/human/gencode/gencode.v47.annotation.gtf'
CLUSTER_ERRORTABLE = ('/projects/Genometechlab/referenceFastaFiles/gIVT/'
                      'GM12878_genomic_IVT_refmatch_9mer_all_threshold.dorado_1.0.tsv.gz')
GLOBIN = ['HBA1', 'HBA2', 'HBB', 'HBD', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ', 'HBE1']
MITO_CONTIGS = ('chrM', 'MT', 'chrMT', 'M')

# the paper's colours (sampled from its figures); 2'-O-Me: darker shades of the base's colour
MODS = [('a', 'm$^6$A', '#398DC3'), ('m', 'm$^5$C', '#D593BD'), ('17802', 'Ψ', '#DB8151'),
        ('17596', 'Ino', '#4BB193'), ('69426', 'Am', '#1F5F8B'), ('19228', 'Cm', '#9C5B8E'),
        ('19227', 'Um', '#9E4F25'), ('19229', 'Gm', '#8A7A3B')]
MOD_NAME = {c: n for c, n, _ in MODS}
MOD_COLOUR = {c: k for c, _, k in MODS}
MAIN_MODS = ['a', 'm', '17802', '17596']
OME_MODS = ['69426', '19228', '19227', '19229']
ALL_MODS = MAIN_MODS + OME_MODS
MOD_INDEX = {c: i for i, c in enumerate(ALL_MODS)}
MOD_BASE = {'a': 'A', '17596': 'A', '69426': 'A', 'm': 'C', '19228': 'C', '17802': 'T', '19227': 'T', '19229': 'G'}
BATCH_COLOURS = ['#5769AB', '#EDB157', '#935193', '#6FBD81', '#5D9FC9', '#E76969', '#AB7B4B', '#63CFE1']
GROUP_MARKERS = ['o', 's', '^', 'D', 'v', 'P']
COMPARTMENT = [('nuclear', 'nuclear (no globin)', '#5769AB'), ('globin', 'globin', '#BD7B7B'),
               ('mito', 'mitochondrial', '#EDB157')]
MATRIX_CMAP = LinearSegmentedColormap.from_list('paper_matrix', ['#3B0000', '#7B3333', '#C9A9A2', '#FBF3EE'])
HEAT_CMAP = LinearSegmentedColormap.from_list('paper_heat', ['#FBF3EE', '#C9A9A2', '#7B3333', '#3B0000'])
PANEL = (3.25, 2.25)                                               # half the style sheet's 6.5 x 4.5
OCC_BINS = list(range(20, 100, 10))                                # 20-30 ... 90-100
SITE_FILTER_DEFAULT = r'(_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-)'

# Known mitochondrial rRNA modifications (rCRS = GRCh38 chrM, 1-based). Ψ1397 = m.3067, m5C841 = m.1488,
# m4C839 = m.1486, m5U429 = m.1076, m6,6A936/937 = m.1592/1593 as given by Van Haute et al. (NAR 2019,
# METTL15); 16S starts at m.1671, so Gm1145 = m.2815, Um1369 = m.3039, Gm1370 = m.3040, m1A947 = m.2617
# (derived with the same offset). The code is the Dorado model expected to call the site, None = no model.
MT_RRNA_SITES = [(1076, None, 'm$^5$U429', '12S', 'T'), (1486, None, 'm$^4$C839', '12S', 'C'),
                 (1488, 'm', 'm$^5$C841', '12S', 'C'), (1592, None, 'm$^6_2$A936', '12S', 'A'),
                 (1593, None, 'm$^6_2$A937', '12S', 'A'), (2617, None, 'm$^1$A947', '16S', 'A'),
                 (2815, '19229', 'Gm1145', '16S', 'G'), (3039, '19227', 'Um1369', '16S', 'T'),
                 (3040, '19229', 'Gm1370', '16S', 'G'), (3067, '17802', 'Ψ1397', '16S', 'T')]

# Marker genes for relative composition scores when no CBC or deconvolution is given. Scores are the mean
# z-score of log2(CPM + 1) over a cell type's markers: relative between samples, not proportions.
CELL_MARKERS = {
    'neutrophil': ['FCGR3B', 'CSF3R', 'CXCR2', 'CXCR1', 'FPR1', 'AQP9', 'CMTM2', 'PROK2'],
    'eosinophil': ['CLC', 'IL5RA', 'CCR3', 'SIGLEC8', 'HRH4'],
    'monocyte': ['CD14', 'CSF1R', 'FCN1', 'VCAN', 'CD68'],
    'T cell': ['CD3E', 'CD3D', 'CD3G', 'CD2', 'LCK', 'CD28'],
    'B cell': ['MS4A1', 'CD79A', 'CD79B', 'CD19', 'BANK1', 'PAX5'],
    'NK cell': ['GNLY', 'NKG7', 'KLRD1', 'PRF1', 'KLRF1'],
    'platelet': ['PPBP', 'PF4', 'ITGA2B', 'GP9', 'GP1BA', 'TUBB1'],
    'reticulocyte': ['ALAS2', 'SLC4A1', 'CA1', 'AHSP', 'EPB42'],
}
SEX_Y = ['RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY', 'ZFY', 'USP9Y']
CATEGORY_SHORT = {'full-splice_match': 'FSM', 'incomplete-splice_match': 'ISM', 'novel_in_catalog': 'NIC',
                  'novel_not_in_catalog': 'NNC', 'full_splice_match': 'FSM', 'incomplete_splice_match': 'ISM'}
FIGURES = ['Fig1', 'Fig2', 'Fig3', 'Fig4', 'Fig5', 'Fig6', 'Fig7', 'Fig8', 'Fig9', 'Fig10', 'Fig11', 'Fig12',
           'Fig13', 'Fig14', 'Fig15', 'Fig16', 'FigS1', 'FigS2']


# ---------------------------------------------------------------- tables
def _header(path, sep):
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        line = fh.readline().rstrip('\n')
    return [c.lstrip('#') if i == 0 else c for i, c in enumerate(line.split(sep))]


def read_table(path, sep='\t', columns=None, numeric=()):
    """A DataFrame of text (empty cell = ''), except the columns in `numeric`, read as float. A '#' before
    the first header name (the pipeline's BED files) is dropped."""
    names = _header(path, sep)
    use = [c for c in (columns or names) if c in names]
    num = set(numeric) if not callable(numeric) else {c for c in use if numeric(c)}
    dtype = {c: (np.float64 if c in num else str) for c in use}
    kw = dict(sep=sep, header=0, names=names, usecols=use, dtype=dtype, keep_default_na=False,
              na_values={c: ['', 'NA', 'nan', 'NaN', '.', 'None'] for c in num}, quoting=csv.QUOTE_NONE)
    try:
        df = pd.read_csv(path, **kw)
    except ValueError:                                              # a numeric column with text in it
        kw['dtype'] = {c: str for c in use}
        kw['na_values'] = {}
        df = pd.read_csv(path, **kw)
        for c in num:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def read_matrix(path, keys):
    """A cohort matrix: the key columns as text, every other column as float."""
    return read_table(path, numeric=lambda c: c not in keys)


def numbers(df, cols):
    out = np.empty((len(df), len(cols)))
    for j, c in enumerate(cols):
        v = df[c]
        out[:, j] = v.to_numpy(dtype=float) if v.dtype.kind == 'f' else pd.to_numeric(v, errors='coerce').to_numpy(dtype=float)
    return out


def text(df, col):
    return [None if (v is None or v == '' or (isinstance(v, float) and v != v)) else v for v in df[col].tolist()]


def strarr(df, col):
    return np.asarray(df[col].fillna('').astype(str).to_numpy(), dtype=object)


def write_rows(path, header, rows):
    with open(path, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for r in rows:
            out.write('\t'.join('NA' if (isinstance(v, (float, np.floating)) and math.isnan(v)) else
                                (f'{v:.6g}' if isinstance(v, (float, np.floating)) else str(v)) for v in r) + '\n')


def opener(path):
    return gzip.open(path, 'rt') if path.endswith(('.gz', '.bgz')) else open(path)


# ---------------------------------------------------------------- statistics
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
    """Welch's t-test per row, B minus A, ignoring missing values; rows with fewer than min_n values on
    either side get NaN. Returns (difference of means, p)."""
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


def wilson(k, n, z=1.96):
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid='ignore', divide='ignore'):
        p = k / n
        d = 1 + z * z / n
        c = (p + z * z / (2 * n)) / d
        h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def trigamma_inverse(x):
    """limma's trigammaInverse: y with trigamma(y) = x (Newton)."""
    x = np.asarray(x, float)
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = polygamma(1, y)
        dif = tri * (1 - tri / x) / polygamma(2, y)
        y = y + dif
        if np.all(np.abs(dif / y) < 1e-8):
            break
    return y


def squeeze_var(s2, df):
    """limma's squeezeVar (Smyth 2004): empirical-Bayes prior (d0, s0^2) from the residual variances and
    each feature's posterior variance."""
    ok = np.isfinite(s2) & (s2 > 0) & (df > 0)
    if ok.sum() < 3:
        return s2, np.inf, np.nan
    z = np.log(s2[ok])
    e = z - digamma(df[ok] / 2) + np.log(df[ok] / 2)
    emean = e.mean()
    evar = ((e - emean) ** 2).sum() / (len(e) - 1) - polygamma(1, df[ok] / 2).mean()
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        s0 = math.exp(emean + digamma(d0 / 2) - math.log(d0 / 2))
    else:
        d0, s0 = np.inf, math.exp(emean)
    if np.isinf(d0):
        post = np.where(np.isfinite(s2), s0, np.nan)
    else:
        with np.errstate(invalid='ignore'):
            post = (d0 * s0 + df * s2) / (d0 + df)
    return post, d0, s0


def moderated_lm(Y, X, coef, min_extra_df=1):
    """limma-style per-row linear model with missing values (NaN rows left out per feature) and
    empirical-Bayes moderated t for one coefficient. Y: features x samples; X: samples x p.
    Returns (estimate, moderated t, p, residual df, prior df)."""
    Y = np.asarray(Y, float)
    W = (~np.isnan(Y)).astype(float)
    Yz = np.nan_to_num(Y)
    G, p = Y.shape[0], X.shape[1]
    A = np.einsum('gn,np,nq->gpq', W, X, X) + 1e-10 * np.eye(p)
    b = np.einsum('gn,np,gn->gp', W, X, Yz)
    k = W.sum(1)
    dfres = k - p
    est = np.full(G, np.nan)
    unscaled = np.full(G, np.nan)
    s2 = np.full(G, np.nan)
    good = dfres >= min_extra_df
    if good.any():
        Ai = np.linalg.inv(A[good])
        beta = np.einsum('gpq,gq->gp', Ai, b[good])
        fit = np.einsum('np,gp->gn', X, beta)
        rss = (W[good] * (Yz[good] - fit) ** 2).sum(1)
        est[good] = beta[:, coef]
        unscaled[good] = Ai[:, coef, coef]
        s2[good] = rss / dfres[good]
    post, d0, _ = squeeze_var(s2, dfres.astype(float))
    with np.errstate(invalid='ignore', divide='ignore'):
        t = est / np.sqrt(post * unscaled)
        dft = dfres + (d0 if np.isfinite(d0) else 1e6)
        pv = 2 * stats.t.sf(np.abs(t), dft)
    pv[~good] = np.nan
    return est, t, pv, dfres, d0


TAU2_GRID = np.concatenate([[0.0], np.logspace(-3, 1.3, 24)])


def re_logit(m, n, X, coef=None, chunk=4000):
    """Random-effects logit model per site, the sampling variance treated as known measurement error:
        y_ij = logit((m_ij + 0.5) / (n_ij + 1)) = x_j' beta_i + u_ij + e_ij,  u ~ N(0, tau_i^2),
        e ~ N(0, v_ij), v_ij = 1/(m_ij + 0.5) + 1/(n_ij - m_ij + 0.5).
    tau^2 by REML on a grid, with an empirical-Bayes prior over the grid shared by all sites (EM), which
    steadies tau^2 at small n. m, n: sites x samples (NaN = not measured); m may be FP-corrected counts.
    Returns dict(tau2, beta (coef or None), se, z (as 't'), p, k)."""
    S, N = m.shape
    p = X.shape[1]
    G = len(TAU2_GRID)
    M = np.isfinite(m) & np.isfinite(n) & (np.nan_to_num(n) > 0)
    mm = np.clip(np.nan_to_num(m), 0, None)
    nn = np.nan_to_num(n)
    mm = np.minimum(mm, nn)
    Y = np.where(M, np.log((mm + 0.5) / (nn - mm + 0.5)), 0.0)
    V = np.where(M, 1 / (mm + 0.5) + 1 / (nn - mm + 0.5), 1.0)
    k = M.sum(1)
    ll = np.full((S, G), -np.inf)
    bc = np.full((S, G), np.nan)
    vc = np.full((S, G), np.nan)
    ok = k >= p + 2
    for s0 in range(0, S, chunk):
        sl = slice(s0, min(S, s0 + chunk))
        idx = np.where(ok[sl])[0] + s0
        if not len(idx):
            continue
        Ys, Vs, Ms = Y[idx], V[idx], M[idx].astype(float)
        for g, t2 in enumerate(TAU2_GRID):
            W = Ms / (t2 + Vs)
            A = np.einsum('sn,np,nq->spq', W, X, X) + 1e-10 * np.eye(p)
            b = np.einsum('sn,np,sn->sp', W, X, Ys)
            Ai = np.linalg.inv(A)
            beta = np.einsum('spq,sq->sp', Ai, b)
            r = Ys - beta @ X.T
            sign, logdet = np.linalg.slogdet(A)
            ll[idx, g] = -0.5 * ((Ms * np.log(t2 + Vs)).sum(1) + logdet + (W * r * r).sum(1))
            if coef is not None:
                bc[idx, g] = beta[:, coef]
                vc[idx, g] = Ai[:, coef, coef]
    prior = np.full(G, 1.0 / G)
    if ok.any():
        L = ll[ok]
        for _ in range(100):
            a = L + np.log(prior)
            a -= a.max(1, keepdims=True)
            post = np.exp(a)
            post /= post.sum(1, keepdims=True)
            new = post.mean(0) + 1e-12
            new /= new.sum()
            if np.abs(new - prior).max() < 1e-7:
                prior = new
                break
            prior = new
    gi = np.argmax(ll + np.log(prior), axis=1)
    tau2 = np.where(ok, TAU2_GRID[gi], np.nan)
    out = dict(tau2=tau2, k=k, prior=prior)
    if coef is not None:
        rows = np.arange(S)
        beta = np.where(ok, bc[rows, gi], np.nan)
        se = np.sqrt(np.where(ok, vc[rows, gi], np.nan))
        with np.errstate(invalid='ignore', divide='ignore'):
            t = beta / se
            # Wald z: with tau^2 steadied by the shared prior, the normal reference held its size in simulation
            # (4.3-5.4% at alpha 0.05 for 6-20 samples); t(k - p) was far too conservative (0.5% at 6)
            pv = 2 * stats.norm.sf(np.abs(t))
        out.update(beta=beta, se=se, t=t, p=np.where(ok, pv, np.nan))
    return out


def cmh_2xk(strata):
    """Generalized Cochran-Mantel-Haenszel test of a 2 x K table (modified / not, by isoform) over strata
    (samples). strata: list of (m_k, n_k) arrays over the same K isoforms (n_k = 0 where absent).
    Returns (statistic, df, p)."""
    K = len(strata[0][0])
    dev = np.zeros(K)
    cov = np.zeros((K, K))
    for m, n in strata:
        m, n = np.asarray(m, float), np.asarray(n, float)
        N, Mt = n.sum(), m.sum()
        if N < 2 or (n > 0).sum() < 2 or Mt == 0 or Mt == N:
            continue
        e = n * Mt / N
        dev += m - e
        cov += Mt * (N - Mt) / (N * N * (N - 1)) * (N * np.diag(n) - np.outer(n, n))
    keep = np.diag(cov) > 0
    if keep.sum() < 2:
        return np.nan, 0, np.nan
    d, c = dev[keep][:-1], cov[np.ix_(keep, keep)][:-1, :-1]
    try:
        q = float(d @ np.linalg.solve(c, d))
    except np.linalg.LinAlgError:
        q = float(d @ np.linalg.pinv(c) @ d)
    df = int(keep.sum() - 1)
    return q, df, float(stats.chi2.sf(q, df))


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


def r2_on(y, covariate):
    """R^2 and p of y on one covariate (numeric: linear; categorical: one-way ANOVA)."""
    kind, vals = covariate
    ok = np.isfinite(y)
    if kind == 'num':
        v = np.asarray(vals, float)
        ok &= np.isfinite(v)
        if ok.sum() < 4 or np.nanstd(v[ok]) == 0:
            return np.nan, np.nan
        r, p = stats.pearsonr(v[ok], y[ok])
        return r * r, p
    v = np.asarray(vals, dtype=object)
    ok &= np.array([x not in (None, '', '.', 'NA') for x in v])
    levels = sorted(set(v[ok]))
    if len(levels) < 2 or ok.sum() < len(levels) + 2:
        return np.nan, np.nan
    groups = [y[ok & (v == lv)] for lv in levels]
    if any(len(g_) == 0 for g_ in groups):
        return np.nan, np.nan
    yy = y[ok]
    ss_tot = ((yy - yy.mean()) ** 2).sum()
    ss_b = sum(len(g_) * (g_.mean() - yy.mean()) ** 2 for g_ in groups)
    if ss_tot == 0:
        return np.nan, np.nan
    dfb, dfw = len(levels) - 1, len(yy) - len(levels)
    if dfw <= 0:
        return ss_b / ss_tot, np.nan
    f = (ss_b / dfb) / ((ss_tot - ss_b) / dfw) if ss_tot > ss_b else np.inf
    return ss_b / ss_tot, float(stats.f.sf(f, dfb, dfw))


def pca(Z, n=5):
    """PCs of the columns (samples) of Z (features x samples), features centred. Returns (scores, fractions)."""
    Zc = Z - Z.mean(1, keepdims=True)
    _, sv, vt = np.linalg.svd(Zc, full_matrices=False)
    return vt.T[:, :n] * sv[:n], (sv ** 2 / (sv ** 2).sum())[:n]


# ---------------------------------------------------------------- style
def register_fonts(dirs):
    """Add every .ttf/.otf in the given directories to matplotlib, so the style sheet's Helvetica (and its
    bold and oblique faces) resolve without a system install."""
    added = []
    for d in dirs:
        if d and os.path.isdir(d):
            for f in sorted(glob.glob(os.path.join(d, '*.[ot]tf')) + glob.glob(os.path.join(d, '*.[OT]TF'))):
                try:
                    font_manager.fontManager.addfont(f)
                    added.append(os.path.basename(f))
                except Exception as e:                              # a damaged font file: skip it
                    log.warning('font %s not added: %s', f, e)
    return added


def apply_style(path, font_dirs):
    notes = []
    dirs = list(font_dirs or [])
    if path:
        sd = os.path.dirname(os.path.abspath(path))
        dirs += [os.path.join(sd, 'fonts'), os.path.join(os.path.dirname(sd), 'fonts')]
    dirs += FONT_DIRS + [os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fonts')]
    added = register_fonts(list(dict.fromkeys(dirs)))
    if added:
        notes.append(f'fonts registered: {", ".join(added)}')
    if path and os.path.exists(path):
        plt.style.use(path)
        notes.append(f'style: {path}')
    elif path:
        notes.append(f'style {path} not found: matplotlib defaults')
    # v4's overrides, kept: screen dpi for layout (the PDFs are vector), lighter titles for long panel titles
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
    notes.append(f'font: {chosen}' + ('' if chosen == wanted else f' ({wanted} is not installed: --font-dir)'))
    return notes


def figure(rows, cols, widths=None, heights=None):
    w = sum(widths) if widths else cols * PANEL[0]
    h = sum(heights) if heights else rows * PANEL[1]
    fig = plt.figure(figsize=(w + 0.4 * cols, h + 0.5 * rows))
    gs = fig.add_gridspec(rows, cols, width_ratios=widths, height_ratios=heights)
    return fig, gs


def letter(ax, s):
    ax.text(-0.02, 1.02, s, transform=ax.transAxes, fontsize=18, fontweight='bold', ha='right', va='bottom')


def stat_box(ax, lines, loc='upper left', fontsize=10):
    """Statistics as plain text in a corner (the lab style has no stats boxes)."""
    x, ha = (0.04, 'left') if 'left' in loc else (0.96, 'right')
    y, va = (0.96, 'top') if 'upper' in loc else (0.04, 'bottom')
    ax.text(x, y, '\n'.join(lines), transform=ax.transAxes, ha=ha, va=va, fontsize=fontsize)


def empty(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha='center', va='center', fontsize=10, wrap=True)
    ax.set_axis_off()


def sample_axis(ax, samples, batch, colours, fontsize=9, short=False):
    labels = [s[-4:] if short else s for s in samples]
    ax.set_xticks(range(len(samples)), labels, rotation=90, fontsize=fontsize)
    for lab, s in zip(ax.get_xticklabels(), samples):
        lab.set_color(colours[batch[s]])


def batch_legend(ax, colours, fontsize=9, **kw):
    handles = [plt.Line2D([], [], marker='o', ls='', color=c, label=b) for b, c in colours.items()]
    ax.legend(handles=handles, title='batch', fontsize=fontsize, title_fontsize=fontsize + 1, **kw)


def fig_batch_legend(fig, colours, extra=None):
    handles = [plt.Line2D([], [], marker='o', ls='', color=c, label=b) for b, c in colours.items()]
    handles += list(extra or [])
    fig.legend(handles=handles, title=None, loc='lower center', bbox_to_anchor=(0.5, 1.0),
               ncol=len(handles), fontsize=11, handletextpad=0.3, columnspacing=1.2)


def log_axis(axis):
    """Labels at powers of ten, as in the paper; no minor-tick labels. A range of less than ~1.5 decades
    gets 1-2-5 ticks with plain numbers instead, so it is never left with a single label."""
    lo, hi = axis.get_view_interval()
    if lo > 0 and hi / lo < 30:
        axis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        axis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    else:
        axis.set_major_locator(LogLocator(base=10))
    axis.set_minor_formatter(NullFormatter())


def hist2d(ax, x, y, label_x, label_y, bins=40, rng=(0, 100), cmap='viridis', vmax=None, stats_text=True):
    """The paper's 2D histograms: log colour scale, y = x dashed, n / r / RMSE in the corner."""
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
    if stats_text:
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


def ecdf(ax, v, **kw):
    v = np.sort(np.asarray(v, float)[np.isfinite(v)])
    if len(v):
        ax.step(v, np.arange(1, len(v) + 1) / len(v), where='post', **kw)


def weighted_ecdf(ax, vals, counts, **kw):
    if len(vals):
        o = np.argsort(vals)
        cum = np.cumsum(np.asarray(counts, float)[o])
        ax.step(np.asarray(vals)[o], cum / cum[-1], where='post', **kw)


def jitter(n, width=0.25, seed=0):
    return (np.random.default_rng(seed).random(n) - 0.5) * width


def plain(code):
    """ASCII modification name for tables and notes (m6A, m5C, Psi, Ino, Am, Cm, Um, Gm)."""
    return MOD_NAME.get(code, code).replace('$', '').replace('^', '').replace('Ψ', 'Psi')


def mod_axis(ax, codes, fontsize=9):
    ax.set_xticks(range(len(codes)), [MOD_NAME[m] for m in codes], fontsize=fontsize)
    for lab, m in zip(ax.get_xticklabels(), codes):
        lab.set_color(MOD_COLOUR[m])


def is_globin(name):
    """True for a catalogue gene whose name (several, comma-joined, for merged genes) includes a globin gene."""
    return bool(set(str(name).split(',')) & set(GLOBIN))


def is_mito_gene(name):
    return any(p.startswith('MT-') for p in str(name).split(','))


def short_name(name):
    """First name of a merged catalogue gene, with how many more it holds."""
    parts = [p_ for p_ in str(name).split(',') if p_]
    return parts[0] + (f' +{len(parts) - 1}' if len(parts) > 1 else '') if parts else str(name)


def hist_median(vals, cnt):
    cum = np.cumsum(cnt)
    return float(vals[np.searchsorted(cum, cum[-1] / 2)])


# ---------------------------------------------------------------- masks
class SiteMask:
    """Sites kept for modification figures: not on contigs matching the regex (unplaced, random, alt...),
    and outside --mask-bed intervals. Counts what each rule removed from the cohort matrix, per mod."""

    def __init__(self, contig_regex, bed):
        self.pattern = contig_regex or ''
        self.rx = re.compile(contig_regex) if contig_regex else None
        self.iv, raw = {}, {}
        if bed:
            with opener(bed) as fh:
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


def interval_hits(starts_by_chrom, chroms, pos, lo=0, hi=0):
    """For each (chrom, pos): is any sorted position of starts_by_chrom[chrom] within [pos+lo, pos+hi]?"""
    chroms = np.asarray(chroms, dtype=object)
    pos = np.asarray(pos, np.int64)
    out = np.zeros(len(pos), bool)
    for c_ in set(chroms.tolist()):
        p_ = starts_by_chrom.get(c_)
        if p_ is None or not len(p_):
            continue
        sel = np.where(chroms == c_)[0]
        a = np.searchsorted(p_, pos[sel] + lo, side='left')
        b = np.searchsorted(p_, pos[sel] + hi, side='right')
        out[sel] = b > a
    return out


def bed_intervals(path, name_filter=None):
    """chrom -> (sorted starts, ends, values) from a BED; value = column 5 as float when numeric."""
    raw = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            if name_filter and (len(f) < 4 or not name_filter(f[3])):
                continue
            try:
                s_, e_ = int(f[1]), int(f[2])
            except ValueError:
                continue
            try:
                v = float(f[4]) if len(f) > 4 else np.nan
            except ValueError:
                v = np.nan
            st = f[5] if len(f) > 5 and f[5] in '+-' else '.'
            raw.setdefault(f[0], []).append((s_, e_, v, st))
    out = {}
    for c_, rows in raw.items():
        rows.sort()
        a = np.array([(r[0], r[1], r[2]) for r in rows], float)
        out[c_] = (a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2], np.array([r[3] for r in rows]))
    return out


def point_in_bed(iv, chroms, pos, strands=None):
    """Indices into the BED (or -1) for each 0-based position (start <= pos < end), strand-aware when both
    sides have a strand. Works for point BEDs (end = start + 1) and for intervals."""
    chroms = np.asarray(chroms, dtype=object)
    pos = np.asarray(pos, np.int64)
    hit = np.full(len(pos), -1, np.int64)
    val = np.full(len(pos), np.nan)
    for c_ in set(chroms.tolist()):
        if c_ not in iv:
            continue
        s_, e_, v_, st_ = iv[c_]
        sel = np.where(chroms == c_)[0]
        j = np.searchsorted(s_, pos[sel], side='right') - 1
        okj = j >= 0
        jj = np.clip(j, 0, None)
        ok = okj & (pos[sel] < e_[jj])
        if strands is not None:
            sd = np.asarray(strands, dtype=object)[sel]
            ok &= (st_[jj] == '.') | (st_[jj] == sd)
        hit[sel[ok]] = jj[ok]
        val[sel[ok]] = v_[jj[ok]]
    return hit, val


# ---------------------------------------------------------------- reference sequence
class Fasta:
    """Random access to an uncompressed, indexed FASTA (memory-mapped through its .fai)."""

    def __init__(self, path):
        self.ok = bool(path) and os.path.exists(path) and os.path.exists(path + '.fai') and not path.endswith('.gz')
        self.path = path
        if self.ok:
            self.fai = {}
            with open(path + '.fai') as fh:
                for line in fh:
                    f = line.split('\t')
                    self.fai[f[0]] = tuple(int(x) for x in f[1:5])
            self.mm = np.memmap(path, dtype=np.uint8, mode='r')

    def context(self, chroms, starts, strands, flank=2):
        """(n, 2*flank+1) uint8 bases centred on each 0-based start, RNA sense (reverse complement on '-');
        rows off a contig end or on an unknown contig are 'N'."""
        n = len(starts)
        w = 2 * flank + 1
        out = np.full((n, w), ord('N'), np.uint8)
        if not self.ok:
            return out
        chroms = np.asarray(chroms, dtype=object)
        starts = np.asarray(starts, np.int64)
        offs = np.arange(-flank, flank + 1)
        for c_ in set(chroms.tolist()):
            if c_ not in self.fai:
                continue
            length, offset, lb, lw = self.fai[c_]
            sel = np.where(chroms == c_)[0]
            p = starts[sel]
            ok = (p >= flank) & (p + flank < length)
            q = p[ok, None] + offs
            out[sel[ok]] = self.mm[offset + (q // lb) * lw + q % lb]
        low = (out >= 97) & (out <= 122)                            # soft-masked bases: upper case
        out[low] -= 32
        minus = np.asarray(strands, dtype=object) == '-'
        comp = np.arange(256, dtype=np.uint8)
        for a_, b_ in zip(b'ACGTN', b'TGCAN'):
            comp[a_] = b_
        out[minus] = comp[out[minus]][:, ::-1]
        return out


def is_drach(ctx):
    """ctx: (n, 5) uint8 RNA-sense bases centred on the A (T written for U). DRACH: D=A/G/T, R=A/G, A, C, H=A/C/T."""
    D = np.isin(ctx[:, 0], np.frombuffer(b'AGT', np.uint8))
    R = np.isin(ctx[:, 1], np.frombuffer(b'AG', np.uint8))
    A = ctx[:, 2] == ord('A')
    C = ctx[:, 3] == ord('C')
    H = np.isin(ctx[:, 4], np.frombuffer(b'ACT', np.uint8))
    return D & R & A & C & H


# ---------------------------------------------------------------- the cohort
class Cohort:
    def __init__(self, args, mask):
        self.args = args
        self.R, self.ref, self.mask, self.out = args.results, args.ref, mask, args.out
        self.tables = os.path.join(self.out, 'tables')
        os.makedirs(self.tables, exist_ok=True)
        self.legends = {}                                           # figure name -> legend text
        self.summary = {}                                           # sample -> {column: value}, Table 1
        self._conc = None
        self.stream = {}                                            # sample -> cov20 summaries (stream pass)
        self.cohort = os.path.join(self.R, 'cohort', self.ref)
        sheet = pd.read_csv(args.samplesheet, dtype=str, keep_default_na=False)
        sheet.columns = [c.strip() for c in sheet.columns]
        if 'reference' in sheet.columns and (sheet['reference'] == self.ref).any():
            sheet = sheet[sheet['reference'] == self.ref].reset_index(drop=True)
        self.sheet = sheet
        ids = sheet['sample_id'].tolist()
        bcol = args.batch_column
        b = sheet[bcol].tolist() if bcol in sheet.columns else ['all'] * len(ids)
        self.batch = dict(zip(ids, b))
        self.samples = sorted(ids, key=lambda s: (self.batch[s], s))
        names = sorted(set(self.batch.values()))
        self.colours = {n: BATCH_COLOURS[i % len(BATCH_COLOURS)] for i, n in enumerate(names)}
        self.group, self.levels = None, None
        if args.group:
            if args.group not in sheet.columns:
                sys.exit(f'cohort_plots.py: no column {args.group!r} in the samplesheet')
            self.group = dict(zip(ids, sheet[args.group].tolist()))
            lv = args.levels or sorted({v for v in self.group.values() if v not in ('', '.', 'NA')})
            self.levels = lv
        self.gmarker = {lv: GROUP_MARKERS[i % len(GROUP_MARKERS)] for i, lv in enumerate(self.levels or [])}
        self._mods = None
        self._counts = None
        self._expr = None
        self._comp = None
        self._gene_names = None
        self.fasta = Fasta(args.reference)

    def path(self, *parts):
        return os.path.join(self.R, *parts)

    def sample_file(self, kind, s, suffix):
        return self.path(kind, self.ref, s, f'{s}{suffix}')

    def table(self, name, header, rows):
        write_rows(os.path.join(self.tables, name), header, rows)
        return name

    def note_sample(self, s, key, value):
        self.summary.setdefault(s, {})[key] = value

    # -------------------------------------------------- modification matrices
    def mods_matrix(self):
        """Masked cohort sites: (rows, samples, percent, mods, n_valid or None)."""
        if self._mods is None:
            keys = ['chrom', 'start', 'end', 'mod', 'strand']
            pm = read_matrix(os.path.join(self.cohort, f'cohort.{self.ref}.mods.pct_fpcorr.tsv.gz'), keys)
            ms = [s for s in self.samples if s in pm.columns]
            mods = strarr(pm, 'mod')
            keep = self.mask.keep(strarr(pm, 'chrom'), numbers(pm, ['start']).ravel(), mods)
            pct = numbers(pm, ms)[keep]
            nv, nm = None, None
            for metric in ('n_valid', 'n_mod'):
                p_ = os.path.join(self.cohort, f'cohort.{self.ref}.mods.{metric}.tsv.gz')
                if not os.path.exists(p_):
                    continue
                nt = read_matrix(p_, keys)
                if len(nt) == len(pm) and all((strarr(nt, k) == strarr(pm, k)).all() for k in ('chrom', 'start', 'mod', 'strand')):
                    v = numbers(nt, ms)[keep]
                    if metric == 'n_valid':
                        nv = v
                    else:
                        nm = v
                else:
                    log.warning('%s and pct_fpcorr matrices list different sites; not used', metric)
            rows = {k: strarr(pm, k)[keep] for k in ('chrom', 'strand')}
            rows['start'] = numbers(pm, ['start']).ravel()[keep].astype(np.int64)
            self._mods = (rows, ms, pct, mods[keep], nv)
            self._nmod = nm
        return self._mods

    def mods_counts(self):
        """(n_valid, n_mod, fp, clipped) aligned to mods_matrix rows: raw counts, the site's gIVT rate
        fp = (raw% - corrected%) / 100, and where gIVT clipped the corrected value to 0 (fp only a lower
        bound there) unless the stream pass recovered the exact rate."""
        if self._counts is None:
            rows, ms, pct, mods, nv = self.mods_matrix()
            nm = self._nmod
            if nv is None or nm is None:
                self._counts = (None, None, None, None)
            else:
                with np.errstate(invalid='ignore', divide='ignore'):
                    raw = 100 * nm / nv
                    fp = np.clip((raw - pct) / 100, 0, 1)
                clipped = (pct == 0) & (raw > 0)
                fp = np.where(np.isnan(pct), np.nan, fp)
                exact = getattr(self, 'cohort_fp', None)
                if exact is not None and exact.shape == fp.shape:
                    have = np.isfinite(exact)
                    fp = np.where(have, exact, fp)
                    clipped &= ~have
                self._counts = (nv, nm, fp, clipped)
        return self._counts

    def site_keys(self):
        rows, ms, pct, mods, nv = self.mods_matrix()
        return rows, mods

    def depth(self):
        qc = read_table(os.path.join(self.cohort, f'cohort.{self.ref}.read_qc.tsv'))
        other = (numbers(qc, ['non_globin_reads']) - numbers(qc, ['mito_reads'])).ravel()
        return dict(zip(text(qc, 'sample'), other))

    def read_qc(self):
        qc = read_table(os.path.join(self.cohort, f'cohort.{self.ref}.read_qc.tsv'))
        row = {s: i for i, s in enumerate(text(qc, 'sample'))}
        v = {k: numbers(qc, [k]).ravel() for k in ('mapped_reads', 'globin_reads', 'mito_reads', 'non_globin_reads')
             if k in qc.columns}
        return row, v

    # -------------------------------------------------- expression
    def expression(self):
        """(samples, log2(CPM + 1) of genes with median CPM > 1, names, gene ids) with CPM renormalised
        over genes that are neither globin nor mitochondrial: the pipeline's CPM divides by all assigned
        reads, so a sample's globin share would otherwise scale every other gene."""
        if self._expr is None:
            cpm = read_matrix(os.path.join(self.cohort, f'cohort.{self.ref}.gene_cpm.tsv.gz'), ['gene_id', 'gene_name'])
            es = [s for s in self.samples if s in cpm.columns]
            m = np.nan_to_num(numbers(cpm, es))
            ids = strarr(cpm, 'gene_id')
            raw = strarr(cpm, 'gene_name')
            names = np.array([n if n not in ('', '.') else g for g, n in zip(ids, raw)], dtype=object)
            drop = np.array([is_globin(n) or is_mito_gene(n) for n in raw])
            m = m[~drop]
            ids, names = ids[~drop], names[~drop]
            tot = m.sum(0)
            m = m / np.where(tot > 0, tot, 1) * 1e6
            keep = np.median(m, axis=1) > 1
            self._expr_all = (es, names, ids, np.log2(m + 1))
            self._expr = (es, np.log2(m[keep] + 1), list(names[keep]), list(ids[keep]))
            self.n_globin_mito_genes = int(drop.sum())
        return self._expr

    def gene_log(self, gene_names):
        """log2(CPM + 1) rows for the named genes (renormalised CPM; NaN rows for missing genes)."""
        self.expression()
        es, names, ids, lm = self._expr_all
        first = {}
        for i, n in enumerate(names):
            for p_ in str(n).split(','):
                first.setdefault(p_, i)
        return es, np.array([lm[first[g]] if g in first else np.full(len(es), np.nan) for g in gene_names])

    def gene_names(self):
        if self._gene_names is None:
            p = os.path.join(self.cohort, f'cohort.{self.ref}.gene_cpm.tsv.gz')
            self._gene_names = {}
            if os.path.exists(p):
                t = read_table(p, columns=['gene_id', 'gene_name'])
                self._gene_names = {g: n for g, n in zip(text(t, 'gene_id'), text(t, 'gene_name')) if n not in ('.', None)}
        return self._gene_names

    # -------------------------------------------------- composition and covariates
    def composition(self):
        """dict(names, values (samples x types, aligned to self.samples), kind, source)."""
        if self._comp is not None:
            return self._comp
        a = self.args
        S = self.samples
        if a.composition:
            t = pd.read_csv(a.composition, sep=None, engine='python', dtype=str, keep_default_na=False)
            scol = next((c for c in t.columns if c.lower() in ('sample', 'sample_id', 'mixture')), t.columns[0])
            cols = [c for c in t.columns if c != scol and pd.to_numeric(t[c], errors='coerce').notna().mean() > 0.9
                    and c.lower() not in ('p-value', 'correlation', 'rmse', 'absolute score (sig.score)')]
            idx = dict(zip(t[scol], range(len(t))))
            vals = np.array([[pd.to_numeric(t[c], errors='coerce').iloc[idx[s]] if s in idx else np.nan for c in cols] for s in S], float)
            self._comp = dict(names=cols, values=vals, kind='measured', source=os.path.basename(a.composition))
        elif a.composition_columns:
            cols = [c for c in a.composition_columns if c in self.sheet.columns]
            idx = dict(zip(self.sheet['sample_id'], range(len(self.sheet))))
            vals = np.array([[pd.to_numeric(self.sheet[c], errors='coerce').iloc[idx[s]] for c in cols] for s in S], float)
            self._comp = dict(names=cols, values=vals, kind='measured', source='samplesheet')
        else:
            names, vals, used = [], [], {}
            try:
                es, _ = self.gene_log(['XIST'])
            except FileNotFoundError:
                self._comp = dict(names=[], values=np.zeros((len(S), 0)), kind='none', source='no gene_cpm matrix')
                return self._comp
            col = {s: i for i, s in enumerate(es)}
            for ct, genes in CELL_MARKERS.items():
                _, lm = self.gene_log(genes)
                ok = np.isfinite(lm).all(1) & (np.nanstd(lm, 1) > 0)
                if ok.sum() < 2:
                    continue
                z = (lm[ok] - lm[ok].mean(1, keepdims=True)) / lm[ok].std(1, keepdims=True)
                score = z.mean(0)
                names.append(ct)
                vals.append([score[col[s]] if s in col else np.nan for s in S])
                used[ct] = [g for g, k in zip(genes, ok) if k]
            self._comp = dict(names=names, values=np.array(vals, float).T if vals else np.zeros((len(S), 0)),
                              kind='marker score (z)', source='marker genes', markers=used)
        return self._comp

    def covariates(self, extra=True):
        """name -> ('num' | 'cat', values aligned to self.samples). Samplesheet columns (no ids or paths,
        or --covariates), then derived per-sample quantities."""
        out = {}
        sh = self.sheet
        idx = dict(zip(sh['sample_id'], range(len(sh))))
        skip = {'sample_id', 'reference', 'data_type', 'bam', 'bai', 'pod5', 'ubam'}
        cols = self.args.covariates or [c for c in sh.columns if c not in skip]
        for c in cols:
            if c not in sh.columns or c in skip:
                continue
            v = [sh[c].iloc[idx[s]] if s in idx else '' for s in self.samples]
            if any(('/' in x) or x.endswith(('.bam', '.bai', '.pod5', '.fa', '.gz')) for x in v):
                continue
            nonempty = [x for x in v if x not in ('', '.', 'NA')]
            if len(set(nonempty)) < 2:
                continue
            num = pd.to_numeric(pd.Series(v), errors='coerce')
            if num.notna().sum() == len(nonempty) and len(nonempty) >= 4:
                out[c] = ('num', num.to_numpy(float))
            elif len(set(nonempty)) < len(nonempty):                # categorical, not an identifier
                out[c] = ('cat', np.array(v, dtype=object))
        if extra:
            try:
                row, v = self.read_qc()
                for key, lab in (('globin_reads', 'globin %'), ('mito_reads', 'mito %')):
                    if key in v and 'mapped_reads' in v:
                        out[lab] = ('num', np.array([100 * v[key][row[s]] / v['mapped_reads'][row[s]] if s in row else np.nan
                                                     for s in self.samples]))
                d = self.depth()
                out['depth (M)'] = ('num', np.array([d.get(s, np.nan) / 1e6 for s in self.samples]))
            except (FileNotFoundError, KeyError):
                pass
            comp = self.composition()
            for j, nm in enumerate(comp['names']):
                out[nm] = ('num', comp['values'][:, j])
            for key in ('complete_5p', 'speed_median', 'qscore_median', 'read_len_median'):
                vals = np.array([self.summary.get(s, {}).get(key, np.nan) for s in self.samples], float)
                if np.isfinite(vals).sum() >= 4:
                    out[key.replace('_', ' ')] = ('num', vals)
            thr = self.thresholds()
            if thr:
                vals = np.array([thr.get(s, {}).get('a', np.nan) for s in self.samples])
                if np.isfinite(vals).sum() >= 4 and np.nanstd(vals) > 0:
                    out['m6A threshold'] = ('num', vals)
        return out

    def thresholds(self):
        """sample -> {mod code: call threshold}, from each sample's givt_flags.txt (the thresholds gIVT's
        rates and modkit used)."""
        if not hasattr(self, '_thr'):
            self._thr = {}
            for s in self.samples:
                p = self.sample_file('mods', s, '.givt_flags.txt')
                if os.path.exists(p):
                    d = {}
                    for tok in open(p).read().split():
                        if ',' in tok:
                            k, v = tok.split(',', 1)
                            try:
                                d[k] = float(v)
                            except ValueError:
                                pass
                    self._thr[s] = d
        return self._thr

    def snvs(self, s):
        """chrom -> sorted 0-based positions of the sample's PASS known/novel SNVs (editing sites excluded),
        and a set of (chrom, pos0, gt) for concordance. None when the variant table is missing."""
        p = self.sample_file('variants', s, '.variants.tsv.gz')
        if not os.path.exists(p):
            return None
        t = read_table(p, columns=['chrom', 'pos', 'ref', 'alt', 'filter', 'gt', 'dp', 'class'])
        ok = np.ones(len(t), bool)
        if 'filter' in t.columns:
            ok &= strarr(t, 'filter') == 'PASS'
        if 'class' in t.columns:
            ok &= np.isin(strarr(t, 'class'), ['known', 'novel'])
        ref, alt = strarr(t, 'ref'), strarr(t, 'alt')
        ok &= np.array([len(a) == 1 and len(b) == 1 for a, b in zip(ref, alt)])
        t = t[ok]
        pos0 = numbers(t, ['pos']).ravel().astype(np.int64) - 1
        ch = strarr(t, 'chrom')
        by = {}
        for c_, p_ in zip(ch, pos0):
            by.setdefault(c_, []).append(p_)
        by = {k: np.array(sorted(v), np.int64) for k, v in by.items()}
        gts = strarr(t, 'gt') if 'gt' in t.columns else np.array(['.'] * len(t), dtype=object)
        cls = strarr(t, 'class') if 'class' in t.columns else np.array(['.'] * len(t), dtype=object)
        calls = {(c_, int(p_), a_): (g_.replace('|', '/'), k_) for c_, p_, a_, g_, k_ in zip(ch, pos0, alt, gts, cls)}
        return by, calls


# ---------------------------------------------------------------- the cov20_fpcorr stream pass (figures 11-15)
# One pass over each <sample>.cov20_fpcorr.bed (every position x mod code with >= 20 valid reads; columns
# 1-18 modkit bedMethyl, 19 = max(0, percent - 100 fp), empty where gIVT could not correct). Summaries are
# small arrays per modification, cached per sample in <out>/.cache.
STREAM_VERSION = 'v5.2'
NCOV_EDGES = np.array([20, 25, 30, 40, 50, 60, 80, 100, 130, 160, 200, 250, 320, 400, 500, 640, 800, 1000,
                       1300, 1600, 2000, 3000, 5000, 1e12])
NCB = len(NCOV_EDGES) - 1
PCT_EDGES = np.linspace(0, 100, 51)
P_LOG_EDGES = np.concatenate([[0.0], np.logspace(-30, -12, 19),    # bin 0: p < 1e-30; 1 per decade to 1e-12,
                              np.logspace(-12, 0, 121)[1:]])         # then 10 per decade, where BH thresholds fall
NPB = len(P_LOG_EDGES) - 1
FP_EDGES = np.concatenate([[-1e-12], np.logspace(-6, 0, 61)])     # bin 0: fp < 1e-6 (zero)
XT_OFFSETS = np.arange(-10, 11)
THIN_FRACS = (0.1, 0.25, 0.5, 0.75)
BED_COLS = [0, 1, 3, 5, 9, 10, 11, 15, 18]                          # chrom start mod strand Nvalid pct Nmod Nfail corr
NM = len(ALL_MODS)


class GivtTable:
    """gIVT's 9-mer false-positive table, looked up exactly as bin/fp_correct.py does (gIVT's mod_to_col
    column choice for the sample's thresholds, 18-bit 9-mer codes). Used for the rows whose corrected value
    was clipped to 0, where the rate cannot be read back from the file."""
    SHIFTS = (2 * np.arange(8, -1, -1)).astype(np.int64)
    CODE = np.full(256, 255, dtype=np.uint8)
    for _i, _b in enumerate(b'ACGT'):
        CODE[_b] = _i

    def __init__(self, path, thresholds):
        from collections import defaultdict
        thr = defaultdict(lambda: 0.7)
        thr.update(thresholds or {})
        with opener(path) as fh:
            header = fh.readline()
        m_to_c = defaultdict(lambda: 2)
        for mod in thr:                                             # gIVT's mod_to_col, unchanged
            for i, l in enumerate(header.strip().split('\t')):
                if 'false_positive' not in l:
                    continue
                if float(l.split('_')[2]) > thr[mod]:
                    m_to_c[mod] = 2 if i == 2 else i - 1
                    break
        used = sorted({2, *m_to_c.values()})
        tab = pd.read_csv(path, sep='\t', header=None, skiprows=1, usecols=[0, 1, *used], dtype={0: str, 1: str},
                          keep_default_na=False, na_values=['nan', 'NaN'], float_precision='round_trip')
        mers = np.asarray(tab[0].to_numpy(dtype=object), dtype='S9').view(np.uint8).reshape(-1, 9)
        codes = (self.CODE[mers].astype(np.int64) << self.SHIFTS).sum(1)
        mods = tab[1].to_numpy(dtype=object)
        self.index = {}
        self.fps = np.zeros((NM, 4 ** 9))
        self.present = np.zeros((NM, 4 ** 9), bool)
        for code in set(mods.tolist()):
            if code not in MOD_INDEX:
                continue
            rows = mods == code
            self.fps[MOD_INDEX[code], codes[rows]] = tab[m_to_c[code]].to_numpy(float)[rows]
            self.present[MOD_INDEX[code], codes[rows]] = True

    def lookup(self, fasta, chroms, starts, strands, midx):
        ctx = fasta.context(chroms, starts, strands, flank=4)
        v = self.CODE[ctx]
        ok = (v != 255).all(1)
        code = (v.astype(np.int64) << self.SHIFTS).sum(1)
        code[~ok] = 0
        out = np.full(len(starts), np.nan)
        hit = ok & self.present[midx, code]
        out[hit] = np.nan_to_num(self.fps[midx[hit], code[hit]])
        return out


def _stream_arrays():
    z = lambda *s: np.zeros(s, np.int64)                           # noqa: E731
    return dict(rows=z(NM), nvalid=z(NM), nfail=z(NM), correctable=z(NM), clipped=z(NM), exact_fp=z(NM),
                called20=z(NM), funnel=z(NM, NCB, 50), fp_hist=z(NM, len(FP_EDGES) - 1),
                plin=z(NM, 50), plog=z(NM, NCB, 2, NPB), thin=z(NM, len(THIN_FRACS)),
                xt_sum=np.zeros((NM, NM, len(XT_OFFSETS))), xt_n=z(NM, NM, len(XT_OFFSETS)),
                bg_sum=np.zeros(NM), bg_n=z(NM), snp_rows=z(NM), snp_called=z(NM),
                hist_snp=z(NM, 50), hist_nosnp=z(NM, 50), masked_contig=z(NM), masked_bed=z(NM))


def stream_one(job):
    """Summarise one sample's cov20_fpcorr.bed (runs in a worker process). Returns the cache path."""
    s, path, cache, key = job['sample'], job['path'], job['cache'], job['key']
    rng = np.random.default_rng(job['seed'])
    A = _stream_arrays()
    rx = re.compile(job['rx']) if job['rx'] else None
    mask_iv = job['mask_iv']
    ckeys, corder = job['cohort_keys'], job['cohort_order']
    chrom_ids = job['chrom_ids']
    snps = job['snps']
    mt_pos = set(p for p, *_ in MT_RRNA_SITES)
    fasta = Fasta(job['reference']) if job['reference'] else None
    table = None
    if job['errortable'] and fasta is not None and fasta.ok:
        try:
            table = GivtTable(job['errortable'], job['thresholds'])
        except Exception as e:                                      # the table is optional: say why, go on
            log.warning('%s: gIVT table not usable (%s); clipped rows keep lower-bound rates', s, e)
    mito_rows, coh = [], []
    reader = pd.read_csv(path, sep='\t', header=None, usecols=BED_COLS, quoting=csv.QUOTE_NONE,
                         dtype={0: str, 1: np.int64, 3: str, 5: str, 9: np.int64, 10: np.float64, 11: np.int64,
                                15: np.int64, 18: np.float64},
                         keep_default_na=False, na_values={18: ['']}, chunksize=job['chunk'])
    for df in reader:
        ch = df[0].to_numpy(dtype=object)
        md_codes = df[3].to_numpy(dtype=object)
        midx = np.array([MOD_INDEX.get(c_, -1) for c_ in md_codes], np.int64) if len(df) < 1000 else \
            pd.Series(md_codes).map(MOD_INDEX).fillna(-1).to_numpy(np.int64)
        keep = midx >= 0
        if rx is not None:
            u = {c_: bool(rx.search(c_)) for c_ in pd.unique(ch)}
            bad = pd.Series(ch).map(u).to_numpy(bool)
            np.add.at(A['masked_contig'], midx[bad & keep], 1)
            keep &= ~bad
        st = df[1].to_numpy(np.int64)
        if mask_iv:
            inb = np.zeros(len(df), bool)
            for c_, (s_, e_) in mask_iv.items():
                sel = np.where(ch == c_)[0]
                if len(sel):
                    j = np.searchsorted(s_, st[sel], side='right') - 1
                    inb[sel[(j >= 0) & (st[sel] < e_[np.clip(j, 0, None)])]] = True
            np.add.at(A['masked_bed'], midx[inb & keep], 1)
            keep &= ~inb
        if not keep.any():
            continue
        ch, st, md = ch[keep], st[keep], midx[keep]
        minus = (df[5].to_numpy(dtype=object)[keep] == '-')
        n = df[9].to_numpy(np.int64)[keep]
        pct = df[10].to_numpy(float)[keep]
        x = df[11].to_numpy(np.int64)[keep]
        nf = df[15].to_numpy(np.int64)[keep]
        corr = df[18].to_numpy(float)[keep]
        A['rows'] += np.bincount(md, minlength=NM)
        A['nvalid'] += np.bincount(md, weights=n, minlength=NM).astype(np.int64)
        A['nfail'] += np.bincount(md, weights=nf, minlength=NM).astype(np.int64)
        cor = np.isfinite(corr)
        A['correctable'] += np.bincount(md[cor], minlength=NM)
        fp = np.where(cor, np.clip((pct - np.nan_to_num(corr)) / 100, 0, 1), np.nan)
        clip_ = cor & (corr == 0)
        A['clipped'] += np.bincount(md[clip_], minlength=NM)
        known = cor & ~clip_
        if table is not None and clip_.any():
            ci = np.where(clip_)[0]
            ex = table.lookup(fasta, ch[ci], st[ci], np.where(minus[ci], '-', '+'), md[ci])
            got = np.isfinite(ex)
            fp[ci[got]] = ex[got]
            known[ci[got]] = True
            A['exact_fp'] += np.bincount(md[ci[got]], minlength=NM)
        cb = np.clip(np.searchsorted(NCOV_EDGES, n, side='right') - 1, 0, NCB - 1)
        pb = np.clip((np.nan_to_num(corr) / 2).astype(np.int64), 0, 49)
        called = cor & (corr >= 20)
        A['called20'] += np.bincount(md[called], minlength=NM)
        A['funnel'] += np.bincount(((md * NCB + cb) * 50 + pb)[cor], minlength=NM * NCB * 50).reshape(NM, NCB, 50)
        fb = np.clip(np.searchsorted(FP_EDGES, fp[known], side='right') - 1, 0, len(FP_EDGES) - 2)
        A['fp_hist'] += np.bincount(md[known] * (len(FP_EDGES) - 1) + fb, minlength=NM * (len(FP_EDGES) - 1)).reshape(NM, -1)
        # null p-values against Binomial(Nvalid, fp): conventional P(X >= x) for counting, randomized PIT
        # (uniform under the null, despite the discreteness) for the calibration histogram
        pge = np.ones(len(n))
        pr = np.full(len(n), np.nan)
        kk = np.where(known)[0]
        if len(kk):
            with np.errstate(all='ignore'):
                sf = stats.binom.sf(x[kk], n[kk], fp[kk])
                pmf = stats.binom.pmf(x[kk], n[kk], fp[kk])
            pge[kk] = np.clip(sf + pmf, 0, 1)
            pr[kk] = np.clip(sf + rng.random(len(kk)) * pmf, 0, 1)
        # without an exact rate a clipped row has x/n <= fp, so P(X >= x) >= ~0.5: never significant (p = 1)
        lb = np.clip(np.searchsorted(P_LOG_EDGES, pge, side='right') - 1, 0, NPB - 1)
        idx = ((md * NCB + cb) * 2 + called.astype(np.int64)) * NPB + lb
        A['plog'] += np.bincount(idx[cor], minlength=NM * NCB * 2 * NPB).reshape(NM, NCB, 2, NPB)
        okp = np.isfinite(pr)
        A['plin'] += np.bincount(md[okp] * 50 + np.clip((pr[okp] * 50).astype(np.int64), 0, 49),
                                 minlength=NM * 50).reshape(NM, 50)
        # thinning (rarefaction): only rows with >= 4 modified calls can still pass 20 reads and 20%
        th = np.where(cor & (x >= 4) & np.isfinite(fp))[0]
        if len(th):
            for fi, f in enumerate(THIN_FRACS):
                n2 = rng.binomial(n[th], f)
                x2 = rng.hypergeometric(x[th], n[th] - x[th], np.maximum(n2, 0)) if len(th) else np.zeros(0)
                with np.errstate(invalid='ignore', divide='ignore'):
                    c2 = 100 * x2 / np.maximum(n2, 1) - 100 * fp[th]
                ok2 = (n2 >= 20) & (c2 >= 20)
                A['thin'][:, fi] += np.bincount(md[th][ok2], minlength=NM)
        # cross-talk: raw occupancy of every model at -10..+10 nt (RNA sense) around high-occupancy sites
        codes_, _ = pd.factorize(ch)
        pkey = (codes_.astype(np.int64) << 34) | (minus.astype(np.int64) << 33) | st
        tgt = {}
        for t in range(NM):
            sel = np.where(md == t)[0]
            if not len(sel):
                continue
            o = np.argsort(pkey[sel], kind='stable')
            tgt[t] = (pkey[sel][o], pct[sel][o])
            A['bg_sum'][t] += pct[sel].sum()
            A['bg_n'][t] += len(sel)
        anchor = cor & (corr >= job['anchor_min'])
        for a in range(NM):
            ai = np.where(anchor & (md == a))[0]
            if not len(ai):
                continue
            if len(ai) > job['max_anchors']:
                ai = rng.choice(ai, job['max_anchors'], replace=False)
            ak, am = pkey[ai], minus[ai]
            for t, (tk, tp) in tgt.items():
                for d in XT_OFFSETS:
                    q = ak + d
                    j = np.searchsorted(tk, q)
                    jj = np.clip(j, 0, len(tk) - 1)
                    f_ = (j < len(tk)) & (tk[jj] == q)
                    if not f_.any():
                        continue
                    o = np.where(am[f_], -d, d) + 10
                    A['xt_sum'][a, t] += np.bincount(o, weights=tp[jj[f_]], minlength=len(XT_OFFSETS))
                    A['xt_n'][a, t] += np.bincount(o, minlength=len(XT_OFFSETS))
        # donor SNVs inside the 9-mer: the reference 9-mer's gIVT rate is the wrong one for those reads
        if snps:
            flag = interval_hits(snps, ch, st, -4, 4)
            A['snp_rows'] += np.bincount(md[flag], minlength=NM)
            A['snp_called'] += np.bincount(md[flag & called], minlength=NM)
            A['hist_snp'] += np.bincount((md * 50 + pb)[flag & called], minlength=NM * 50).reshape(NM, 50)
            A['hist_nosnp'] += np.bincount((md * 50 + pb)[~flag & called], minlength=NM * 50).reshape(NM, 50)
        # mitochondrial rRNA positive controls
        isM = np.isin(ch, MITO_CONTIGS)
        if isM.any():
            mi = np.where(isM & np.isin(st + 1, list(mt_pos)))[0]
            for i in mi:
                mito_rows.append((int(st[i]) + 1, int(md[i]), int(minus[i]), int(n[i]), int(x[i]), float(pct[i]),
                                  float(corr[i]) if cor[i] else np.nan, float(fp[i]) if np.isfinite(fp[i]) else np.nan))
        # raw counts and exact rates at the cohort's sites
        if len(ckeys):
            cid = pd.Series(ch).map(chrom_ids).fillna(-1).to_numpy(np.int64)
            ok = cid >= 0
            gk = (((cid * 16 + md) * 2 + minus.astype(np.int64)) << 30) | st
            j = np.searchsorted(ckeys, gk)
            jj = np.clip(j, 0, len(ckeys) - 1)
            hit = ok & (j < len(ckeys)) & (ckeys[jj] == gk)
            if hit.any():
                coh.append(np.column_stack([corder[jj[hit]], n[hit], x[hit], np.nan_to_num(fp[hit], nan=-1)]))
    A['mito'] = np.array(mito_rows, float).reshape(-1, 8)
    A['cohort'] = np.vstack(coh) if coh else np.zeros((0, 4))
    tmp = cache + f'.{os.getpid()}.partial.npz'
    np.savez_compressed(tmp, key=np.array(key), **A)
    os.replace(tmp, cache)
    return cache


def stream_pass(c, threads, notes):
    """Run (or load from cache) the stream pass for every sample; fills c.stream and c.cohort_fp."""
    if c.stream:
        return c.stream
    a = c.args
    cache_dir = os.path.join(c.out, '.cache')
    os.makedirs(cache_dir, exist_ok=True)
    try:
        rows, ms, pct, mods, nv = c.mods_matrix()
        chroms = sorted(set(rows['chrom'].tolist()))
        chrom_ids = {k: i for i, k in enumerate(chroms)}
        cid = np.array([chrom_ids[k] for k in rows['chrom']], np.int64)
        mid = np.array([MOD_INDEX.get(m, 15) for m in mods], np.int64)
        gk = (((cid * 16 + mid) * 2 + (rows['strand'] == '-').astype(np.int64)) << 30) | rows['start']
        order = np.argsort(gk, kind='stable')
        ckeys, corder = gk[order], order
    except FileNotFoundError:
        ms, ckeys, corder, chrom_ids, pct = [], np.zeros(0, np.int64), np.zeros(0, np.int64), {}, np.zeros((0, 0))
    mask_iv = {k: (s_.astype(np.int64), e_.astype(np.int64)) for k, (s_, e_) in c.mask.iv.items()}
    sig_base = hashlib.md5((STREAM_VERSION + c.mask.pattern + json.dumps(c.mask.intervals()) +
                            hashlib.md5(ckeys.tobytes()).hexdigest() + str(a.errortable) + str(a.reference) +
                            str(a.seed)).encode()).hexdigest()
    thr = c.thresholds()
    jobs, paths = [], {}
    for s in c.samples:
        p = c.sample_file('mods', s, '.cov20_fpcorr.bed')
        if not os.path.exists(p):
            continue
        st_ = os.stat(p)
        snp = c.snvs(s)
        snps = snp[0] if snp else None
        et = a.errortable if (a.errortable and os.path.exists(a.errortable)) else None
        extra = json.dumps([st_.st_size, int(st_.st_mtime), sorted((thr.get(s) or {}).items()),
                            hashlib.md5(json.dumps({k: v.tolist() for k, v in (snps or {}).items()}).encode()).hexdigest(),
                            et and int(os.stat(et).st_mtime)])
        key = hashlib.md5((sig_base + os.path.abspath(p) + extra).encode()).hexdigest()
        cp = os.path.join(cache_dir, f'{s}.stream.npz')
        paths[s] = cp
        if os.path.exists(cp):
            try:
                if str(np.load(cp)['key']) == key:
                    continue
            except Exception:
                pass
        jobs.append(dict(sample=s, path=p, cache=cp, key=key, rx=c.mask.pattern, mask_iv=mask_iv,
                         cohort_keys=ckeys, cohort_order=corder, chrom_ids=chrom_ids, snps=snps,
                         errortable=et, reference=a.reference if c.fasta.ok else None,
                         thresholds=thr.get(s), seed=int(hashlib.md5(s.encode()).hexdigest()[:8], 16) ^ a.seed,
                         chunk=a.chunk_rows, anchor_min=50.0, max_anchors=5000))
    if jobs:
        log.info('stream pass: %d cov20_fpcorr.bed file(s), %d at a time', len(jobs), max(1, threads))
        if threads > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=min(threads, len(jobs))) as pool:
                list(pool.map(stream_one, jobs))
        else:
            for j in jobs:
                stream_one(j)
    for s, cp in paths.items():
        z = np.load(cp, allow_pickle=False)
        c.stream[s] = {k: z[k] for k in z.files}
    if not c.stream:
        notes.append('stream pass: no <sample>.cov20_fpcorr.bed found; figures 11-15 lose their calibration panels')
        return c.stream
    # exact rates at the cohort's sites, and a check that the stream saw the matrices' counts
    if len(ckeys):
        fpm = np.full(pct.shape, np.nan)
        col = {s: j for j, s in enumerate(ms)}
        mismatch = 0
        nvm = c._mods[4]
        for s, d in c.stream.items():
            if s not in col or not len(d['cohort']):
                continue
            ci = d['cohort'][:, 0].astype(np.int64)
            f_ = d['cohort'][:, 3]
            fpm[ci, col[s]] = np.where(f_ >= 0, f_, np.nan)
            if nvm is not None:
                mismatch += int(np.sum(np.abs(nvm[ci, col[s]] - d['cohort'][:, 1]) > 0))
        c.cohort_fp = fpm
        c._counts = None
        if mismatch:
            notes.append(f'stream pass: {mismatch:,} cohort cells whose n_valid differs from the cov20 file (the matrix wins)')
    tot = {k: sum(int(d[k].sum()) for d in c.stream.values()) for k in ('rows', 'clipped', 'exact_fp')}
    notes.append(f'stream pass: {len(c.stream)} samples, {tot["rows"]:,} rows; {tot["clipped"]:,} rows clipped to 0 by gIVT, '
                 f'exact rate recovered for {tot["exact_fp"]:,}' +
                 ('' if tot['exact_fp'] else ' (give --errortable and --reference for exact null p-values)'))
    for s, d in c.stream.items():
        for code in ('a', '17802', 'm', '17596'):
            i = MOD_INDEX[code]
            tot_calls = d['nvalid'][i] + d['nfail'][i]
            if tot_calls:
                c.note_sample(s, f'fail_frac_{plain(code)}', round(float(d['nfail'][i] / tot_calls), 4))
    return c.stream


def bh_from_hist(counts, q):
    """BH threshold on P(X >= x) from a histogram over P_LOG_EDGES: the largest bin upper edge t with
    #(p <= t) >= t * m / q. Returns (t, rows significant)."""
    counts = np.asarray(counts, float)
    m = counts.sum()
    if m == 0:
        return 0.0, 0
    cum = np.cumsum(counts)
    upper = P_LOG_EDGES[1:]
    ok = cum >= upper * m / q
    if not ok.any():
        return 0.0, 0
    j = np.where(ok)[0].max()
    return float(upper[j]), int(cum[j])


def null_stats(plin):
    """(lambda_GC, pi0) from a 50-bin randomized-PIT histogram."""
    tot = plin.sum()
    if tot == 0:
        return np.nan, np.nan
    frac = plin / tot
    cum = np.cumsum(frac)
    j = int(np.searchsorted(cum, 0.5))                              # median, interpolated within its 2% bin
    below = cum[j - 1] if j > 0 else 0.0
    pmed = (j + (0.5 - below) / max(frac[j], 1e-12)) / 50
    lam = stats.chi2.isf(pmed, 1) / stats.chi2.isf(0.5, 1)
    pi0 = min(1.0, plin[25:].sum() / (0.5 * tot))
    return float(lam), float(pi0)


# ---------------------------------------------------------------- figure 1: read accuracy
def superplot(ax, c, values, ylabel, fmt='{:.4f}'):
    """One point per sample, grouped by batch, with the batch median as a line: the per-sample view next to
    a pooled-read panel (SuperPlots; Lord et al. 2020, J Cell Biol)."""
    batches = sorted(set(c.batch[s] for s in values))
    for i, b in enumerate(batches):
        v = np.array([values[s] for s in values if c.batch[s] == b])
        ax.scatter(i + jitter(len(v), 0.35, i), v, color=c.colours[b], s=22, zorder=3, lw=0)
        if len(v):
            ax.hlines(np.median(v), i - 0.3, i + 0.3, color='black', lw=1.2, zorder=4)
            ax.text(i + 0.32, np.median(v), fmt.format(np.median(v)), fontsize=8, va='center')
    ax.set_xticks(range(len(batches)), batches, fontsize=9, rotation=30)
    for lab, b in zip(ax.get_xticklabels(), batches):
        lab.set_color(c.colours[b])
    ax.set_xlim(-0.6, len(batches) - 0.2)
    ax.set_ylabel(ylabel)


def fig1_reads(c, rs_dir, notes):
    batches = sorted(set(c.batch.values()))
    ident, subs = {b: [] for b in batches}, {b: np.zeros((4, 4)) for b in batches}
    per_med, per_indel, found = {}, {}, 0
    for s in c.samples:
        ip = os.path.join(rs_dir, f'{s}.read_identity.tsv.gz')
        sp = os.path.join(rs_dir, f'{s}.substitutions.tsv')
        if os.path.exists(ip):
            t = read_table(ip, numeric=('identity', 'aligned_len', 'ins_rate', 'del_rate'))
            v = numbers(t, ['identity', 'aligned_len', 'ins_rate', 'del_rate'])
            ident[c.batch[s]].append(v)
            per_med[s] = float(np.nanmedian(v[:, 0]))
            per_indel[s] = float(np.nanmedian(v[:, 2] + v[:, 3]))
            c.note_sample(s, 'identity_median', round(per_med[s], 5))
            found += 1
        if os.path.exists(sp):
            t = read_table(sp, numeric=('count',))
            idx = {b: i for i, b in enumerate('ACGU')}
            for r_, b_, n_ in zip(text(t, 'ref_base'), text(t, 'read_base'), numbers(t, ['count']).ravel()):
                r_, b_ = r_.replace('T', 'U'), b_.replace('T', 'U')
                if r_ in idx and b_ in idx:
                    subs[c.batch[s]][idx[r_], idx[b_]] += n_
    if not found:
        notes.append(f'figure 1 skipped: no read_stats outputs in {rs_dir} (run scripts/read_stats.sh per sample)')
        return None
    ident = {b: np.vstack(v) for b, v in ident.items() if v}
    nb = len(ident)
    ncol = max(4, nb + 1)
    fig, gs = figure(2, ncol)
    ax = fig.add_subplot(gs[0, 0])                                  # A identity, pooled per batch
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
    ax = fig.add_subplot(gs[0, 1])                                  # B per-sample median identity
    superplot(ax, c, per_med, 'median identity, per sample')
    letter(ax, 'B')
    ax = fig.add_subplot(gs[0, 2])                                  # C indels, pooled per batch
    grid = np.linspace(0, 0.10, 201)
    top = 0
    for b, v in ident.items():
        h, e = np.histogram(v[:, 2] + v[:, 3], bins=grid, density=True)
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
    ax = fig.add_subplot(gs[0, 3])                                  # D per-sample median indel rate
    superplot(ax, c, per_indel, 'median indels per base, per sample')
    letter(ax, 'D')

    ax = fig.add_subplot(gs[1, 0])                                  # E identity by aligned length
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
    letter(ax, 'E')
    for j, b in enumerate([b for b in batches if b in ident]):     # F... substitution matrices
        ax = fig.add_subplot(gs[1, j + 1])
        m = subs[b]
        rate = 100 * m / np.maximum(m.sum(1, keepdims=True), 1)
        show = np.log10(np.clip(rate, 0.01, 100))                    # log colour: the matrix spans 0.1-99%
        ax.pcolormesh(np.arange(5) - 0.5, np.arange(5) - 0.5, show, cmap=MATRIX_CMAP, vmin=-1.3, vmax=2,
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
        letter(ax, 'FGHIJKLM'[j])
    fig_batch_legend(fig, {b: c.colours[b] for b in ident})
    nreads = sum(len(v) for v in ident.values())
    notes.append(f'figure 1: {nreads:,} sampled reads (>= 200 aligned bases), {found} samples')
    c.legends['Fig1_read_qc'] = f"""Read-level accuracy, as the paper's Fig. 1 ({found} samples, {nreads:,} reads).
Reads are a random 2% of each sample's primary alignments with MAPQ >= 40 and >= 200 aligned bases
(scripts/read_stats.sh). A, C POOLED within each batch (every read counts once, so deeper samples weigh more);
B, D PER SAMPLE (one point per sample, line = batch median), which is the unit for any comparison.
A, identity, matches / (matches + mismatches + insertions + deletions); dashed lines and numbers: medians.
C, insertions plus deletions per aligned base. E, identity by aligned length, boxes per batch (median,
quartiles; whiskers 1.5 IQR). F onwards, substitution matrix per batch: % of reference bases (rows, RNA sense)
read as each base (columns), pooled over positions of the read_stats.sh region with >= 20 reads."""
    return fig


# ---------------------------------------------------------------- figure 2: cohort overview
def clustered_order(r):
    d = np.clip(1 - r, 0, None)
    np.fill_diagonal(d, 0)
    d = (d + d.T) / 2
    if len(r) < 3:
        return np.arange(len(r))
    return leaves_list(linkage(squareform(d, checks=False), 'average'))


def fig2_overview(c, notes):
    row, v = c.read_qc()
    ss = [s for s in c.samples if s in row]
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
    ax.set_ylabel('non-globin, non-mito\nreads (M)')
    sample_axis(ax, ss, c.batch, c.colours)
    letter(ax, 'B')

    ax = fig.add_subplot(gs[0, 2])                                  # C sites against depth
    try:
        _, mss, mpct, _, _ = c.mods_matrix()
        sites = dict(zip(mss, (np.nan_to_num(mpct) >= 20).sum(0).astype(float)))
    except FileNotFoundError:
        sites = {}
    ys = np.array([sites.get(s, np.nan) for s in ss]) / 1e3
    ok = ~np.isnan(ys)
    ax.scatter(other[ok] / 1e6, ys[ok], c=[c.colours[c.batch[s]] for s, k in zip(ss, ok) if k], s=30, zorder=3)
    for s, a_, b_ in zip(ss, other / 1e6, ys):
        if not np.isnan(b_):
            ax.annotate(s[-2:], (a_, b_), fontsize=7, xytext=(3, 2), textcoords='offset points')
    if ok.sum() >= 3:
        k, b0 = np.polyfit(other[ok] / 1e6, ys[ok], 1)
        xx = np.linspace(0, (other[ok] / 1e6).max() * 1.05, 10)
        ax.plot(xx, k * xx + b0, ls='--', color='0.5', lw=0.8)
        stat_box(ax, [f'n = {int(ok.sum())}', f'r = {np.corrcoef(other[ok], ys[ok])[0, 1]:.3f}'])
    ax.set_xlabel('non-globin reads (M)')
    ax.set_ylabel('modification sites (k)')
    letter(ax, 'C')

    es, lx, _, _ = c.expression()
    ax = fig.add_subplot(gs[1, 0:2])                                # D PCA
    keep = np.argsort(lx.var(axis=1))[::-1][:min(2000, len(lx))]
    pcs, frac = pca(lx[keep], 2)
    for i, s in enumerate(es):
        mk = c.gmarker.get(c.group.get(s), 'o') if c.group else 'o'
        ax.scatter(pcs[i, 0], pcs[i, 1], color=c.colours[c.batch[s]], s=45, zorder=3, marker=mk)
        ax.annotate(s, (pcs[i, 0], pcs[i, 1]), fontsize=8, xytext=(4, 3), textcoords='offset points')
    ax.set_xlabel(f'PC1 ({100 * frac[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({100 * frac[1]:.1f}%)')
    letter(ax, 'D')

    ax = fig.add_subplot(gs[1, 2])                                  # E correlation matrix, clustered
    r = np.corrcoef(lx, rowvar=False)
    o = clustered_order(r)
    r = r[np.ix_(o, o)]
    eo = [es[i] for i in o]
    off = r[~np.eye(len(r), dtype=bool)]
    shown = r.copy()
    np.fill_diagonal(shown, np.nan)
    cmap = MATRIX_CMAP.copy()
    cmap.set_bad('#FBF3EE')
    n = len(es)
    im = ax.pcolormesh(np.arange(n + 1) - 0.5, np.arange(n + 1) - 0.5, np.ma.masked_invalid(shown), cmap=cmap,
                       vmin=np.nanmin(off), vmax=np.nanmax(off), edgecolors='face', linewidth=0.3)
    if n <= 12:
        for i in range(n):
            for j in range(n):
                if i != j:
                    ax.text(j, i, f'{r[i, j]:.3f}'.lstrip('0'), ha='center', va='center', fontsize=5.5,
                            color='black' if r[i, j] > (np.nanmin(off) + np.nanmax(off)) / 2 else 'white')
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect('equal')
    sample_axis(ax, eo, c.batch, c.colours, fontsize=7)
    ax.set_yticks(range(n), eo, fontsize=7)
    for lab, s in zip(ax.get_yticklabels(), eo):
        lab.set_color(c.colours[c.batch[s]])
    ax.tick_params(length=0)
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(False)
    colourbar(fig, im, ax, label='Pearson r, log$_2$(CPM + 1)')
    letter(ax, 'E')
    extra = [plt.Line2D([], [], marker=c.gmarker[lv], ls='', color='0.3', label=lv) for lv in (c.levels or [])]
    fig_batch_legend(fig, c.colours, extra)
    notes.append(f'figure 2: {lx.shape[0]:,} genes with median CPM > 1 after removing {c.n_globin_mito_genes} globin and '
                 'chrM genes and renormalising CPM')
    for i, s in enumerate(ss):
        c.note_sample(s, 'mapped_reads', int(g['mapped_reads'][i]))
        c.note_sample(s, 'globin_pct', round(100 * g['globin_reads'][i] / g['mapped_reads'][i], 1))
        c.note_sample(s, 'mito_pct', round(100 * g['mito_reads'][i] / g['mapped_reads'][i], 1))
        c.note_sample(s, 'non_globin_non_mito_reads', int(other[i]))
    gcp = os.path.join(c.cohort, f'cohort.{c.ref}.gene_counts.tsv.gz')
    if os.path.exists(gcp):
        gct = read_matrix(gcp, ['gene_id', 'gene_name'])
        gs_ = [s for s in ss if s in gct.columns]
        for s, d in zip(gs_, (np.nan_to_num(numbers(gct, gs_)) >= 10).sum(0)):
            c.note_sample(s, 'genes_ge10_reads', int(d))
    c.legends['Fig2_cohort_overview'] = f"""Cohort overview, one value per sample throughout ({len(ss)} samples, coloured by batch{', marker = ' + c.args.group if c.group else ''}).
A, mapped reads in globin genes, on chrM, and everything else, as % of the sample's mapped reads; the number
above each bar is the sample's mapped reads (millions). B, mapped reads outside globin genes and chrM: the
depth the rest of the analysis works with; dashed: cohort median. C, modification sites per sample (all eight
modifications; >= 20 reads and >= 20% after gIVT correction; masked regions left out) against B; dashed:
least-squares line. Site counts are mostly a depth readout (figure 14 thins each sample to show how far from
saturation it is). D, PCA of log2(CPM + 1) over the 2,000 most variable of {lx.shape[0]:,} genes with median CPM
> 1. CPM is renormalised over genes that are neither globin nor mitochondrial ({c.n_globin_mito_genes} removed):
the pipeline's CPM divides by all assigned reads, so a sample's globin share would scale every other gene and
become PC1. E, Pearson r between samples over the same genes, hierarchically clustered (average linkage on
1 - r){'; values printed' if n <= 12 else ''}. Per-sample numbers: tables/Table1_cohort_summary.tsv."""
    return fig


# ---------------------------------------------------------------- figure 3: gIVT correction
# Each <sample>.cov20_fpcorr.bed has a row for every position x mod code with >= 20 reads: tens of millions of
# rows. awk streams each one into 2%-wide bins of Dorado percent (column 11) x corrected percent (column 19;
# -1 where gIVT could not correct), applying the same masks as the other modification figures. Bins are cached
# in <out>/.cache, keyed by the file's size and time and by the masks, so a rerun only reads changed files.
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
        p = c.sample_file('mods', s, '.cov20_fpcorr.bed')
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


def paired_points(ax, i, before, after, colour, w=0.18):
    """Per-sample before -> after as two columns of points joined by grey lines, around x = i."""
    before, after = np.maximum(np.asarray(before, float), 0.8), np.maximum(np.asarray(after, float), 0.8)
    for a_, b_ in zip(before, after):
        ax.plot([i - w, i + w], [a_, b_], color='0.75', lw=0.5, zorder=1)
    ax.scatter(np.full(len(before), i - w), before, s=12, facecolor='white', edgecolor=colour, lw=0.8, zorder=3)
    ax.scatter(np.full(len(after), i + w), after, s=12, color=colour, lw=0, zorder=3)
    if len(before):
        ax.hlines([np.median(before), np.median(after)], [i - w - 0.1, i + w - 0.1], [i - w + 0.1, i + w + 0.1],
                  color='black', lw=1.2, zorder=4)


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
        for code in ALL_MODS:
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
    for k, code in enumerate(ALL_MODS):
        ax = fig.add_subplot(gs[k // 4, k % 4])
        axes.append(ax)
        h = np.ma.masked_equal(hists[code].T, 0)
        if h.count():
            last = ax.pcolormesh(edges, edges, h, cmap='viridis', norm=LogNorm(vmin=1, vmax=vmax), edgecolors='face', linewidth=0.2)
        ax.axvline(20, color='#E00000', ls='--', lw=0.9)
        ax.axhline(20, color='#E00000', ls='--', lw=0.9)
        hh = hists[code]
        kept, removed, raised = hh[10:, 10:].sum(), hh[10:, :10].sum(), hh[:10, 10:].sum()
        ax.text(23, 97, f'kept {kept:,.0f}', ha='left', va='top', fontsize=8)
        ax.text(98, 2, f'removed {removed:,.0f}', ha='right', va='bottom', fontsize=8)
        if raised:
            ax.text(2, 98, f'raised {raised:,.0f}', ha='left', va='top', fontsize=8)
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
    codes = [m for m in ALL_MODS if any(r[1] == plain(m) for r in table)]
    for i, code in enumerate(codes):
        before = np.array([r[3] for r in table if r[1] == plain(code)])
        after = np.array([r[4] for r in table if r[1] == plain(code)])
        paired_points(ax, i, before, after, MOD_COLOUR[code])
        drop = 100 * (1 - np.median(after) / max(np.median(before), 1))
        ax.text(i, max(before.max(), 1) * 1.8, f'−{drop:.0f}%', ha='center', fontsize=9)
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    mod_axis(ax, codes)
    ax.set_xlim(-0.6, len(codes) - 0.4)
    ax.set_ylabel('sites per sample (≥ 20%)')
    ax.set_title('per sample: before (open) and after (filled) gIVT subtraction; lines = medians', fontsize=11)
    letter(ax, 'I')
    unc = sum(r[8] for r in table)
    notes.append(f'figure 3: {len(bins)} samples; {unc:,} site rows >= 20% that gIVT could not correct are left out of A-H '
                 f'(and out of every final site set); per-sample counts in tables/{tname}')
    c.legends['Fig3_ivt_correction'] = f"""Effect of gIVT false-positive subtraction, as the paper's Fig. 3 and Table 4 ({len(bins)} samples).
A-H, POOLED over samples: every position x modification with >= 20 reads in a sample is one entry, so a position
covered in all {len(bins)} samples appears {len(bins)} times and deeper samples contribute more. x: Dorado occupancy
(modkit percent modified at the sample's thresholds); y: max(0, x - 100 x the gIVT false-positive rate of its 9-mer).
2%-wide bins, log colour scale shared by all panels; red dashed lines at 20%. kept: >= 20% before and after;
removed: >= 20% before, < 20% after. {unc:,} entries >= 20% that gIVT could not correct (no 9-mer, base
mismatch, 9-mer not in the table) are not drawn; they are not in any final site set either.
I, PER SAMPLE: sites >= 20% before (open) and after (filled) correction, one pair of points per sample joined by
a line; black ticks: medians; the number is the change in the median. Per-sample counts: tables/{tname}."""
    for kind in ('contig', 'bed'):
        parts = [f'{plain(m)} {masked[(kind, m)]:,}' for m in ALL_MODS if (kind, m) in masked]
        if parts:
            notes.append(f'figure 3, masked by {"contig" if kind == "contig" else "--mask-bed"}: site rows ' + ', '.join(parts))
    if short:
        notes.append(f'figure 3: {short:,} rows with fewer than 19 columns skipped')
    return fig


# ---------------------------------------------------------------- figure 4: modification calls
def gene_sites(c):
    """mod code -> {gene_id: set of (chrom, start, strand)}, union over samples, masked."""
    out = {}
    for s in c.samples:
        p = c.sample_file('mods', s, '.gene_mods.bed')
        if not os.path.exists(p):
            continue
        t = read_table(p, columns=['chrom', 'start', 'mod_code', 'gene_id', 'strand'], numeric=('start',))
        if not len(t):
            continue
        ch, st, codes = strarr(t, 'chrom'), numbers(t, ['start']).ravel(), strarr(t, 'mod_code')
        keep = c.mask.keep(ch, st, None)
        for code, g, c_, s_, sd, k in zip(codes, strarr(t, 'gene_id'), ch, st, strarr(t, 'strand'), keep):
            if k:
                out.setdefault(code, {}).setdefault(g, set()).add((c_, int(s_), sd))
    return out


def isolated(chroms, starts, strands, mods, min_gap=5):
    """True for sites with no site of a different modification within min_gap nt on the same chrom and
    strand (adjacent calls by two models can echo one another)."""
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


def strip_points(ax, xpos, vals, colour, width=0.18, s=8, floor=0.8):
    """One point per sample at x = xpos (jittered within width), black tick at the median."""
    vals = np.asarray(vals, float)
    if floor is not None:
        vals = np.maximum(vals, floor)
    vals = vals[np.isfinite(vals)]
    if not len(vals):
        return
    ax.scatter(xpos + jitter(len(vals), 0.6 * width, int(abs(xpos) * 7) % 97), vals, s=s, color=colour, lw=0,
               alpha=0.85, zorder=3)
    ax.hlines(np.median(vals), xpos - width / 2, xpos + width / 2, color='black', lw=1.0, zorder=4)


def fig4_mods(c, notes):
    rows, ms, pct, mods, nv = c.mods_matrix()
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
            for bi in range(len(OCC_BINS)):
                strip_points(ax, bi - 0.4 + w * (j + 0.5), counts[bi], MOD_COLOUR[code], width=w * 0.9, s=6)
            ax.scatter([], [], color=MOD_COLOUR[code], label=MOD_NAME[code], s=12)
        ax.set_yscale('log')
        log_axis(ax.yaxis)
        ax.set_xticks(range(len(OCC_BINS)), [str(b) for b in OCC_BINS])
        ax.set_xlabel('modification occupancy (%, bin start)')
        ax.set_ylabel('sites per sample')
        ax.legend(fontsize=9, ncol=2, loc='upper right', handlelength=1)
        letter(ax, lab)

    ax = fig.add_subplot(gs[0, 2])                                  # C consistency among measuring samples
    meas = np.isfinite(pct) if nv is None else (np.nan_to_num(nv) >= 20) & np.isfinite(pct)
    called = np.nan_to_num(pct) >= 20
    need = max(3, math.ceil(0.75 * len(ms)))
    fb = np.linspace(0, 1, 11)
    for code in [m for m in ALL_MODS if (mods == m).any()]:
        sel = (mods == code) & (meas.sum(1) >= need)
        if sel.sum() < 5:
            continue
        frac = (called[sel] & meas[sel]).sum(1) / meas[sel].sum(1)
        h, _ = np.histogram(frac, bins=fb)
        ax.plot((fb[:-1] + fb[1:]) / 2, 100 * h / h.sum(), marker='o', ms=3.5, lw=1, color=MOD_COLOUR[code],
                ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_xlabel('fraction of measuring samples calling it')
    ax.set_ylabel('sites (%)')
    ax.set_title(f"sites measured in ≥ {need} of {len(ms)}; 2'-O-Me dashed", fontsize=9)
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

    ax = fig.add_subplot(gs[1, 2])                                  # E co-occurrence vs expression-matched null
    allsites = [(m, g, s_) for m in MAIN_MODS if m in genes for g, ss in genes[m].items() for s_ in ss]
    if allsites:
        iso = isolated([s_[0] for _, _, s_ in allsites], [s_[1] for _, _, s_ in allsites], [s_[2] for _, _, s_ in allsites],
                       [m for m, _, _ in allsites])
        gsets = {}
        for (m, g, _), ok in zip(allsites, iso):
            if ok:
                gsets.setdefault(m, set()).add(g)
        # strata: deciles of mean log2 CPM over the expressed universe (genes with median CPM > 1 or a site)
        try:
            c.expression()
            es_, names_, ids_, lm_ = c._expr_all
            mean_expr = dict(zip(ids_, lm_.mean(1)))
        except FileNotFoundError:
            mean_expr = {}
        universe = {g for g, v in mean_expr.items() if v >= 1} | {g for v in gsets.values() for g in v}
        ev = np.array([mean_expr.get(g, 0.0) for g in sorted(universe)])
        cuts = np.quantile(ev, np.linspace(0, 1, 11)[1:-1]) if len(ev) > 20 else []
        stratum = dict(zip(sorted(universe), np.searchsorted(cuts, ev)))
        present = [m for m in MAIN_MODS if gsets.get(m)]
        if len(present) >= 2:
            w = 0.8 / (len(present) - 1)
            for i, obs in enumerate(present):
                for j, giv in enumerate([g for g in present if g != obs]):
                    both = len(gsets[obs] & gsets[giv])
                    exp_ = 0.0
                    for st_ in set(stratum.values()):
                        in_s = {g for g, v in stratum.items() if v == st_}
                        n_s = len(in_s)
                        if n_s:
                            exp_ += len(gsets[giv] & in_s) * len(gsets[obs] & in_s) / n_s
                    pobs, pexp = both / len(gsets[giv]), exp_ / len(gsets[giv])
                    co_rows.append((plain(obs), plain(giv), len(gsets[giv]), both, pobs, pexp, pobs / pexp if pexp else np.nan))
                    xpos = i - 0.4 + w * (j + 0.5)
                    ax.bar(xpos, pobs, width=w, color=MOD_COLOUR[giv])
                    ax.hlines(pexp, xpos - w / 2, xpos + w / 2, color='black', lw=1.4)
                    ax.text(xpos, 0.03, f'{both:,}', rotation=90, ha='center', va='bottom', fontsize=7, color='white')
            ax.legend([plt.Rectangle((0, 0), 1, 1, color=MOD_COLOUR[m]) for m in present] +
                      [plt.Line2D([], [], color='black', lw=1.4)], [MOD_NAME[m] for m in present] + ['expected'],
                      title='given', fontsize=9, title_fontsize=10, loc='upper left', bbox_to_anchor=(1.0, 1.0), handlelength=1)
            ax.set_xticks(range(len(present)), [MOD_NAME[m] for m in present])
            ax.set_ylim(0, 1.0)
        notes.append(f'figure 4E: {int(iso.sum()):,} of {len(iso):,} gene sites have no other modification within 5 nt')
    ax.set_xlabel('observed modification')
    ax.set_ylabel('P(observed | given)')
    ax.set_title('bars: observed; black: expected within expression deciles', fontsize=9)
    letter(ax, 'E')
    fig_batch_legend(fig, c.colours)
    t1 = c.table('Fig4_sites_per_occupancy_bin.tsv', ['sample', 'mod', 'occupancy_bin', 'sites'], bin_rows)
    t2 = c.table('Fig4_cooccurrence.tsv', ['observed', 'given', 'genes_with_given', 'genes_with_both', 'p_observed_given',
                                           'p_expected_expression_matched', 'observed_over_expected'], co_rows)
    ngenes = len({g for m in genes for g in genes[m]})
    c.legends['Fig4_modification_calls'] = f"""Modification calls in the cohort, as the paper's Fig. 5A, B, D ({len(ms)} samples, {pct.shape[0]:,} cohort sites).
A, B, PER SAMPLE: sites with >= 20 reads and >= 20% gIVT-corrected occupancy in that sample, in 10% occupancy
bins (bin start shown; 90 is 90-100%); one point per sample, black tick = median (A, m6A, m5C, Ψ, inosine;
B, 2'-O-methylations). Points, not bars: a bar on a log axis has no meaningful baseline. tables/{t1}.
C, among sites measured (>= 20 reads, corrected value) in at least {need} of {len(ms)} samples, the fraction of
those measuring samples that call the site: constitutive sites sit at 1, sampling noise near the threshold
spreads sites across the middle. Unlike counting calling samples, this does not mix depth into the answer.
D, POOLED (union over samples): gene-level sites from each sample's gene_mods.bed, each site counted once per
gene however many samples call it; genes by number of sites (> 50 in the last bar), {ngenes:,} genes.
E, POOLED: of the genes carrying a site of the given modification (bar colour), the fraction that also carry the
observed modification (x); a site counts only if no site of another modification lies within 5 nt. Black lines:
the fraction expected if the two were independent within deciles of gene expression: highly expressed genes
carry sites of every kind, so the raw fraction (the paper's Fig. 5D) mostly measures expression. The 5-nt rule
also removes true neighbours (mt 16S Um1369/Gm1370 are adjacent; figure 13 measures the echo directly).
tables/{t2}."""
    return fig


# ---------------------------------------------------------------- figure 5 (ICC) and S2 (Jaccard): agreement between donors
def jaccard_matrix(pct, nv, min_cov=20, min_pct=20):
    """Pairwise Jaccard of called sites among sites both samples measure (>= min_cov reads and a corrected
    percent). Returns (Jaccard, sites called in either)."""
    k = pct.shape[1]
    cov = (np.nan_to_num(nv) >= min_cov) & ~np.isnan(pct)
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


def pairwise_icc(pct, nv, selection='others', min_cov=20, min_pct=20, min_sites=10):
    """ICC(2,1) for every pair of samples over the sites both measure (>= min_cov reads and a corrected value
    in each). selection 'others' (default): sites called (>= min_pct) in at least one OTHER sample, so the
    choice of sites does not depend on the pair's own values; 'either': called in either of the pair (v4,
    conditions on the outcome). Returns (ICC, sites used); NaN below min_sites."""
    k = pct.shape[1]
    cov = (np.nan_to_num(nv) >= min_cov) & ~np.isnan(pct)
    called = (np.nan_to_num(pct) >= min_pct) & cov
    ncalled = called.sum(1)
    out, used = np.full((k, k), np.nan), np.zeros((k, k), int)
    for i in range(k):
        for j in range(i + 1, k):
            both = cov[:, i] & cov[:, j]
            if selection == 'either':
                sel = both & (called[:, i] | called[:, j])
            else:
                sel = both & ((ncalled - called[:, i] - called[:, j]) >= 1)
            used[i, j] = used[j, i] = int(sel.sum())
            if sel.sum() >= min_sites:
                out[i, j] = out[j, i] = icc21(pct[sel][:, [i, j]])
    return out, used


def concordance(c):
    """Pairwise ICC and Jaccard per modification, cohort ICC, and each pair's smaller depth; computed once
    for figures 5 and S2."""
    if c._conc is not None:
        return c._conc
    rows, ms, pct, mods, nv = c.mods_matrix()
    if nv is None:
        c._conc = {}
        return c._conc
    n = len(ms)
    sel = c.args.icc_selection
    groups = np.array([c.batch[s] for s in ms])
    depth = c.depth()
    d = np.array([depth.get(s, np.nan) for s in ms]) / 1e6
    out = {'ms': ms, 'n': n, 'groups': groups, 'same': groups[:, None] == groups[None, :], 'offd': ~np.eye(n, dtype=bool),
           'mindepth': np.minimum(d[:, None], d[None, :]), 'icc': {}, 'icc_n': {}, 'jac': {}, 'jac_n': {}, 'icc_all': {}}
    for code in ALL_MODS:
        s_ = mods == code
        if s_.sum() < 3:
            continue
        out['icc'][code], out['icc_n'][code] = pairwise_icc(pct[s_], nv[s_], sel)
        out['jac'][code], out['jac_n'][code] = jaccard_matrix(pct[s_], nv[s_])
        full = ((np.nan_to_num(nv[s_]) >= 20) & ~np.isnan(pct[s_])).all(1) & (np.nan_to_num(pct[s_]) >= 20).any(1)
        out['icc_all'][code] = (icc21(pct[s_][full]) if full.sum() >= 3 else np.nan, int(full.sum()))
    summary, pairs = [], []
    same, offd = out['same'], out['offd']
    w = lambda m, msk: float(np.nanmean(m[msk])) if np.isfinite(m[msk]).any() else np.nan   # noqa: E731
    for code in out['icc']:
        im, jm = out['icc'][code], out['jac'][code]
        summary.append((plain(code), out['icc_all'][code][1], out['icc_all'][code][0], w(im, offd), w(im, same & offd),
                        w(im, ~same), w(jm, offd), w(jm, same & offd), w(jm, ~same)))
        for a in range(n):
            for b in range(a + 1, n):
                pairs.append((plain(code), ms[a], ms[b], int(same[a, b]), out['mindepth'][a, b], im[a, b],
                              int(out['icc_n'][code][a, b]), jm[a, b], int(out['jac_n'][code][a, b])))
    out['t1'] = c.table('Fig5_between_donor_agreement.tsv',
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
    for k, code in enumerate(ALL_MODS):
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
            ax.set_title(title_fmt.format(name=MOD_NAME[code], mean=np.nanmean(m) if np.isfinite(m).any() else np.nan),
                         color=MOD_COLOUR[code], fontsize=12, fontweight='bold')
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
    codes = [m for m in ALL_MODS if m in mats]
    for i, code in enumerate(codes):
        m = mats[code]
        wv = np.nanmean(m[same & offd]) if np.isfinite(m[same & offd]).any() else np.nan
        bv = np.nanmean(m[~same]) if np.isfinite(m[~same]).any() else np.nan
        ax.plot([i - 0.15, i + 0.15], [wv, bv], color='0.6', lw=0.8, zorder=1)
        ax.scatter([i - 0.15], [wv], color=MOD_COLOUR[code], marker='o', s=36, zorder=3)
        ax.scatter([i + 0.15], [bv], facecolor='white', edgecolor=MOD_COLOUR[code], marker='o', s=36, zorder=3, lw=1.2)
    mod_axis(ax, codes)
    ax.set_ylabel(ylabel)
    ax.set_title('same batch (filled), different batches (open)', fontsize=9)


def depth_check(ax, c, conc, mats, ylabel):
    """One point per sample pair: its agreement averaged over the modifications, against the pair's smaller
    depth. Same-batch pairs filled in the batch colour, different-batch pairs open."""
    ms, n, same = conc['ms'], conc['n'], conc['same']
    stack = np.array([mats[m] for m in mats])
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(stack, 0)
    xs, ys = [], []
    for a in range(n):
        for b in range(a + 1, n):
            if np.isfinite(mean[a, b]) and np.isfinite(conc['mindepth'][a, b]):
                xs.append(conc['mindepth'][a, b])
                ys.append(mean[a, b])
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
    codes = [m for m in ALL_MODS if m in ia and not np.isnan(ia[m][0])]
    ax.bar(range(len(codes)), [ia[m][0] for m in codes], color=[MOD_COLOUR[m] for m in codes], width=0.7)
    for i, m in enumerate(codes):
        ax.text(i, ia[m][0] + 0.01, f'{ia[m][0]:.2f}', ha='center', va='bottom', fontsize=7)
        ax.text(i, 0.02, f'{ia[m][1]:,}', rotation=90, ha='center', va='bottom', fontsize=6, color='white')
    mod_axis(ax, codes, fontsize=8)
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
    rule = ('called (>= 20%) in at least one OTHER sample (leave-pair-out)' if c.args.icc_selection == 'others'
            else '>= 20% in either sample of the pair (v4 rule: selects on the outcome)')
    c.legends['Fig5_between_donor_agreement'] = f"""Agreement of modification occupancy BETWEEN DONORS ({n} samples). Different donors are not replicates:
this agreement mixes measurement error with real donor differences. Figure 14 separates the two (split-half
reliability at matched depth). A-H, ICC(2,1) for every PAIR of samples (two-way random effects, absolute
agreement, single measurement: a sample that reads systematically higher counts as disagreement, so a batch
offset lowers it). Each pair uses the cohort sites both samples measure (>= 20 reads and a gIVT-corrected value
in each) and that are {rule}. Pairs with fewer than 10 such sites are grey, as is the diagonal. Each panel has its
own colour scale (2nd to 98th percentile of its pairs); ICC also rises with how widely occupancy varies between
sites, so compare modifications with care. Black lines separate batches: a batch effect shows as lighter blocks.
I, ICC(2,1) over all {n} samples at the sites measured in every sample and >= 20% in at least one (the paper's
statistic for its HEK293 replicates); the shallowest sample decides which sites qualify.
J, mean pairwise ICC over pairs from the same batch (filled) and from different batches (open).
K, one point per pair, its ICC averaged over the modifications, against the smaller of the two depths. Within-batch
pairs above between-batch pairs of the same depth: batch effect; all pairs following depth together: depth.
Tables: tables/{conc['t1']} (per modification), tables/{conc['t2']} (every pair)."""
    notes.append(f'figure 5: pairwise ICC(2,1), site selection "{c.args.icc_selection}"')
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
A-H, Jaccard index for every pair: among the cohort sites both samples measure (>= 20 reads and a gIVT-corrected
value in each), the sites called (>= 20%) in both divided by the sites called in either. Each panel has its own
colour scale; black lines separate batches; grey: the sample itself. Jaccard is lower than ICC (figure 5) mostly
because of the cutoff: most sites lie at 20-30% (figure 4A), where read sampling alone moves a site across 20%
from one sample to the next (figure 12 shows the same cutoff against the binomial null across coverage).
I, mean Jaccard over pairs from the same batch (filled) and different batches (open).
J, as figure 5K, with Jaccard. Numbers: tables/{conc['t2']}."""
    return fig


# ---------------------------------------------------------------- figure 6: isoform-specific modifications
def _delta_null(m, n, reps, rng, cap=20):
    """Per row (a site in one sample, isoforms padded with n = 0): observed max - min isoform occupancy (%),
    and under the null of one shared occupancy (pooled), P(delta >= cap) and the null mean histogram."""
    valid = n > 0
    with np.errstate(invalid='ignore', divide='ignore'):
        p = np.where(valid, 100 * m / np.maximum(n, 1), np.nan)
    obs = np.nanmax(p, 1) - np.nanmin(p, 1)
    p0 = m.sum(1) / n.sum(1)
    sim = rng.binomial(n[:, :, None], p0[:, None, None], size=(n.shape[0], n.shape[1], reps))
    with np.errstate(invalid='ignore', divide='ignore'):
        ps = np.where(valid[:, :, None], 100 * sim / np.maximum(n[:, :, None], 1), np.nan)
    d = np.nanmax(ps, 1) - np.nanmin(ps, 1)
    edges = np.linspace(0, 100, 41)
    nh = np.array([np.histogram(d[:, r], bins=edges)[0] for r in range(reps)]).mean(0)
    return obs, (d >= cap).mean(1), nh


def fig6_isoforms(c, notes, min_delta=20, reps=100):
    rng = np.random.default_rng(c.args.seed)
    per_sample, deltas, null_hist, sites = {}, {}, {}, {}
    gene_name, found = {}, 0
    for s in c.samples:
        p = c.sample_file('mods', s, '.feature_mods.tsv.gz')
        if not os.path.exists(p):
            continue
        found += 1
        t = read_table(p, columns=['level', 'feature', 'gene_id', 'ref_gene_name', 'chrom', 'pos', 'strand', 'mod',
                                   'n_valid', 'percent_fpcorr'], numeric=('pos', 'n_valid', 'percent_fpcorr'))
        lv = strarr(t, 'level')
        for g_, n_ in zip(strarr(t, 'gene_id')[lv == 'gene'], strarr(t, 'ref_gene_name')[lv == 'gene']):
            if n_ not in ('.', ''):
                gene_name[g_] = n_
        v = numbers(t, ['pos', 'n_valid', 'percent_fpcorr'])
        iso = (lv == 'isoform') & (np.nan_to_num(v[:, 1]) >= 20) & np.isfinite(v[:, 2])
        keep = c.mask.keep(strarr(t, 'chrom'), v[:, 0], None) & iso
        if not keep.any():
            per_sample[s] = {}
            continue
        d = pd.DataFrame({'gene': strarr(t, 'gene_id')[keep], 'chrom': strarr(t, 'chrom')[keep],
                          'pos': v[keep, 0].astype(np.int64), 'strand': strarr(t, 'strand')[keep],
                          'mod': strarr(t, 'mod')[keep], 'feature': strarr(t, 'feature')[keep],
                          'n': v[keep, 1], 'm': np.round(v[keep, 1] * v[keep, 2] / 100)})
        d = d.sort_values('n', ascending=False)
        d['k'] = d.groupby(['gene', 'chrom', 'pos', 'strand', 'mod']).cumcount()
        d = d[d['k'] < 8]
        size = d.groupby(['gene', 'chrom', 'pos', 'strand', 'mod'])['k'].transform('size')
        d = d[size >= 2]
        if not len(d):
            per_sample[s] = {}
            continue
        keys = d[['gene', 'chrom', 'pos', 'strand', 'mod']].drop_duplicates().reset_index(drop=True)
        keys['row'] = np.arange(len(keys))
        d = d.merge(keys, on=['gene', 'chrom', 'pos', 'strand', 'mod'])
        M = np.zeros((len(keys), 8))
        Nn = np.zeros((len(keys), 8))
        M[d['row'], d['k']] = d['m']
        Nn[d['row'], d['k']] = d['n']
        mods_ = keys['mod'].to_numpy(dtype=object)
        per_sample[s] = {}
        for code in set(mods_.tolist()):
            sel = mods_ == code
            Ms, Ns = M[sel], Nn[sel].astype(np.int64)
            ol, pl_, nh = [], [], np.zeros(40)
            for s0 in range(0, len(Ms), 5000):
                o_, p_, h_ = _delta_null(Ms[s0:s0 + 5000], Ns[s0:s0 + 5000], reps, rng, min_delta)
                ol.append(o_)
                pl_.append(p_)
                nh += h_
            obs, pn = np.concatenate(ol), np.concatenate(pl_)
            per_sample[s][code] = (int((obs >= min_delta).sum()), float(pn.sum()))
            deltas.setdefault(code, []).append(obs)
            null_hist[code] = null_hist.get(code, 0) + nh
        for row in d.itertuples(index=False):
            key = (row.gene, row.chrom, row.pos, row.strand, row.mod)
            sites.setdefault(key, {}).setdefault(s, {})[row.feature] = (row.m, row.n)
    if not found:
        notes.append('figure 6 skipped: no <sample>.feature_mods.tsv.gz')
        return None
    # cohort-level test per site: generalized CMH over samples (strata), isoforms as columns
    trows, sig = [], {}
    for key, by_s in sites.items():
        feats = sorted({f for v in by_s.values() for f in v})
        strata = [(np.array([v.get(f, (0, 0))[0] for f in feats]), np.array([v.get(f, (0, 0))[1] for f in feats]))
                  for v in by_s.values()]
        q, df, pv = cmh_2xk(strata)
        mt = np.sum([m for m, _ in strata], 0)
        nt = np.sum([n for _, n in strata], 0)
        okf = nt >= 20
        if okf.sum() < 2:
            continue
        pk = 100 * mt[okf] / nt[okf]
        fk = [f for f, k in zip(feats, okf) if k]
        trows.append([key[0], gene_name.get(key[0], '.'), key[1], key[2], key[3], plain(key[4]), key[4], len(by_s),
                      float(pk.max() - pk.min()), fk[int(np.argmax(pk))], fk[int(np.argmin(pk))], q, df, pv])
    pvals = np.array([r[13] for r in trows], float)
    qv = np.full(len(trows), np.nan)
    mod_col = np.array([r[6] for r in trows], dtype=object)
    for code in set(mod_col.tolist()):
        sel = mod_col == code
        qv[sel] = bh(pvals[sel])
    for r, q_ in zip(trows, qv):
        r.append(q_)
    trows.sort(key=lambda r: (np.nan_to_num(r[14], nan=2), -r[8]))
    tname = c.table('Fig6_isoform_specific_sites.tsv',
                    ['gene_id', 'gene_name', 'chrom', 'pos', 'strand', 'mod', 'mod_code', 'samples', 'pooled_delta_pct',
                     'most_modified_isoform', 'least_modified_isoform', 'cmh_statistic', 'cmh_df', 'cmh_p', 'cmh_q'],
                    [r[:6] + r[7:] for r in trows])

    fig, gs = figure(2, 3, widths=[PANEL[0] * 1.3, PANEL[0] * 1.2, PANEL[0] * 1.2])
    ax = fig.add_subplot(gs[0, 0])                                  # A per sample: observed vs null-expected
    ss = [s for s in c.samples if s in per_sample]
    for code in [m for m in MAIN_MODS if m in deltas]:
        o = [per_sample[s].get(code, (0, 0))[0] for s in ss]
        e = [per_sample[s].get(code, (0, 0))[1] for s in ss]
        ax.plot(range(len(ss)), o, marker='o', ms=4, lw=1, color=MOD_COLOUR[code], label=MOD_NAME[code])
        ax.plot(range(len(ss)), e, marker='o', ms=4, lw=0.8, ls='--', color=MOD_COLOUR[code], mfc='white')
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylabel(f'sites with isoform Δ ≥ {min_delta}%')
    sample_axis(ax, ss, c.batch, c.colours)
    ax.legend(fontsize=7, ncol=4, loc='lower left', bbox_to_anchor=(0, 1.0), handlelength=1.5)
    ax.text(1.0, 1.12, 'solid: observed; dashed: binomial null', transform=ax.transAxes, ha='right', fontsize=7)
    letter(ax, 'A')
    ax = fig.add_subplot(gs[0, 1])                                  # B pooled delta vs null
    edges = np.linspace(0, 100, 41)
    for code in [m for m in MAIN_MODS if m in deltas]:
        h, _ = np.histogram(np.concatenate(deltas[code]), bins=edges)
        ax.step(edges[:-1], np.maximum(h, 0.8), where='post', color=MOD_COLOUR[code], lw=1.3, label=MOD_NAME[code])
        if code in null_hist:
            ax.step(edges[:-1], np.maximum(np.ravel(null_hist[code]), 0.8), where='post', color=MOD_COLOUR[code], lw=0.8, ls='--')
    ax.axvline(min_delta, color='#E00000', ls='--', lw=0.9)
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    ax.set_xlabel('max occupancy difference between isoforms (%)')
    ax.set_ylabel('site-sample pairs')
    ax.legend(fontsize=8)
    letter(ax, 'B')
    ax = fig.add_subplot(gs[0, 2])                                  # C volcano: pooled delta vs CMH p
    if trows:
        dlt = np.array([r[8] for r in trows])
        lp = -np.log10(np.clip(np.nan_to_num(pvals, nan=1), 1e-300, 1))
        isq = np.nan_to_num(qv, nan=1) < c.args.fdr
        ax.scatter(dlt[~isq], lp[~isq], s=4, color='0.75', lw=0)
        for code in ALL_MODS:
            s_ = isq & (mod_col == code)
            if s_.any():
                ax.scatter(dlt[s_], lp[s_], s=7, color=MOD_COLOUR[code], lw=0, label=MOD_NAME[code])
        stat_box(ax, [f'{int(isq.sum()):,} of {int(np.isfinite(pvals).sum()):,}', f'FDR < {c.args.fdr:g}'])
        if isq.any():
            ax.legend(fontsize=7, loc='lower right', ncol=2, handlelength=1)
    ax.set_xlabel('pooled isoform Δ (%)')
    ax.set_ylabel('-log$_{10}$ p (CMH over samples)')
    letter(ax, 'C')
    rec = {}
    for r, q_ in zip(trows, qv):
        if np.nan_to_num(q_, nan=1) < c.args.fdr:
            rec[r[0]] = rec.get(r[0], 0) + 1
    top = sorted(rec, key=lambda g: -rec[g])[:2]
    for j, g in enumerate(top):                                    # D, E example genes
        ax = fig.add_subplot(gs[1, j])
        pts = sorted([(r[3], r[8], r[6], np.nan_to_num(q_, nan=1)) for r, q_ in zip(trows, qv) if r[0] == g])
        xs_ = np.array([p_[0] for p_ in pts])
        ys_ = np.array([p_[1] for p_ in pts])
        ax.scatter(xs_, ys_, c=[MOD_COLOUR.get(p_[2], '0.5') for p_ in pts], s=28, zorder=3)
        big = np.array([p_[3] < c.args.fdr for p_ in pts])
        ax.scatter(xs_[big], ys_[big], facecolor='none', edgecolor='#E00000', s=80, lw=1, zorder=4)
        ax.set_ylim(0, 105)
        ax.set_xlabel('genomic position')
        ax.set_ylabel('pooled isoform Δ (%)')
        ax.set_title(gene_name.get(g, g), fontstyle='italic', fontsize=12)
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.tick_params(axis='x', labelsize=9, rotation=20)
        letter(ax, 'DE'[j])
    if not top:
        ax = fig.add_subplot(gs[1, 0:2])
        empty(ax, f'no site with an isoform difference at FDR < {c.args.fdr:g}')
    ax = fig.add_subplot(gs[1, 2])                                  # F significant sites by number of samples
    if trows:
        nsamp = np.array([r[7] for r in trows])
        isq = np.nan_to_num(qv, nan=1) < c.args.fdr
        ks = np.arange(1, max(nsamp) + 1)
        ax.bar(ks, [np.sum(nsamp == k) for k in ks], color='0.8', width=0.8, label='tested')
        ax.bar(ks, [np.sum((nsamp == k) & isq) for k in ks], color='#9E2A2B', width=0.8, label=f'FDR < {c.args.fdr:g}')
        ax.set_yscale('symlog', linthresh=1)
        ax.legend(fontsize=7, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2)
    ax.set_xlabel('samples measuring the site on ≥ 2 isoforms')
    ax.set_ylabel('sites')
    letter(ax, 'F')
    fig_batch_legend(fig, c.colours)
    nsig = int(np.sum(np.nan_to_num(qv, nan=1) < c.args.fdr))
    c.legends['Fig6_isoform_modifications'] = f"""Isoform-specific modification, as the paper's Fig. 5C, tested ({found} samples).
Within each sample, for every site measured on two or more isoforms of a gene (>= 20 reads of each isoform, unique
IsoQuant assignments, gIVT-corrected counts), Δ is the largest minus the smallest isoform occupancy. Δ is the range
of noisy estimates, so it is large by chance: at 20 reads per isoform and 30% occupancy, two isoforms with the same
true value differ by >= 20 points about 17% of the time, and more isoforms make it more often.
A, PER SAMPLE: sites with Δ >= {min_delta} (solid) against the number expected if every isoform at the site had the
pooled occupancy (dashed; {reps} binomial draws per site at the observed read numbers).
B, POOLED site-sample pairs by Δ (solid) and the null expectation (dashed); red line at {min_delta}.
C, one point per site: pooled Δ over samples against a generalized Cochran-Mantel-Haenszel test of modified /
unmodified by isoform, stratified by sample (samples are the strata, so a deep sample cannot carry the test
alone); coloured: BH FDR < {c.args.fdr:g} within the modification ({nsig:,} sites). Most real differences are
expected at alternative 3' ends (Nelson et al. 2026, Modulator), so read them with the isoforms' structures.
D, E, the two genes with the most such sites; circled: FDR < {c.args.fdr:g}. F, tested and significant sites by how
many samples measure them on two isoforms. Every tested site: tables/{tname}."""
    notes.append(f'figure 6: {len(trows):,} sites tested (CMH over samples), {nsig:,} at FDR < {c.args.fdr:g}')
    return fig


# ---------------------------------------------------------------- figure 7: poly(A)
def polya_data(c, notes, isoforms=()):
    """Read-level tail histograms per class and sample, plus per-read tails of the chosen isoforms."""
    out = {k: {} for k, _, _ in COMPARTMENT}
    iso_pt = {i: [] for i in isoforms}
    dropped = [0, 0]
    for s in c.samples:
        rp = c.sample_file('polya', s, '.polya_reads.tsv.gz')
        fp = c.sample_file('polya', s, '.polya_features.tsv.gz')
        if not os.path.exists(rp):
            continue
        globin_ids = set()
        if os.path.exists(fp):
            f = read_table(fp, columns=['level', 'feature_id', 'ref_gene_name'])
            globin_ids = {fid for lv, fid, nm in zip(strarr(f, 'level'), strarr(f, 'feature_id'), strarr(f, 'ref_gene_name'))
                          if lv == 'gene' and is_globin(nm)}
        cols = ['compartment', 'pt', 'status', 'gene'] + (['isoform'] if isoforms else [])
        r = read_table(rp, columns=cols, numeric=('pt',))
        used = strarr(r, 'status') == 'used'
        mito = strarr(r, 'compartment') == 'mito'
        gene = strarr(r, 'gene')
        glob = np.isin(gene, sorted(globin_ids)) & ~mito
        unassigned = np.isin(gene, ['.', ''])
        pt = numbers(r, ['pt']).ravel()
        used &= np.isfinite(pt)
        nuclear = ~mito & ~glob & ~unassigned
        dropped[0] += int(np.sum(used & ~mito & unassigned))
        dropped[1] += int(np.sum(used & ~mito))
        for k, m in (('mito', mito), ('globin', glob), ('nuclear', nuclear)):
            v = np.round(pt[used & m]).astype(int)
            if len(v):
                out[k][s] = np.unique(v, return_counts=True)
        if isoforms and 'isoform' in r.columns:
            isoc = strarr(r, 'isoform')
            hit = np.isin(isoc, list(isoforms)) & used
            for i_, p_ in zip(isoc[hit], pt[hit]):
                iso_pt[i_].append(p_)
    if dropped[1]:
        notes.append(f'figure 7: {dropped[0]:,} of {dropped[1]:,} non-mitochondrial tails ({100 * dropped[0] / dropped[1]:.1f}%) '
                     'are on reads assigned to no gene (MAPQ < 40, or outside the catalogue) and are left out')
    return out, iso_pt


def isoform_candidates(c, n=100, min_tails=10):
    """Genes whose two best-covered isoforms differ most in median tail, consistently over samples. Each:
    (|difference|, gene_id, name, isoform a, isoform b, label a, label b, median a, median b, samples a, samples b,
    tails a, tails b, {sample: (median a, median b)})."""
    per, names, labels = {}, {}, {}
    for s in c.samples:
        p = c.sample_file('polya', s, '.isoform_polya.bed')
        if not os.path.exists(p):
            continue
        t = read_table(p, numeric=('n_tails', 'median_pt'))
        if not len(t):
            continue
        v = numbers(t, ['n_tails', 'median_pt'])
        tx = strarr(t, 'ref_tx_id') if 'ref_tx_id' in t.columns else np.array(['.'] * len(t), dtype=object)
        cat = strarr(t, 'structural_category') if 'structural_category' in t.columns else np.array(['.'] * len(t), dtype=object)
        for fid, gid, nm, rt, ct, (nt, md) in zip(strarr(t, 'feature_id'), strarr(t, 'gene_id'), strarr(t, 'ref_gene_name'),
                                                  tx, cat, v):
            if nt >= min_tails and not is_globin(nm):
                per.setdefault(gid, {}).setdefault(fid, {})[s] = (nt, md)
                if nm not in ('.', ''):
                    names[gid] = nm
                base = rt if rt not in ('.', '') else fid
                labels[fid] = f'{base} ({CATEGORY_SHORT.get(ct, ct)})' if ct not in ('.', '') else base
    best = []
    for gid, isos in per.items():
        if len(isos) < 2:
            continue
        a, b = sorted(isos, key=lambda f: -sum(x[0] for x in isos[f].values()))[:2]
        ma = float(np.median([x[1] for x in isos[a].values()]))
        mb = float(np.median([x[1] for x in isos[b].values()]))
        la, lb = (labels[a], labels[b]) if labels[a] != labels[b] else (a, b)
        both = {s: (isos[a][s][1], isos[b][s][1]) for s in isos[a] if s in isos[b]}
        if len(both) >= max(3, len(c.samples) // 3):
            d_ = np.array([y - x for x, y in both.values()])
            score = abs(np.median(d_)) * math.sqrt(len(d_)) / (np.std(d_) + 5)   # consistent across samples
            best.append((abs(ma - mb), gid, names.get(gid, gid), a, b, la, lb, ma, mb, len(isos[a]), len(isos[b]),
                         int(sum(x[0] for x in isos[a].values())), int(sum(x[0] for x in isos[b].values())), both, score))
    best.sort(key=lambda x: -x[14])
    return best[:n]


def paired_test(diffs):
    d = np.asarray(diffs, float)
    d = d[np.isfinite(d) & (d != 0)]
    if len(d) >= 6:
        return 'Wilcoxon', float(stats.wilcoxon(d).pvalue)
    if len(d) >= 1:
        return 'sign test', float(stats.binomtest(int((d > 0).sum()), len(d)).pvalue)
    return 'no test', np.nan


def fig7_polya(c, notes):
    cands = isoform_candidates(c)
    isoforms = [x for b in cands[:2] for x in (b[3], b[4])]
    hists, iso_pt = polya_data(c, notes, isoforms)
    if not any(hists[k] for k in hists):
        notes.append('figure 7 skipped: no <sample>.polya_reads.tsv.gz')
        return None
    fig, gs = figure(2, 4)
    grid = np.arange(0, 301)

    def density(vals, cnt):
        d = np.zeros(len(grid))
        inside = vals <= 300
        d[np.clip(vals[inside], 0, 300)] = cnt[inside]
        return np.convolve(d / max(cnt.sum(), 1), np.ones(5) / 5, mode='same')

    ax = fig.add_subplot(gs[0, 0])                                  # A pooled reads
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
    have_gene = os.path.exists(mp) and os.path.exists(ntp)
    if have_gene:
        keys = ['feature_id', 'gene_id', 'ref_gene_name', 'chrom', 'start', 'end', 'strand']
        med, ntt = read_matrix(mp, keys), read_matrix(ntp, keys)
        ss = [s for s in c.samples if s in med.columns]
        mv, nv = numbers(med, ss), np.nan_to_num(numbers(ntt, ss))
        mv[nv < 10] = np.nan
        names = strarr(med, 'ref_gene_name')
        chroms = strarr(med, 'chrom')
        with np.errstate(invalid='ignore'):
            gmed = np.nanmedian(mv, 1)
        cls = np.array(['mito' if ch in MITO_CONTIGS else 'globin' if is_globin(nm) else 'nuclear'
                        for ch, nm in zip(chroms, names)])
        top2, meds2 = 0, []
        for k, label, col in COMPARTMENT:
            v = gmed[(cls == k) & ~np.isnan(gmed)]
            if len(v) < 3:
                continue
            h, e = np.histogram(v, bins=np.arange(0, 305, 5), density=True)
            ax.fill_between(e[:-1] + 2.5, h, color=col, alpha=0.3, lw=0)
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
    for (k, label, col), mk in zip(COMPARTMENT, ('s', '^', 'D')):
        xs = [i for i, s in enumerate(c.samples) if s in hists[k]]
        ax.plot(xs, [hist_median(*hists[k][c.samples[i]]) for i in xs], marker=mk, ms=4, lw=0.8, color=col, label=label)
    sample_axis(ax, c.samples, c.batch, c.colours)
    ax.set_ylabel('median poly(A), per read (nt)')
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, loc='best', handlelength=1.5)
    letter(ax, 'C')

    ax = fig.add_subplot(gs[0, 3])                                  # D per-sample ECDFs
    for s in c.samples:
        if s in hists['nuclear']:
            weighted_ecdf(ax, *hists['nuclear'][s], color=c.colours[c.batch[s]], lw=0.8, alpha=0.8)
        if s in hists['mito']:
            weighted_ecdf(ax, *hists['mito'][s], color=c.colours[c.batch[s]], lw=0.6, alpha=0.6, ls='--')
    ax.set_xlim(0, 300)
    ax.set_xlabel('poly(A) length (nt)')
    ax.set_ylabel('cumulative fraction of reads')
    ax.set_title('per sample: nuclear (solid), chrM (dashed)', fontsize=9)
    letter(ax, 'D')

    for j, cand in enumerate(cands[:2]):                             # E, F isoform tails: per-sample paired medians
        ax = fig.add_subplot(gs[1, j])
        _, gid, gname, a, b, la, lb = cand[:7]
        both = cand[13]
        data = [np.array(iso_pt.get(a, [])), np.array(iso_pt.get(b, []))]
        if all(len(d) >= 5 for d in data):
            parts = ax.violinplot(data, positions=[0, 1], showextrema=False, widths=0.8)
            for body in parts['bodies']:
                body.set_facecolor('0.88')
                body.set_alpha(0.7)
                body.set_edgecolor('none')
        for s, (ma, mb) in both.items():
            ax.plot([0, 1], [ma, mb], color=c.colours[c.batch[s]], lw=0.8, marker='o', ms=4, zorder=3)
        test, p = paired_test([mb - ma for ma, mb in both.values()])
        stat_box(ax, [f'{len(both)} samples', f'{test} p = {p:.2g}'], loc='upper left')
        ax.set_xticks([0, 1], [la, lb], fontsize=8, rotation=12)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel('poly(A) length (nt)')
        ax.set_ylim(bottom=0)
        ax.set_title(gname, fontsize=11, fontstyle='italic')
        letter(ax, 'EF'[j])
    if not cands:
        ax = fig.add_subplot(gs[1, 0:2])
        empty(ax, 'no gene with two isoforms of >= 10 tails in enough samples')

    ax = fig.add_subplot(gs[1, 2])                                  # G long / short tailed genes
    if have_gene:
        want = [g for g in ('DDX17', 'RPS24') if g in set(names)]
        if len(want) < 2:
            ok = (nv >= 10).all(1) & ~np.isnan(gmed) & (cls == 'nuclear')
            if ok.any():
                idx = np.where(ok)[0]
                want = [names[idx[np.argmax(gmed[idx])]], names[idx[np.argmin(gmed[idx])]]]
        nl = list(names)
        for i, gname in enumerate(want[:2]):
            y_ = mv[nl.index(gname)]
            xs = np.arange(len(ss))
            col, mk = ('0.15', 'o') if i == 0 else ('0.6', 's')
            ax.scatter(xs[~np.isnan(y_)], y_[~np.isnan(y_)], color=col, marker=mk, s=26, label=gname, zorder=3)
            ax.axhline(np.nanmedian(y_), color=col, ls='--', lw=0.8)
        sample_axis(ax, ss, c.batch, c.colours)
        ax.set_ylabel('gene median poly(A) (nt)')
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.0, 1.0), handletextpad=0.2)
    letter(ax, 'G')

    ax = fig.add_subplot(gs[1, 3])                                  # H consistency of isoform differences
    show = cands[:15]
    if not show:
        empty(ax, 'no isoform pair to show')
    for i, cand in enumerate(show):
        diffs = np.array([mb - ma for ma, mb in cand[13].values()])
        ax.scatter(diffs, np.full(len(diffs), i) + jitter(len(diffs), 0.3, i), s=8, color='0.35', lw=0)
        ax.vlines(np.median(diffs), i - 0.35, i + 0.35, color='black', lw=1.2)
    if show:
        ax.axvline(0, color='0.5', lw=0.6)
        ax.set_yticks(range(len(show)), [short_name(x[2]) for x in show], fontsize=7, fontstyle='italic')
        ax.set_ylim(len(show) - 0.5, -0.5)
        ax.set_xlabel('isoform b - isoform a, median tail (nt)')
        ax.set_title('per sample, top isoform pairs', fontsize=9)
    letter(ax, 'H')
    fig_batch_legend(fig, c.colours)
    tname = c.table('Fig7_isoform_polya_pairs.tsv', ['abs_difference', 'gene_id', 'gene_name', 'isoform_a', 'isoform_b',
                                                      'label_a', 'label_b', 'median_of_sample_medians_a',
                                                      'median_of_sample_medians_b', 'samples_a', 'samples_b', 'tails_a',
                                                      'tails_b', 'paired_test', 'paired_p', 'consistency_score'],
                    [list(x[:13]) + list(paired_test([mb - ma for ma, mb in x[13].values()])) + [x[14]] for x in cands])
    for s in c.samples:
        if s in hists['nuclear']:
            c.note_sample(s, 'polya_median_per_read_nuclear', hist_median(*hists['nuclear'][s]))
        if s in hists['mito']:
            c.note_sample(s, 'polya_median_per_read_mito', hist_median(*hists['mito'][s]))
    if have_gene:
        for j, s in enumerate(ss):
            v = mv[(cls == 'nuclear') & ~np.isnan(mv[:, j]), j]
            if len(v):
                c.note_sample(s, 'polya_median_of_gene_medians_nuclear', float(np.median(v)))
    pair = '; '.join(f'{x[2]} ({x[5]} vs {x[6]})' for x in cands[:2]) or 'none'
    c.legends['Fig7_polya'] = f"""Poly(A) tails, as the paper's Fig. S2.
A, POOLED READS: every tail from every sample on a read assigned to a gene (nuclear without globin; globin) or on
chrM, smoothed over 5 nt; dashed lines and numbers: medians. Pooling weights deep samples and abundant genes.
Reads assigned to no gene are left out (below MAPQ 40, which gene assignment requires, among them the duplicated
HBA1 and HBA2, and reads outside the catalogue).
B, ONE VALUE PER GENE: each gene's median tail in each sample with >= 10 tails, then the median over samples.
C, PER SAMPLE: the per-read median of each class. D, PER SAMPLE: the whole distribution (ECDF), nuclear solid,
chrM dashed: a sample whose tails are shifted (handling, degradation, a calibration change of the tail
estimator, which Chang et al., GigaScience 2025, found to be sample-dependent) stands out here.
E, F, the two genes whose two best-covered isoforms differ most consistently across samples (|median paired
difference| x sqrt(samples) / (SD + 5 nt), >= 3 samples): {pair}. Grey
violins: pooled reads (context only); lines: each sample's two isoform medians (isoform_polya.bed), coloured by
batch; the test is across samples (paired Wilcoxon, or a sign test below 6 samples), not across reads.
G, PER SAMPLE: gene median tail of DDX17 and RPS24 (the paper's long- and short-tailed examples; otherwise the
longest- and shortest-tailed well-covered genes); dashed: median over samples.
H, for the 15 top isoform pairs, each sample's difference (points) and the median (line).
Tables: tables/{tname}; per-sample medians: Table 1."""
    notes.append(f'figure 7: isoform pairs: {pair}')
    return fig


# ---------------------------------------------------------------- figure 8: group comparison
def composition_pc1(c):
    comp = c.composition()
    V = comp['values']
    if V.shape[1] == 0:
        return None
    ok = np.isfinite(V).all(0) & (np.nanstd(V, 0) > 0)
    if ok.sum() == 0:
        return None
    Z = (V[:, ok] - np.nanmean(V[:, ok], 0)) / np.nanstd(V[:, ok], 0)
    if not np.isfinite(Z).all():
        return None
    if Z.shape[1] == 1:
        return Z[:, 0]
    u, sv, vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    pc = u[:, 0] * sv[0]
    j = int(np.argmax(np.abs(vt[0])))                               # sign: PC1 rises with its strongest type
    return pc * np.sign(vt[0, j])


def design_matrix(c, samples, covariate_names):
    """Intercept, group (B vs A), then covariates (numeric standardised, categorical one-hot). Columns that do
    not add rank (e.g. batch nested in group) are dropped and reported."""
    a, b = c.levels[:2]
    covs = c.covariates()
    pc1 = composition_pc1(c)
    if pc1 is not None:
        covs['composition PC1'] = ('num', pc1)
    pos = {s: i for i, s in enumerate(c.samples)}
    cols, names, notes = [np.ones(len(samples)), np.array([1.0 if c.group[s] == b else 0.0 for s in samples])], ['intercept', f'{b} vs {a}'], []
    for nm in covariate_names:
        if nm not in covs:
            notes.append(f'{nm}: not a covariate here')
            continue
        kind, vals = covs[nm]
        v = [vals[pos[s]] for s in samples]
        if kind == 'num':
            v = np.array(v, float)
            if not np.isfinite(v).all() or np.std(v) == 0:
                notes.append(f'{nm}: missing or constant, left out')
                continue
            cols.append((v - v.mean()) / v.std())
            names.append(nm)
        else:
            lv = sorted(set(v))
            for l_ in lv[1:]:
                cols.append(np.array([1.0 if x == l_ else 0.0 for x in v]))
                names.append(f'{nm}={l_}')
    X = np.column_stack(cols[:2])
    keep_names = names[:2]
    for col, nm in zip(cols[2:], names[2:]):
        trial = np.column_stack([X, col])
        if np.linalg.matrix_rank(trial) > np.linalg.matrix_rank(X) and trial.shape[1] < len(samples) - 1:
            X, keep_names = trial, keep_names + [nm]
        else:
            notes.append(f'{nm}: collinear with the group or too many terms for {len(samples)} samples, left out')
    return X, keep_names, notes


def default_model_covariates(c):
    out = []
    covs = c.covariates()
    if c.args.batch_column in covs:
        out.append(c.args.batch_column)
    if composition_pc1(c) is not None:
        out.append('composition PC1')
    return out


def volcano(ax, effect, p, q, labels, xlabel, colour, q_max, top=8, ylabel='-log$_{10}$ p'):
    lp = -np.log10(np.clip(np.nan_to_num(p, nan=1.0), 1e-300, 1))
    sig = np.nan_to_num(q, nan=1.0) < q_max
    ax.scatter(effect[~sig], lp[~sig], s=5, color='0.75', lw=0)
    ax.scatter(effect[sig], lp[sig], s=7, color=colour, lw=0)
    for i in np.argsort(np.nan_to_num(p, nan=1.0))[:top]:
        ax.annotate(labels[i], (effect[i], lp[i]), fontsize=7, xytext=(3, 2), textcoords='offset points')
    ax.axvline(0, color='0.5', lw=0.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    stat_box(ax, [f'{int(sig.sum()):,} of {int(np.sum(~np.isnan(p))):,}', f'FDR < {q_max:g}'], loc='upper left')


def isoform_rows_at(c, s, keyset):
    """{(chrom, pos, strand, mod): {isoform: (corrected m, n)}} for the given sites, from feature_mods."""
    p = c.sample_file('mods', s, '.feature_mods.tsv.gz')
    out = {}
    if not os.path.exists(p) or not keyset:
        return out
    t = read_table(p, columns=['level', 'feature', 'chrom', 'pos', 'strand', 'mod', 'n_valid', 'percent_fpcorr'],
                   numeric=('pos', 'n_valid', 'percent_fpcorr'))
    lv = strarr(t, 'level')
    v = numbers(t, ['pos', 'n_valid', 'percent_fpcorr'])
    ch, sd, md, ft = strarr(t, 'chrom'), strarr(t, 'strand'), strarr(t, 'mod'), strarr(t, 'feature')
    for i in np.where(lv == 'isoform')[0]:
        k = (ch[i], int(v[i, 0]), sd[i], md[i])
        if k in keyset and np.isfinite(v[i, 2]) and v[i, 1] > 0:
            out.setdefault(k, {})[ft[i]] = (v[i, 1] * v[i, 2] / 100, v[i, 1])
    return out


def shift_share(c, sites, samples, a, b, min_n=10):
    """Site-level occupancy change between groups split into within-isoform change and isoform-usage change:
    sum_k w_k pi_k (w = share of the site's reads on isoform k, pi = isoform occupancy, both averaged over a
    group's samples); B - A = sum w_bar dpi + sum pi_bar dw, exactly."""
    keyset = set(sites)
    per = {s: isoform_rows_at(c, s, keyset) for s in samples}
    out = []
    for k in sites:
        grp = {}
        for g_ in (a, b):
            ws, ps = {}, {}
            for s in samples:
                if c.group[s] != g_ or k not in per[s]:
                    continue
                iso = per[s][k]
                tot = sum(n for _, n in iso.values())
                for f, (m, n) in iso.items():
                    ws.setdefault(f, []).append(n / tot)
                    if n >= min_n:
                        ps.setdefault(f, []).append(m / n)
            grp[g_] = ({f: np.mean(v) for f, v in ws.items()}, {f: np.mean(v) for f, v in ps.items()})
        feats = [f for f in grp[a][1] if f in grp[b][1]]
        if len(feats) < 1:
            continue
        wa = np.array([grp[a][0].get(f, 0) for f in feats])
        wb = np.array([grp[b][0].get(f, 0) for f in feats])
        if wa.sum() == 0 or wb.sum() == 0:
            continue
        wa, wb = wa / wa.sum(), wb / wb.sum()
        pa = np.array([grp[a][1][f] for f in feats])
        pb = np.array([grp[b][1][f] for f in feats])
        within = float(np.sum((wa + wb) / 2 * (pb - pa)))
        usage = float(np.sum((pa + pb) / 2 * (wb - wa)))
        out.append((k, 100 * within, 100 * usage, len(feats)))
    return out


def fig8_groups(c, fdr, notes):
    if not c.group or len(c.levels or []) < 2:
        sys.exit(f'cohort_plots.py: --group {c.args.group} needs two levels; pick them with --levels A B')
    a, b = c.levels[:2]
    if len(c.levels) > 2:
        notes.append(f'figure 8: --group has {len(c.levels)} levels; comparing {b} against {a}')
    label, tag = f'{b} vs {a}', f'{c.args.group}.{b}_vs_{a}'
    covnames = c.args.model_covariates if c.args.model_covariates is not None else default_model_covariates(c)
    fig, gs = figure(3, 4, heights=[PANEL[1] * 1.25, PANEL[1] * 1.25, PANEL[1] * 1.15])
    S = [s for s in c.samples if c.group.get(s) in (a, b)]
    ax = fig.add_subplot(gs[0, 0])                                  # A design balance
    batches = sorted(set(c.batch[s] for s in S))
    tab = np.array([[sum(1 for s in S if c.group[s] == g_ and c.batch[s] == bt) for bt in batches] for g_ in (a, b)])
    ax.pcolormesh(np.arange(len(batches) + 1) - 0.5, np.arange(3) - 0.5, tab, cmap=HEAT_CMAP, vmin=0,
                  vmax=max(1, tab.max()), edgecolors='white', linewidth=1)
    for i in range(2):
        for j in range(len(batches)):
            ax.text(j, i, str(tab[i, j]), ha='center', va='center', fontsize=10,
                    color='white' if tab[i, j] > tab.max() / 2 else 'black')
    ax.set_xticks(range(len(batches)), batches, rotation=30, fontsize=9)
    for lab, bt in zip(ax.get_xticklabels(), batches):
        lab.set_color(c.colours[bt])
    ax.set_yticks([0, 1], [a, b])
    ax.set_ylim(1.5, -0.5)
    ax.set_xlabel('batch')
    ax.tick_params(length=0)
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(False)
    confounded = all((tab[:, j] > 0).sum() == 1 for j in range(len(batches)))
    ax.set_title('samples per group and batch' + ('  (group nested in batch!)' if confounded else ''), fontsize=9,
                 color='#E00000' if confounded else 'black')
    letter(ax, 'A')

    ax = fig.add_subplot(gs[0, 1])                                  # B composition by group
    comp = c.composition()
    if comp['names']:
        pos = {s: i for i, s in enumerate(c.samples)}
        show = comp['names'][:6]
        for j, nm in enumerate(show):
            vals = comp['values'][:, comp['names'].index(nm)]
            for k, g_ in enumerate((a, b)):
                v = [vals[pos[s]] for s in S if c.group[s] == g_]
                strip_points(ax, j + (k - 0.5) * 0.35, v, ['#5769AB', '#E76969'][k], width=0.3, s=12, floor=None)
        ax.set_xticks(range(len(show)), show, rotation=40, fontsize=8)
        ax.set_ylabel(comp['kind'])
        ax.legend([plt.Line2D([], [], marker='o', ls='', color=col) for col in ('#5769AB', '#E76969')], [a, b], fontsize=8)
    else:
        empty(ax, 'no composition: give --composition or gene CPM')
    letter(ax, 'B')

    X, xnames, dnotes = design_matrix(c, S, covnames)
    notes.extend(f'figure 8 model: {n_}' for n_ in dnotes)
    model_txt = ' + '.join(xnames[1:])
    es, lx, names, gids = c.expression()
    ei = [es.index(s) for s in S if s in es]
    Xe, _, _ = design_matrix(c, [es[i] for i in ei], covnames)
    est, t, p, dfr, d0 = moderated_lm(lx[:, ei], Xe, 1)
    q = bh(p)
    ax = fig.add_subplot(gs[0, 2])                                  # C expression
    volcano(ax, est, p, q, names, f'log$_2$ fold change, {label}', '#E76969', fdr, ylabel='-log$_{10}$ p (moderated t)')
    letter(ax, 'C')
    c.table(f'expression.{tag}.tsv', ['gene', 'log2fc', 't', 'p', 'q'],
            sorted(zip(names, est, t, p, q), key=lambda r: (math.isnan(r[3]), r[3])))
    ax = fig.add_subplot(gs[0, 3])                                  # D top genes, z-scores
    top = [i for i in np.argsort(np.nan_to_num(p, nan=1.0))[:30]]
    order = [i for i in ei if c.group[es[i]] == a] + [i for i in ei if c.group[es[i]] == b]
    z = lx[np.ix_(top, order)]
    z = (z - z.mean(1, keepdims=True)) / (z.std(1, keepdims=True) + 1e-9)
    im = ax.pcolormesh(np.arange(len(order) + 1) - 0.5, np.arange(len(top) + 1) - 0.5, z, cmap='RdBu_r', vmin=-2.5,
                       vmax=2.5, edgecolors='face', linewidth=0.3)
    ax.set_ylim(len(top) - 0.5, -0.5)
    sample_axis(ax, [es[i] for i in order], c.batch, c.colours, fontsize=7)
    ax.set_yticks(range(len(top)), [names[i] for i in top], fontsize=6)
    na_ = sum(1 for i in order if c.group[es[i]] == a)
    ax.axvline(na_ - 0.5, color='black', lw=1)
    ax.tick_params(length=0)
    colourbar(fig, im, ax, label='z-score')
    letter(ax, 'D')

    rows, ms, pct, mods, nv = c.mods_matrix()
    nvc, nmc, fpc, _ = c.mods_counts()
    mi = [ms.index(s) for s in S if s in ms]
    msel = [ms[i] for i in mi]
    ga = np.array([c.group[s] == a for s in msel])
    gb = ~ga
    ax = fig.add_subplot(gs[1, 0])                                  # E modifications
    res, top_sites = None, []
    if nv is not None and len(mi) >= 4:
        P, Nv = pct[:, mi], nv[:, mi]
        meas = np.isfinite(P) & (np.nan_to_num(Nv) >= 20)
        enough = (meas[:, ga].sum(1) >= min(3, ga.sum())) & (meas[:, gb].sum(1) >= min(3, gb.sum()))
        idx = np.where(enough)[0]
        if len(idx) >= 10:
            Xm, _, _ = design_matrix(c, msel, covnames)
            mcorr = np.where(meas[idx], P[idx] * Nv[idx] / 100, np.nan)
            res = re_logit(mcorr, np.where(meas[idx], Nv[idx], np.nan), Xm, coef=1)
            with np.errstate(invalid='ignore'):
                dpct = np.nanmean(np.where(meas[idx][:, gb], P[idx][:, gb], np.nan), 1) - \
                    np.nanmean(np.where(meas[idx][:, ga], P[idx][:, ga], np.nan), 1)
            pv = res['p']
            qv = np.full(len(idx), np.nan)
            for code in set(mods[idx].tolist()):
                s_ = mods[idx] == code
                qv[s_] = bh(pv[s_])
            lab = [f'{ch}:{int(st) + 1} {plain(mo)}' for ch, st, mo in zip(rows['chrom'][idx], rows['start'][idx], mods[idx])]
            sig = np.nan_to_num(qv, nan=1.0) < fdr
            lp = -np.log10(np.clip(np.nan_to_num(pv, nan=1.0), 1e-300, 1))
            ax.scatter(dpct[~sig], lp[~sig], s=5, color='0.75', lw=0)
            for code in ALL_MODS:
                s_ = sig & (mods[idx] == code)
                if s_.any():
                    ax.scatter(dpct[s_], lp[s_], s=8, color=MOD_COLOUR[code], lw=0, label=MOD_NAME[code])
            for i in np.argsort(np.nan_to_num(pv, nan=1.0))[:6]:
                ax.annotate(lab[i], (dpct[i], lp[i]), fontsize=6, xytext=(3, 2), textcoords='offset points')
            ax.axvline(0, color='0.5', lw=0.6)
            stat_box(ax, [f'{int(sig.sum()):,} of {int(np.isfinite(pv).sum()):,}', f'FDR < {fdr:g}'])
            if sig.any():
                ax.legend(fontsize=7, loc='upper right', ncol=2, handlelength=1)
            ax.set_xlabel(f'Δ occupancy (%), {label}')
            ax.set_ylabel('-log$_{10}$ p (random-effects logit)')
            c.table(f'modifications.{tag}.tsv', ['site', 'chrom', 'start', 'strand', 'mod', 'delta_pct', 'logit_effect',
                                                 'se', 'tau2', 'samples', 'p', 'q'],
                    sorted(zip(lab, rows['chrom'][idx], rows['start'][idx], rows['strand'][idx], [plain(m) for m in mods[idx]],
                               dpct, res['beta'], res['se'], res['tau2'], res['k'], pv, qv),
                           key=lambda r: (math.isnan(r[10]), r[10])))
            top_sites = [idx[i] for i in np.argsort(np.nan_to_num(pv, nan=1.0))[:6]]
            res['idx'], res['q'], res['dpct'] = idx, qv, dpct
        else:
            empty(ax, 'fewer than 10 sites measured in\n≥ 3 samples of each group')
    else:
        empty(ax, 'needs the n_valid matrix and ≥ 4 samples')
    letter(ax, 'E')

    ax = fig.add_subplot(gs[1, 1])                                  # F forest plot of the top sites
    if top_sites and nvc is not None:
        for r_, si in enumerate(top_sites):
            for k, s in enumerate(msel):
                j = ms.index(s)
                n_, m_ = nvc[si, j], nmc[si, j]
                if not (np.isfinite(n_) and n_ > 0 and np.isfinite(pct[si, j])):
                    continue
                lo, hi = wilson(m_, n_)
                f_ = 0 if not np.isfinite(fpc[si, j]) else fpc[si, j]
                lo, hi = max(0, 100 * (lo - f_)), max(0, 100 * (hi - f_))
                y = r_ + (0.18 if c.group[s] == b else -0.18) + jitter(1, 0.16, k * 31 + r_)[0]
                ax.plot([lo, hi], [y, y], color=c.colours[c.batch[s]], lw=0.7)
                ax.scatter([pct[si, j]], [y], color=c.colours[c.batch[s]], marker=c.gmarker[c.group[s]], s=14, zorder=3)
        ax.set_yticks(range(len(top_sites)), [f'{rows["chrom"][si]}:{rows["start"][si] + 1}\n{plain(mods[si])}' for si in top_sites], fontsize=7)
        ax.set_ylim(len(top_sites) - 0.5, -0.5)
        ax.set_xlim(0, 100)
        ax.set_xlabel('corrected occupancy (%)')
        ax.set_title(f'per sample, 95% Wilson interval; {a} above, {b} below', fontsize=8)
    else:
        empty(ax, 'no site to show')
    letter(ax, 'F')

    ax = fig.add_subplot(gs[1, 2])                                  # G poly(A)
    mp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.median_pt.tsv.gz')
    ntp = os.path.join(c.cohort, f'cohort.{c.ref}.polya.gene.n_tails.tsv.gz')
    if os.path.exists(mp) and os.path.exists(ntp):
        keys = ['feature_id', 'gene_id', 'ref_gene_name', 'chrom', 'start', 'end', 'strand']
        med, ntt = read_matrix(mp, keys), read_matrix(ntp, keys)
        ss = [s for s in S if s in med.columns]
        pv_ = numbers(med, ss)
        pv_[np.nan_to_num(numbers(ntt, ss)) < 10] = np.nan
        gg = np.array([c.group[s] for s in ss])
        ok = ((np.isfinite(pv_[:, gg == a]).sum(1) >= 3) & (np.isfinite(pv_[:, gg == b]).sum(1) >= 3))
        if ok.sum() >= 10:
            Xp, _, _ = design_matrix(c, ss, covnames)
            est2, t2, p2, _, _ = moderated_lm(pv_[ok], Xp, 1)
            q2 = bh(p2)
            gn = [n if n not in ('.', '') else g_ for g_, n in zip(strarr(med, 'gene_id')[ok], strarr(med, 'ref_gene_name')[ok])]
            volcano(ax, est2, p2, q2, gn, f'Δ median poly(A) (nt), {label}', '#5769AB', fdr, ylabel='-log$_{10}$ p (moderated t)')
            c.table(f'polya.{tag}.tsv', ['gene', 'delta_median_pt', 't', 'p', 'q'],
                    sorted(zip(gn, est2, t2, p2, q2), key=lambda r: (math.isnan(r[3]), r[3])))
        else:
            empty(ax, 'fewer than 10 genes with ≥ 10 tails\nin ≥ 3 samples of each group')
    else:
        empty(ax, 'no gene poly(A) matrix')
    letter(ax, 'G')

    ax = fig.add_subplot(gs[1, 3])                                  # H power
    na, nb = int(ga.sum()), int(gb.sum())
    cov = np.logspace(np.log10(20), 3, 50)
    p0, alpha, power = 0.3, 1e-3, 0.8
    zz = stats.norm.isf(alpha / 2) + stats.norm.isf(1 - power)
    for code in MAIN_MODS:
        if res is None:
            break
        s_ = mods[res['idx']] == code
        if s_.sum() < 10:
            continue
        t2m = float(np.nanmedian(res['tau2'][s_]))
        for mult, ls in ((1, '-'), (2, ':')):
            se = np.sqrt((t2m + 1 / (cov * p0 * (1 - p0))) * (1 / (na * mult) + 1 / (nb * mult)))
            d = 100 * (expit(logit(p0) + zz * se) - p0)
            ax.plot(cov, d, color=MOD_COLOUR[code], ls=ls, lw=1.2, label=f'{MOD_NAME[code]} (τ² {t2m:.2f})' if mult == 1 else None)
    ax.set_xscale('log')
    log_axis(ax.xaxis)
    ax.set_xlabel('reads per sample at the site')
    ax.set_ylabel(f'detectable Δ (points, from {int(100 * p0)}%)')
    ax.set_title(f'{na} vs {nb} samples (solid), twice that (dotted)', fontsize=8)
    ax.legend(fontsize=7)
    letter(ax, 'H')

    ax = fig.add_subplot(gs[2, 0:2])                                # I site change: within isoform vs usage
    dec = []
    if res is not None:
        order_ = np.argsort(np.nan_to_num(res['p'], nan=1.0))[:40]
        keys_ = [(rows['chrom'][res['idx'][i]], int(rows['start'][res['idx'][i]]), rows['strand'][res['idx'][i]],
                  mods[res['idx'][i]]) for i in order_]
        dec = shift_share(c, keys_, msel, a, b)
    if dec:
        dec = dec[:20]
        y = np.arange(len(dec))
        within = np.array([d_[1] for d_ in dec])
        usage = np.array([d_[2] for d_ in dec])
        ax.barh(y, within, color='#3B0000', height=0.7, label='within isoforms')
        ax.barh(y, usage, left=np.where(np.sign(usage) == np.sign(within), within, 0), color='#C9A9A2', height=0.7,
                label='isoform usage')
        ax.axvline(0, color='0.4', lw=0.6)
        ax.set_yticks(y, [f'{k[0]}:{k[1] + 1} {plain(k[3])} ({n_} iso)' for k, _, _, n_ in dec], fontsize=6)
        ax.set_ylim(len(dec) - 0.5, -0.5)
        ax.set_xlabel(f'Δ occupancy (points), {label}')
        ax.legend(fontsize=8, loc='lower right')
        c.table(f'modifications_shift_share.{tag}.tsv', ['chrom', 'start', 'strand', 'mod', 'within_isoform_pts', 'usage_pts', 'isoforms'],
                [(k[0], k[1], k[2], plain(k[3]), w_, u_, n_) for k, w_, u_, n_ in dec])
    else:
        empty(ax, 'no top site with isoform-level counts in both groups')
    letter(ax, 'I')

    ax = fig.add_subplot(gs[2, 2:])                                 # J top site against composition
    pc1 = composition_pc1(c)
    if top_sites and pc1 is not None:
        pos = {s: i for i, s in enumerate(c.samples)}
        si = top_sites[0]
        for s in msel:
            j = ms.index(s)
            if np.isfinite(pct[si, j]):
                ax.scatter(pc1[pos[s]], pct[si, j], color=c.colours[c.batch[s]], marker=c.gmarker[c.group[s]], s=30, zorder=3)
        ok = np.array([np.isfinite(pct[si, ms.index(s)]) for s in msel])
        xs = np.array([pc1[pos[s]] for s in msel])[ok]
        ys = np.array([pct[si, ms.index(s)] for s in msel])[ok]
        if len(xs) >= 4:
            r_, p_ = stats.spearmanr(xs, ys)
            stat_box(ax, [f'Spearman ρ = {r_:.2f}', f'p = {p_:.2g}'], loc='upper right')
        ax.set_xlabel('composition PC1 (' + (c.composition()['kind']) + ')')
        ax.set_ylabel('corrected occupancy (%)')
        ax.set_title(f'top site {rows["chrom"][si]}:{rows["start"][si] + 1} {plain(mods[si])}: does it follow composition?', fontsize=9)
    else:
        empty(ax, 'needs a significant site and composition')
    letter(ax, 'J')
    extra = [plt.Line2D([], [], marker=c.gmarker[lv], ls='', color='0.3', label=lv) for lv in (a, b)]
    fig_batch_legend(fig, c.colours, extra)
    notes.append(f'figure 8: {label}, {len(S)} samples; model ~ {model_txt}')
    c.legends['Fig8_group_comparison'] = f"""{c.args.group}: {b} against {a}. Samples are the units throughout; model: y ~ {model_txt}.
A, samples per group and batch: if the group is nested in batch, no model can separate them (flagged in red).
B, composition by group ({c.composition()['kind']}, {c.composition()['source']}): whole-blood expression, and the
occupancy of sites in cell-type-specific transcripts, follow leukocyte composition, so composition is in the model.
C, expression: per gene, linear model on log2(CPM + 1) (CPM without globin and chrM genes) with limma-style
empirical-Bayes moderated t (Smyth 2004); BH FDR; coloured FDR < {fdr:g}. Use DESeq2 / edgeR / dream on the counts
before claiming a gene. D, z-scores of the 30 genes with the lowest p.
E, modification occupancy: per site, logit of the gIVT-corrected fraction with its binomial sampling variance
treated as known measurement error, plus a between-sample variance tau^2 (REML on a grid, empirical-Bayes prior
over the grid shared by sites); Wald z on the group coefficient (size 4-5% at alpha 0.05 in simulation with 6-20
samples). Sites measured in >= 3
samples of each group. x: difference of group means of corrected occupancy. For the reference implementation of
the same idea use DSS (DMLfit.multiFactor, smoothing off) on n_mod / n_valid. F, the six lowest-p sites: each
sample's corrected occupancy and 95% Wilson interval (raw counts minus the gIVT rate), {a} above and {b} below.
G, gene median poly(A) (>= 10 tails), the same moderated linear model. H, smallest change in occupancy (from
30%) detectable at alpha = 1e-3 with 80% power for these group sizes (solid) and twice them (dotted), using the
median tau^2 of each modification: below the curve, a null result says little.
I, for the 20 lowest-p sites with isoform-level counts: the group difference split exactly into change within
isoforms (sum of mean isoform shares x change in isoform occupancy) and change in isoform usage (sum of mean
isoform occupancy x change in shares). A site whose change is usage is an isoform result, not a modification
result (cf. Nelson et al. 2026, Modulator). J, the top site's occupancy against composition PC1.
Tables: tables/expression.{tag}.tsv, modifications.{tag}.tsv, polya.{tag}.tsv, modifications_shift_share.{tag}.tsv."""
    c.fig8_res = res
    return fig


# ---------------------------------------------------------------- figure 9: isoform usage between samples
def fig9_usage(c, fdr, notes, min_reads=20):
    base = c.path('isoforms', c.ref, 'high_confidence')
    cp = os.path.join(base, f'cohort.{c.ref}.counts.tsv')
    clp = os.path.join(base, f'cohort.{c.ref}.high_confidence.classification.tsv')
    if not (os.path.exists(cp) and os.path.exists(clp)):
        notes.append(f'figure 9 skipped: needs {cp} and {clp}')
        return None
    cl = read_table(clp)
    col = (lambda k: strarr(cl, k) if k in cl.columns else np.array(['.'] * len(cl), dtype=object))
    info = {i: (rt, cat) for i, rt, cat in zip(strarr(cl, 'isoform_id'), col('ref_tx_id'), col('structural_category'))}
    cnt = read_matrix(cp, ['transcript_id', 'gene_id'])
    ss = [s for s in c.samples if s in cnt.columns]
    n = len(ss)
    keep = np.array([t in info for t in strarr(cnt, 'transcript_id')])
    X = np.nan_to_num(numbers(cnt, ss))[keep]
    tx = list(strarr(cnt, 'transcript_id')[keep])
    gid = list(strarr(cnt, 'gene_id')[keep])
    names = c.gene_names()
    label = {t: (info[t][0] if info[t][0] not in ('.', '') else t) for t in tx}
    by_gene = {}
    for i, g in enumerate(gid):
        by_gene.setdefault(g, []).append(i)
    min_samples = max(3, math.ceil(0.75 * n))
    grp = None
    if c.group and len(c.levels or []) >= 2:
        grp = (c.levels[:2], np.array([c.group.get(s) for s in ss], dtype=object))
    F_all = np.full(X.shape, np.nan)
    rows, sw, seen, examples, fm_by_gene = [], np.zeros(n), np.zeros(n), [], {}
    globin_skipped = 0
    for g, idx in by_gene.items():
        if len(idx) < 2:
            continue
        if is_globin(names.get(g, '')):                             # near-identical copies: mapping, not usage
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
            continue
        F_all[idx] = F
        fm = F[major, ok]
        fm_by_gene[g] = F[major]
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
               label[tx[idx[major]]], float(pooled[major]), float(fm.min()), float(fm.max()), obs_sd, exp_sd, obs_sd / exp_sd, p,
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
              'major_isoform', 'major_label', 'major_fraction_pooled', 'major_fraction_min', 'major_fraction_max',
              'sd_observed', 'sd_expected_binomial', 'sd_ratio', 'chi2_p', 'chi2_q', 'samples_major_not_dominant']
    if grp:
        (la, lb), _ = grp
        qg = bh([r[-1] for r in rows])
        for r, q_ in zip(rows, qg):
            r.append(q_)
        header += [f'delta_major_fraction_{lb}_vs_{la}', 'welch_p', 'welch_q']
    # composition: how strongly does each gene's major-isoform fraction track composition PC1?
    pc1 = composition_pc1(c)
    pos = {s: i for i, s in enumerate(c.samples)}
    rcomp = {}
    if pc1 is not None:
        xs_ = np.array([pc1[pos[s]] for s in ss])
        for r in rows:
            f_ = fm_by_gene[r[0]]
            ok = np.isfinite(f_)
            if ok.sum() >= 5:
                rcomp[r[0]] = stats.spearmanr(xs_[ok], f_[ok])[0]
    for r in rows:
        r.append(rcomp.get(r[0], np.nan))
    header.append('spearman_major_fraction_vs_composition_pc1')
    rows.sort(key=lambda r: -r[13])
    t1 = c.table('Fig9_isoform_usage_genes.tsv', header, rows)
    used = {r[0] for r in rows}
    frows = [[t, g, names.get(g, g), info[t][0], info[t][1]] + [F_all[i, j] for j in range(n)]
             for i, (t, g) in enumerate(zip(tx, gid)) if g in used]
    t2 = c.table('Fig9_isoform_fractions.tsv', ['transcript_id', 'gene_id', 'gene_name', 'ref_tx_id', 'structural_category'] + ss, frows)
    cats = {}
    for t in info:
        cats[info[t][1]] = cats.get(info[t][1], 0) + 1
    t3 = c.table('Fig9_isoform_categories.tsv', ['structural_category', 'high_confidence_isoforms'],
                 sorted(cats.items(), key=lambda kv: -kv[1]))
    for j, s in enumerate(ss):
        if seen[j]:
            c.note_sample(s, 'pct_genes_major_isoform_not_dominant', round(100 * sw[j] / seen[j], 1))

    fig, gs = figure(2, 4, widths=[PANEL[0] * 1.1, PANEL[0] * 1.1, PANEL[0] * 1.2, PANEL[0] * 1.1])
    ax = fig.add_subplot(gs[0, 0])                                  # A observed vs expected variation
    obs, exp = np.array([r[11] for r in rows]), np.array([r[12] for r in rows])
    ratio, rng_ = obs / exp, np.array([r[10] - r[9] for r in rows])
    hit = (np.nan_to_num(np.array([r[15] for r in rows], float), nan=1) < fdr) & (ratio >= 2) & (rng_ >= 0.2)
    lo_ = max(1e-3, min(exp.min(), obs[obs > 0].min() if (obs > 0).any() else 1e-3))
    hi_ = max(exp.max(), obs.max())
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
        topv = np.argsort(Z.var(1))[::-1][:min(500, len(Z))]
        pcs, fr = pca(Z[topv], 2)
        for i, s in enumerate(ss):
            mk = c.gmarker.get(c.group.get(s), 'o') if c.group else 'o'
            ax.scatter(pcs[i, 0], pcs[i, 1], color=c.colours[c.batch[s]], s=40, zorder=3, marker=mk)
            ax.annotate(s[-2:], (pcs[i, 0], pcs[i, 1]), fontsize=8, xytext=(3, 2), textcoords='offset points')
        ax.set_xlabel(f'PC1 ({100 * fr[0]:.1f}%)')
        ax.set_ylabel(f'PC2 ({100 * fr[1]:.1f}%)')
        ax.set_title(f'isoform fractions, {len(topv)} most variable isoforms', fontsize=10)
        c.isoform_pcs = (ss, pcs)
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
    ax = fig.add_subplot(gs[0, 3])                                  # G usage variation against composition
    if rcomp:
        rr = np.array([abs(rcomp.get(r[0], np.nan)) for r in rows])
        okr = np.isfinite(rr)
        ax.scatter(ratio[okr & ~hit], rr[okr & ~hit], s=4, color='0.7', lw=0)
        ax.scatter(ratio[okr & hit], rr[okr & hit], s=7, color='#9E2A2B', lw=0)
        ax.set_xscale('log')
        log_axis(ax.xaxis)
        ax.set_xlabel('SD ratio (observed / expected)')
        ax.set_ylabel('|Spearman ρ| with composition PC1')
        if (okr & hit).any():
            stat_box(ax, [f'median |ρ|: flagged {np.median(rr[okr & hit]):.2f}', f'others {np.median(rr[okr & ~hit]):.2f}'], loc='upper left', fontsize=8)
    else:
        empty(ax, 'no composition')
    letter(ax, 'G')
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
        ax = fig.add_subplot(gs[1, 0:3])
        empty(ax, 'no gene with >= 100 reads in every sample and a major-isoform range >= 0.2')
    fig_batch_legend(fig, c.colours)
    shown = '; '.join(short_name(names.get(e[1], e[1])) for e in examples[:3]) or 'none'
    c.legends['Fig9_isoform_usage'] = f"""Isoform usage compared BETWEEN samples ({n} samples, {len(rows):,} genes; {globin_skipped} globin genes left out,
since reads split between near-identical globin copies by mapping, not by isoform use).
For every gene with >= 2 high-confidence isoforms, the fraction of the gene's reads on each isoform, per sample
(IsoQuant requantification counts). A sample is used for a gene where the gene has >= {min_reads} reads; genes need
>= {min_samples} such samples. The major isoform has the largest pooled share.
A, ONE POINT PER GENE: SD of the major isoform's fraction across samples against the SD read sampling alone would
give (binomial). Dashed: equal; dotted: twice. Dark red: >= 2x, range >= 0.2 and chi-square FDR < {fdr:g} (the
chi-square assumes binomial counts and finds almost everything: rank by the SD ratio).
B, PCA of samples on isoform fractions (500 most variable isoforms of genes measured in every sample).
C, PER SAMPLE: % of measured genes whose major isoform is not that sample's most used isoform.
G, the SD ratio against how strongly each gene's major fraction follows composition PC1 ({c.composition()['kind']}):
in whole blood many isoform-usage differences are cell-type composition (neutrophil against lymphocyte isoforms).
D-F, the three genes with >= 100 reads in every sample, a major-fraction range >= 0.2 and the highest SD ratio
({shown}). Exploratory; for tests between groups use DRIMSeq or satuRn on
isoforms/{c.ref}/high_confidence/cohort.{c.ref}.counts.tsv. Tables: {t1}, {t2}, {t3}."""
    notes.append(f'figure 9: {len(rows):,} genes ({globin_skipped} globin genes left out); {int(hit.sum()):,} with SD >= 2x; examples {shown}')
    return fig


# ---------------------------------------------------------------- figure S1: one batch against the rest
def figS1_batch(c, compare, notes):
    _, ms, pct, mods, _ = c.mods_matrix()
    names = sorted(set(c.batch.values()))
    if len(names) < 2:
        notes.append('figure S1 skipped: one batch only')
        return None
    a = compare or names[-1]
    if a not in names:
        sys.exit(f'cohort_plots.py: --compare-batch {a} is not a batch ({", ".join(names)})')
    fig, gs = figure(3, 4)
    ia = np.array([i for i, s in enumerate(ms) if c.batch[s] == a])
    ib = np.array([i for i, s in enumerate(ms) if c.batch[s] != a])
    with np.errstate(invalid='ignore'):
        meas_a, meas_b = (~np.isnan(pct[:, ia])).sum(1), (~np.isnan(pct[:, ib])).sum(1)
        mean_a, mean_b = np.nanmean(pct[:, ia], 1), np.nanmean(pct[:, ib], 1)
    keep = (meas_a >= min(2, len(ia))) & (meas_b >= min(2, len(ib)))
    vmax = 10
    for code in ALL_MODS:
        sel = keep & (mods == code) & ~np.isnan(mean_a) & ~np.isnan(mean_b)
        if sel.any():
            vmax = max(vmax, np.histogram2d(mean_b[sel], mean_a[sel], bins=40, range=[(0, 100), (0, 100)])[0].max())
    last, axes, rvals = None, [], {}
    for k, code in enumerate(ALL_MODS):
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
    ax = fig.add_subplot(gs[2, :2])
    codes = [m for m in ALL_MODS if not np.isnan(rvals.get(m, (np.nan, 0))[0])]
    ax.bar(range(len(codes)), [rvals[m][0] for m in codes], color=[MOD_COLOUR[m] for m in codes], width=0.7)
    for i, m in enumerate(codes):
        ax.text(i, rvals[m][0] + 0.01, f'{rvals[m][0]:.2f}', ha='center', va='bottom', fontsize=9)
        ax.text(i, 0.02, f'{rvals[m][1]:,}', ha='center', va='bottom', fontsize=7, color='white', rotation=90)
    mod_axis(ax, codes)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel(f'r, {a} vs others')
    ax.set_title('per modification (sites in bars)', fontsize=11)
    letter(ax, 'I')
    notes.append(f'figure S1: {a} ({len(ia)} samples) vs the other batches ({len(ib)} samples)')
    c.legends['FigS1_batch_comparison'] = f"""One batch against the rest: {a} ({len(ia)} samples) against the other batches ({len(ib)}).
A-H, per cohort site measured in >= 2 samples on each side: the mean corrected occupancy over the other batches'
samples (x) and over {a}'s (y); log colour, one scale for all panels; dashed y = x; n sites, Pearson r, RMSE.
I, r per modification, sites in the bars. Figure 5 makes the same comparison for every batch at once."""
    return fig


# ---------------------------------------------------------------- per-sample QC from cohort_extract.py
def load_sample_qc(c, qc_dir, notes):
    """Reads cohort_extract.py outputs (<s>.integrity.tsv, <s>.read_categories.tsv, <s>.read_metrics.tsv.gz)
    and records per-sample summaries for Table 1 and the covariates."""
    out = {'integrity': {}, 'categories': {}, 'metrics': {}}
    if not qc_dir or not os.path.isdir(qc_dir):
        return out
    for s in c.samples:
        p = os.path.join(qc_dir, f'{s}.integrity.tsv')
        if os.path.exists(p):
            t = read_table(p, numeric=('length_lo', 'length_hi', 'reads', 'reached_5p', 'median_covered_fraction'))
            out['integrity'][s] = numbers(t, ['length_lo', 'length_hi', 'reads', 'reached_5p', 'median_covered_fraction'])
            v = out['integrity'][s]
            sel = (v[:, 0] >= 1000) & (v[:, 1] <= 2000)
            if sel.any() and v[sel, 2].sum() > 0:
                c.note_sample(s, 'complete_5p', round(float(np.average(v[sel, 3], weights=v[sel, 2])), 4))
        p = os.path.join(qc_dir, f'{s}.read_categories.tsv')
        if os.path.exists(p):
            t = read_table(p, numeric=('reads',))
            out['categories'][s] = dict(zip(strarr(t, 'category'), numbers(t, ['reads']).ravel()))
        p = os.path.join(qc_dir, f'{s}.read_metrics.tsv.gz')
        if os.path.exists(p):
            t = read_table(p, numeric=('read_len', 'speed_bps', 'qscore', 'mapq'))
            cols = [k for k in ('read_len', 'speed_bps', 'qscore', 'mapq') if k in t.columns]
            out['metrics'][s] = {k: numbers(t, [k]).ravel() for k in cols}
            for k, lab in (('read_len', 'read_len_median'), ('speed_bps', 'speed_median'), ('qscore', 'qscore_median')):
                if k in out['metrics'][s]:
                    v = out['metrics'][s][k]
                    if np.isfinite(v).any():
                        c.note_sample(s, lab, round(float(np.nanmedian(v)), 3))
    got = {k: len(v) for k, v in out.items() if v}
    if got:
        notes.append('sample QC (cohort_extract.py): ' + ', '.join(f'{k} for {v} samples' for k, v in got.items()))
    return out


# ---------------------------------------------------------------- figure 10: what drives each data layer
def occupancy_pcs(c, code='a', n=3):
    try:
        rows, ms, pct, mods, nv = c.mods_matrix()
    except FileNotFoundError:
        return None
    if nv is None:
        return None
    sel = (mods == code) & ((np.nan_to_num(nv) >= 20) & np.isfinite(pct)).all(1)
    if sel.sum() < 50:
        return None
    y = logit(np.clip(pct[sel] / 100, 0.01, 0.99))
    sc, fr = pca(y, n)
    return ms, sc, fr, int(sel.sum())


def fig10_drivers(c, qc, notes):
    S = c.samples
    pos = {s: i for i, s in enumerate(S)}
    covs = c.covariates()
    pc1 = composition_pc1(c)
    if pc1 is not None:
        covs['composition PC1'] = ('num', pc1)
    layers = {}                                                     # name -> values aligned to S

    def put(name, samples, vals):
        v = np.full(len(S), np.nan)
        for s_, x_ in zip(samples, vals):
            if s_ in pos:
                v[pos[s_]] = x_
        layers[name] = v
    try:
        es, lx, _, _ = c.expression()
        keep = np.argsort(lx.var(axis=1))[::-1][:min(2000, len(lx))]
        sc, fr = pca(lx[keep], 5)
        for k in range(sc.shape[1]):
            put(f'expression PC{k + 1} ({100 * fr[k]:.0f}%)', es, sc[:, k])
    except FileNotFoundError:
        pass
    if getattr(c, 'isoform_pcs', None):
        ss_, pcs_ = c.isoform_pcs
        for k in range(min(2, pcs_.shape[1])):
            put(f'isoform PC{k + 1}', ss_, pcs_[:, k])
    op = occupancy_pcs(c)
    if op:
        ms, sc, fr, nsites = op
        for k in range(sc.shape[1]):
            put(f'm$^6$A occupancy PC{k + 1} ({100 * fr[k]:.0f}%)', ms, sc[:, k])
    for key, lab in (('sites_all', 'sites called'), ('polya_median_per_read_nuclear', 'poly(A) median'),
                     ('identity_median', 'read identity'), ('fail_frac_m6A', 'm$^6$A fail fraction')):
        v = np.array([c.summary.get(s, {}).get(key, np.nan) for s in S], float)
        if key == 'sites_all':
            try:
                _, mss, mpct, _, _ = c.mods_matrix()
                d = dict(zip(mss, (np.nan_to_num(mpct) >= 20).sum(0)))
                v = np.array([d.get(s, np.nan) for s in S], float)
            except FileNotFoundError:
                pass
        if np.isfinite(v).sum() >= 4:
            layers[lab] = v
    if c.stream:
        lam = np.array([null_stats(c.stream[s]['plin'][MOD_INDEX['a']])[0] if s in c.stream else np.nan for s in S])
        if np.isfinite(lam).sum() >= 4:
            layers['m$^6$A null λ'] = lam
    cov_names = list(covs)
    R2 = np.full((len(layers), len(cov_names)), np.nan)
    P = np.full_like(R2, np.nan)
    for i, (ln, lv) in enumerate(layers.items()):
        for j, cn in enumerate(cov_names):
            R2[i, j], P[i, j] = r2_on(lv, covs[cn])
    Q = bh(P.ravel()).reshape(P.shape) if np.isfinite(P).any() else P
    rows_out = [(ln, cn, R2[i, j], P[i, j], Q[i, j]) for i, ln in enumerate(layers) for j, cn in enumerate(cov_names)]
    tname = c.table('Fig10_drivers.tsv', ['layer', 'covariate', 'r2', 'p', 'q'],
                    [(ln.replace('$^6$', '6'), cn, r, p_, q_) for ln, cn, r, p_, q_ in rows_out])

    ncov = max(1, len(cov_names))
    wa, wb = PANEL[0] * 2.2, PANEL[0] * max(1.4, 0.19 * ncov)
    htop, hbot = PANEL[1] * max(1.5, 0.17 * len(layers)), PANEL[1] * 1.35
    fig = plt.figure(figsize=(max(wa + wb, 4 * PANEL[0]) + 1.8, htop + hbot + 1.2))
    outer = fig.add_gridspec(2, 1, height_ratios=[htop, hbot])
    top = outer[0].subgridspec(1, 2, width_ratios=[wa, wb])
    bot = outer[1].subgridspec(1, 4)
    ax = fig.add_subplot(top[0, 0])                                 # A composition
    comp = c.composition()
    if comp['names']:
        V = comp['values']
        if comp['kind'] == 'measured' and np.nanmax(np.nansum(V, 1)) <= 1.05 + 99 * (np.nanmax(V) > 1.5):
            scale = 100 if np.nanmax(V) <= 1.05 else 1
            bottom = np.zeros(len(S))
            cols = plt.cm.tab10(np.linspace(0, 1, max(10, len(comp['names']))))
            for j, nm in enumerate(comp['names']):
                ax.bar(range(len(S)), np.nan_to_num(V[:, j]) * scale, bottom=bottom, color=cols[j % len(cols)], width=0.8, label=nm)
                bottom += np.nan_to_num(V[:, j]) * scale
            ax.set_ylabel('% of leukocytes')
            ax.legend(fontsize=7, ncol=4, loc='lower left', bbox_to_anchor=(0, 1.02), handlelength=1)
            sample_axis(ax, S, c.batch, c.colours, fontsize=7)
        else:
            im = ax.pcolormesh(np.arange(len(S) + 1) - 0.5, np.arange(len(comp['names']) + 1) - 0.5, V.T, cmap='RdBu_r',
                               vmin=-2, vmax=2, edgecolors='face', linewidth=0.3)
            ax.set_yticks(range(len(comp['names'])), comp['names'], fontsize=8)
            ax.set_ylim(len(comp['names']) - 0.5, -0.5)
            sample_axis(ax, S, c.batch, c.colours, fontsize=7)
            ax.tick_params(length=0)
            colourbar(fig, im, ax, label='marker score (z)')
        ax.set_title(f'composition: {comp["kind"]} ({comp["source"]})', fontsize=9)
    else:
        empty(ax, 'no composition (no gene CPM, --composition or --composition-columns)')
    letter(ax, 'A')

    ax = fig.add_subplot(top[0, 1])                                 # B drivers
    if len(layers) and len(cov_names):
        im = ax.pcolormesh(np.arange(len(cov_names) + 1) - 0.5, np.arange(len(layers) + 1) - 0.5, np.nan_to_num(R2),
                           cmap=HEAT_CMAP, vmin=0, vmax=1, edgecolors='white', linewidth=0.6)
        for i in range(len(layers)):
            for j in range(len(cov_names)):
                if np.isfinite(R2[i, j]) and (R2[i, j] >= 0.25 or (np.isfinite(Q[i, j]) and Q[i, j] < 0.1)):
                    star = '**' if np.isfinite(Q[i, j]) and Q[i, j] < 0.1 else ('*' if P[i, j] < 0.05 else '')
                    ax.text(j, i, f'{R2[i, j]:.2f}{star}', ha='center', va='center', fontsize=6,
                            color='white' if R2[i, j] > 0.55 else 'black')
        ax.set_xticks(range(len(cov_names)), cov_names, rotation=60, ha='right', fontsize=7)
        ax.set_yticks(range(len(layers)), list(layers), fontsize=7)
        ax.set_ylim(len(layers) - 0.5, -0.5)
        ax.tick_params(length=0)
        for side in ('left', 'bottom'):
            ax.spines[side].set_visible(False)
        colourbar(fig, im, ax, label='R$^2$')
        ax.set_title('what each data layer follows (* p < 0.05, ** FDR < 0.1)', fontsize=9)
    else:
        empty(ax, 'no layers or covariates')
    letter(ax, 'B')

    ax = fig.add_subplot(bot[0, 0])                                  # C sex check
    mism = []
    try:
        es_, xi = c.gene_log(['XIST'])
        _, ylog = c.gene_log(SEX_Y)
        yv = np.nanmean(ylog, 0) if np.isfinite(ylog).any() else np.full(len(es_), np.nan)
        scol = next((col for col in c.sheet.columns if col.lower() in ('sex', 'gender')), None)
        rec = dict(zip(c.sheet['sample_id'], c.sheet[scol])) if scol else {}
        for i, s in enumerate(es_):
            r_ = str(rec.get(s, '')).upper()[:1]
            colr = {'F': '#D593BD', 'M': '#398DC3'}.get(r_, '0.5')
            ax.scatter(xi[0, i], yv[i], color=colr, s=28, zorder=3)
            pred = 'M' if yv[i] > xi[0, i] else 'F'
            if r_ in ('F', 'M') and r_ != pred:
                mism.append(s)
                ax.annotate(s, (xi[0, i], yv[i]), fontsize=7, color='#E00000', xytext=(3, 2), textcoords='offset points')
        lim = [np.nanmin([np.nanmin(xi), np.nanmin(yv)]) - 0.5, np.nanmax([np.nanmax(xi), np.nanmax(yv)]) + 0.5]
        ax.plot(lim, lim, ls=':', color='0.5', lw=0.8)
        ax.set_xlabel('XIST, log$_2$(CPM + 1)')
        ax.set_ylabel('chrY genes, mean log$_2$(CPM + 1)')
        ax.set_title('recorded sex: female pink, male blue' if scol else 'no sex column in the samplesheet', fontsize=8)
        if mism:
            notes.append(f'figure 10: recorded sex disagrees with XIST/chrY expression for {", ".join(mism)}')
    except FileNotFoundError:
        empty(ax, 'no gene CPM matrix')
    letter(ax, 'C')

    ax = fig.add_subplot(bot[0, 1])                                  # D genotype concordance
    calls = {}
    for s in S:
        r_ = c.snvs(s)
        if r_:
            cl_ = r_[1]
            known = {k: v for k, v in cl_.items() if v[1] == 'known'}
            calls[s] = known if len(known) >= 100 else cl_
    related = []
    if len(calls) >= 2:
        ss = [s for s in S if s in calls]
        n = len(ss)
        Mx = np.full((n, n), np.nan)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                A_, B_ = set(calls[ss[i]]), set(calls[ss[j]])
                if min(len(A_), len(B_)):
                    Mx[i, j] = len(A_ & B_) / min(len(A_), len(B_))
                    if i < j and Mx[i, j] > 0.8:
                        related.append((ss[i], ss[j], Mx[i, j]))
        cmap = HEAT_CMAP.copy()
        cmap.set_bad('0.85')
        im = ax.pcolormesh(np.arange(n + 1) - 0.5, np.arange(n + 1) - 0.5, np.ma.masked_invalid(Mx), cmap=cmap,
                           vmin=0, vmax=1, edgecolors='face', linewidth=0.2)
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_aspect('equal')
        ax.set_xticks(range(n), [s[-3:] for s in ss], fontsize=6, rotation=90)
        ax.set_yticks(range(n), [s[-3:] for s in ss], fontsize=6)
        ax.tick_params(length=0)
        colourbar(fig, im, ax, label='shared SNVs / smaller set')
        ax.set_title('genotype sharing, known SNVs', fontsize=8)
        c.table('Fig10_genotype_sharing.tsv', ['sample_a', 'sample_b', 'shared_fraction'],
                [(ss[i], ss[j], Mx[i, j]) for i in range(n) for j in range(i + 1, n)])
        if related:
            notes.append('figure 10: possible same-donor pairs (shared SNV fraction > 0.8): ' +
                         ', '.join(f'{a_}/{b_} ({v:.2f})' for a_, b_, v in related[:10]))
    else:
        empty(ax, 'needs variants/<ref>/<s>/<s>.variants.tsv.gz\nfor ≥ 2 samples')
    letter(ax, 'D')

    ax = fig.add_subplot(bot[0, 2])                                  # E expression PC1 vs its strongest covariate
    if len(layers) and len(cov_names) and np.isfinite(R2[0]).any():
        j = int(np.nanargmax(R2[0]))
        ln = list(layers)[0]
        kind, vals = covs[cov_names[j]]
        y = layers[ln]
        if kind == 'num':
            for i, s in enumerate(S):
                ax.scatter(vals[i], y[i], color=c.colours[c.batch[s]], s=28, zorder=3)
            ax.set_xlabel(cov_names[j])
        else:
            lv = sorted(set(vals))
            for k, l_ in enumerate(lv):
                sel = np.array([v_ == l_ for v_ in vals])
                ax.scatter(k + jitter(sel.sum(), 0.3, k), y[sel], color=[c.colours[c.batch[s]] for s, k_ in zip(S, sel) if k_], s=28)
            ax.set_xticks(range(len(lv)), lv, rotation=30, fontsize=8)
            ax.set_xlabel(cov_names[j])
        ax.set_ylabel(ln)
        stat_box(ax, [f'R$^2$ = {R2[0, j]:.2f}'], loc='upper left')
    else:
        empty(ax, 'nothing to show')
    letter(ax, 'E')

    ax = fig.add_subplot(bot[0, 3])                                  # F design balance
    if c.group:
        batches = sorted(set(c.batch.values()))
        lv = c.levels
        tab = np.array([[sum(1 for s in S if c.group.get(s) == g_ and c.batch[s] == b_) for b_ in batches] for g_ in lv])
        ax.pcolormesh(np.arange(len(batches) + 1) - 0.5, np.arange(len(lv) + 1) - 0.5, tab, cmap=HEAT_CMAP, vmin=0,
                      vmax=max(1, tab.max()), edgecolors='white', linewidth=1)
        for i in range(len(lv)):
            for j in range(len(batches)):
                ax.text(j, i, str(tab[i, j]), ha='center', va='center', fontsize=9,
                        color='white' if tab[i, j] > tab.max() / 2 else 'black')
        ax.set_xticks(range(len(batches)), batches, rotation=30, fontsize=8)
        ax.set_yticks(range(len(lv)), lv, fontsize=8)
        ax.set_ylim(len(lv) - 0.5, -0.5)
        ax.tick_params(length=0)
        for side in ('left', 'bottom'):
            ax.spines[side].set_visible(False)
        ax.set_title(f'{c.args.group} by batch', fontsize=9)
    else:
        batches = sorted(set(c.batch.values()))
        ax.bar(range(len(batches)), [sum(1 for s in S if c.batch[s] == b_) for b_ in batches],
               color=[c.colours[b_] for b_ in batches])
        ax.set_xticks(range(len(batches)), batches, rotation=30, fontsize=8)
        ax.set_ylabel('samples')
        ax.set_title('samples per batch (no --group)', fontsize=9)
    letter(ax, 'F')
    fig_batch_legend(fig, c.colours)
    mk = comp.get('markers', {})
    mtxt = '; '.join(f'{k}: {", ".join(v)}' for k, v in mk.items()) if mk else ''
    c.legends['Fig10_cohort_drivers'] = f"""What drives each data layer ({len(S)} samples). Whole-blood expression is dominated by the leukocyte
differential (O'Connell et al., Genomics 2023), and any layer measured in bulk (isoform usage, modification
occupancy, tails) inherits that mixture, so composition, depth and handling are shown against every layer here.
A, composition: {comp['kind']} from {comp['source']}.{' Marker scores are the mean z-score of log2(CPM + 1) over each type' + "'" + 's markers (' + mtxt + '): relative between samples, not proportions; replace them with CBC or deconvolution (CIBERSORTx LM22, etc.) through --composition.' if comp['kind'].startswith('marker') else ''}
B, R^2 of each layer (rows) on each covariate (columns): linear for numeric covariates, one-way ANOVA for
categorical ones. Rows: expression PCs (2,000 most variable genes, CPM without globin/chrM), isoform-fraction PCs
(figure 9B), m6A occupancy PCs (logit, sites measured in every sample), and per-sample summaries; the m6A null
lambda is figure 12's calibration of the gIVT null. * p < 0.05; ** BH FDR < 0.1 over the whole matrix. With
{len(S)} samples, R^2 below ~{max(0.1, 4 / max(len(S), 1)):.2f} is indistinguishable from noise.
C, sex check: XIST against the mean of chrY genes ({", ".join(SEX_Y)}); samples whose recorded sex is on the
wrong side of the diagonal are named in red. D, genotype sharing between samples: shared Clair3-RNA SNVs (PASS,
gnomAD-known when >= 100 per sample; RNA-editing sites excluded) over the smaller set; a pair near 1 is one
donor sequenced twice (or a swap). E, expression PC1 against its strongest covariate. F, design balance.
Table: tables/{tname}."""
    return fig


# ---------------------------------------------------------------- figure 11: integrity and signal
def fig11_integrity(c, qc, notes):
    S = c.samples
    fig, gs = figure(2, 4)
    ax = fig.add_subplot(gs[0, 0])                                  # A 5' completeness
    if qc['integrity']:
        for s, v in qc['integrity'].items():
            mid = np.sqrt(np.maximum(v[:, 0], 1) * np.minimum(v[:, 1], 2e4))
            ok = v[:, 2] >= 20
            ax.plot(mid[ok], 100 * v[ok, 3], color=c.colours[c.batch[s]], lw=0.9, marker='o', ms=2.5)
        ax.set_xscale('log')
        log_axis(ax.xaxis)
        ax.set_ylim(0, 100)
        ax.set_xlabel('transcript length (nt)')
        ax.set_ylabel("reads reaching the 5' end (%)")
        ax.set_title("within 50 nt of the isoform's 5' end", fontsize=9)
    else:
        empty(ax, "needs cohort_extract.py integrity\n(read assignments + catalogue GTF)")
    letter(ax, 'A')
    ax = fig.add_subplot(gs[0, 1])                                  # B read lengths
    drawn = False
    for s in S:
        v = qc['metrics'].get(s, {}).get('read_len')
        if v is None:
            p = os.path.join(c.args.read_stats or c.path('read_stats'), f'{s}.read_identity.tsv.gz')
            if os.path.exists(p):
                v = numbers(read_table(p, columns=['aligned_len'], numeric=('aligned_len',)), ['aligned_len']).ravel()
        if v is not None and len(v):
            ecdf(ax, v, color=c.colours[c.batch[s]], lw=0.8)
            drawn = True
    if drawn:
        ax.set_xscale('log')
        log_axis(ax.xaxis)
        ax.set_xlabel('read length (nt)')
        ax.set_ylabel('cumulative fraction')
    else:
        empty(ax, 'no read lengths')
    letter(ax, 'B')
    ax = fig.add_subplot(gs[0, 2])                                  # C read categories
    if qc['categories']:
        cats = sorted({k for d in qc['categories'].values() for k in d}, key=lambda k: -sum(d.get(k, 0) for d in qc['categories'].values()))
        cats = cats[:7]
        ss = [s for s in S if s in qc['categories']]
        bottom = np.zeros(len(ss))
        shades = ['#3B0000', '#7B3333', '#B4766F', '#C9A9A2', '#E8CFC4', '#5769AB', '#EDB157']
        for k, cat in enumerate(cats):
            v = np.array([qc['categories'][s].get(cat, 0) / max(1, sum(qc['categories'][s].values())) for s in ss]) * 100
            ax.bar(range(len(ss)), v, bottom=bottom, color=shades[k % len(shades)], width=0.8, label=CATEGORY_SHORT.get(cat, cat))
            bottom += v
        sample_axis(ax, ss, c.batch, c.colours, fontsize=7)
        ax.set_ylabel('reads (%)')
        ax.legend(fontsize=6.5, ncol=4, loc='lower left', bbox_to_anchor=(0, 1.02), handlelength=1)
    else:
        empty(ax, 'needs cohort_extract.py isoquant')
    letter(ax, 'C')
    for k, (key, lab, letr) in enumerate((('speed_bps', 'translocation speed (nt/s)', 'D'), ('qscore', 'basecall qscore', 'E'))):
        ax = fig.add_subplot(gs[0, 3] if k == 0 else gs[1, 0])
        ss = [s for s in S if key in qc['metrics'].get(s, {})]
        if ss:
            data = [qc['metrics'][s][key][np.isfinite(qc['metrics'][s][key])] for s in ss]
            bp = ax.boxplot([d if len(d) else [np.nan] for d in data], positions=range(len(ss)), widths=0.7, showfliers=False,
                            patch_artist=True, medianprops={'color': 'black', 'lw': 0.8})
            for patch, s in zip(bp['boxes'], ss):
                patch.set_facecolor(c.colours[c.batch[s]])
            sample_axis(ax, ss, c.batch, c.colours, fontsize=7)
            ax.set_ylabel(lab)
            if key == 'speed_bps':
                ax.axhline(130, color='0.4', ls=':', lw=0.8)
        else:
            empty(ax, f'needs cohort_extract.py bam ({key})')
        letter(ax, letr)
    ax = fig.add_subplot(gs[1, 1])                                  # F modkit thresholds
    thr = c.thresholds()
    if thr:
        ss = [s for s in S if s in thr]
        for k, code in enumerate(ALL_MODS):
            v = [thr[s].get(code, np.nan) for s in ss]
            strip_points(ax, k, v, MOD_COLOUR[code], width=0.5, s=12, floor=None)
        mod_axis(ax, ALL_MODS, fontsize=8)
        ax.set_ylabel('call threshold (sample-probs)')
        spread = {plain(m): np.nanmax([thr[s].get(m, np.nan) for s in ss]) - np.nanmin([thr[s].get(m, np.nan) for s in ss])
                  for m in ALL_MODS if any(m in thr[s] for s in ss)}
        big = {k: v for k, v in spread.items() if v > 0.03}
        if big:
            notes.append('figure 11: per-sample call thresholds differ by > 0.03 for ' + ', '.join(f'{k} ({v:.2f})' for k, v in big.items()))
    else:
        empty(ax, 'no <sample>.givt_flags.txt')
    letter(ax, 'F')
    ax = fig.add_subplot(gs[1, 2])                                  # G fail fraction
    if c.stream:
        ss = [s for s in S if s in c.stream]
        for k, code in enumerate(ALL_MODS):
            i = MOD_INDEX[code]
            v = [c.stream[s]['nfail'][i] / max(1, c.stream[s]['nfail'][i] + c.stream[s]['nvalid'][i]) for s in ss]
            strip_points(ax, k, 100 * np.array(v), MOD_COLOUR[code], width=0.5, s=12, floor=None)
        mod_axis(ax, ALL_MODS, fontsize=8)
        ax.set_ylabel('calls below threshold (%)')
    else:
        empty(ax, 'needs the stream pass (cov20_fpcorr.bed)')
    letter(ax, 'G')
    ax = fig.add_subplot(gs[1, 3])                                  # H model versions
    vcols = [k for k in ('basecall_model', 'mod_models', 'dorado_version') if k in c.sheet.columns]
    if vcols:
        idx = dict(zip(c.sheet['sample_id'], range(len(c.sheet))))
        vals = [[c.sheet[k].iloc[idx[s]] if s in idx else '' for k in vcols] for s in S]
        mat = np.zeros((len(S), len(vcols)))
        legend_txt = []
        for j, k in enumerate(vcols):
            u = sorted({v[j] for v in vals})
            for i, v in enumerate(vals):
                mat[i, j] = u.index(v[j])
            legend_txt.append(f'{k}: ' + ' | '.join(f'{n}={x}' for n, x in enumerate(u)))
        ax.pcolormesh(np.arange(len(vcols) + 1) - 0.5, np.arange(len(S) + 1) - 0.5, mat, cmap='Set2', vmin=0, vmax=7,
                      edgecolors='white', linewidth=0.5)
        for i in range(len(S)):
            for j in range(len(vcols)):
                ax.text(j, i, str(int(mat[i, j])), ha='center', va='center', fontsize=6)
        ax.set_xticks(range(len(vcols)), vcols, rotation=30, fontsize=7)
        ax.set_yticks(range(len(S)), S, fontsize=6)
        for lab, s in zip(ax.get_yticklabels(), S):
            lab.set_color(c.colours[c.batch[s]])
        ax.set_ylim(len(S) - 0.5, -0.5)
        ax.tick_params(length=0)
        mixed = [k for j, k in enumerate(vcols) if len(set(mat[:, j])) > 1]
        ax.set_title('basecalling versions' + (' (MIXED: ' + ', '.join(mixed) + ')' if mixed else ' (uniform)'), fontsize=8,
                     color='#E00000' if mixed else 'black')
        if mixed:
            notes.append('figure 11: samples were basecalled with different ' + ', '.join(mixed) +
                         '; the gIVT table is only right for the model it was built with')
        c._versions_text = '; '.join(legend_txt)
    else:
        empty(ax, 'no basecall_model / mod_models / dorado_version columns')
    letter(ax, 'H')
    fig_batch_legend(fig, c.colours)
    c.legends['Fig11_integrity_signal'] = f"""RNA integrity and signal-level QC, one line or box per sample (coloured by batch).
A, fraction of reads whose 5' end lies within 50 nt of their isoform's annotated 5' end, by transcript length (as
Nelson et al. 2026, Modulator, Fig. 1a): direct RNA reads run 3' to 5', so degradation and pore stalls show as
fewer complete reads on long transcripts; this also sets how many sites near 5' ends are measurable. B, read
length ECDF. C, IsoQuant read classifications (SQANTI categories; ISM includes 5'-truncated fragments).
D, translocation speed from the move table (bases / (moves x stride / 4 kHz)); dotted: the 130 nt/s of
rna004_130bps models. E, Dorado mean qscore per read. F, modkit call thresholds per sample and modification
(sample-probs 10th percentile, from givt_flags.txt): gIVT's false-positive rate is read at the threshold, so
samples at different thresholds are corrected with different tables. G, % of calls below the threshold
(N_fail / (N_valid + N_fail)) in the cov20 rows. H, basecaller, modification models and Dorado version per sample
(numbers index the values: {getattr(c, '_versions_text', 'n/a')}).
A-E need cohort_extract.py (bam: speed, qscore, length; isoquant: integrity, categories) in --sample-qc."""
    return fig


# ---------------------------------------------------------------- figure 12: calibration of the gIVT null
def null_threshold(n, fp, alpha):
    """Corrected occupancy (%) a site needs at n reads to reach P(X >= x | n, fp) <= alpha."""
    k = stats.binom.isf(alpha, n, fp) + 1
    return 100 * (k / n - fp)


def fp_quantile(hist, q):
    cum = np.cumsum(hist)
    if cum[-1] == 0:
        return np.nan
    j = np.searchsorted(cum, q * cum[-1])
    return float(FP_EDGES[min(j + 1, len(FP_EDGES) - 1)])


def fig12_calibration(c, notes, fdr=0.05):
    if not c.stream:
        notes.append('figure 12 skipped: no cov20_fpcorr.bed (stream pass)')
        return None
    S = [s for s in c.samples if s in c.stream]
    D = c.stream
    exact = sum(int(D[s]['exact_fp'].sum()) for s in S)
    clipped = sum(int(D[s]['clipped'].sum()) for s in S)
    approx = clipped > 0 and exact == 0
    fig, gs = figure(3, 4, heights=[PANEL[1] * 1.2, PANEL[1] * 1.1, PANEL[1] * 1.1])
    xe = NCOV_EDGES.astype(float).copy()
    xe[-1] = 1e4
    pooled_fun = sum(D[s]['funnel'] for s in S)
    pooled_fp = sum(D[s]['fp_hist'] for s in S)
    nn = np.unique(np.round(np.logspace(np.log10(20), 4, 120)).astype(np.int64))
    vmax = max(10, max(pooled_fun[MOD_INDEX[m]].max() for m in MAIN_MODS))
    last, axes = None, []
    for k, code in enumerate(MAIN_MODS):                            # A-D funnels
        ax = fig.add_subplot(gs[0, k])
        axes.append(ax)
        i = MOD_INDEX[code]
        h = np.ma.masked_equal(pooled_fun[i].T, 0)
        if h.count():
            last = ax.pcolormesh(xe, PCT_EDGES, h, cmap='viridis', norm=LogNorm(vmin=1, vmax=vmax), edgecolors='face', linewidth=0.1)
        f50, f90 = fp_quantile(pooled_fp[i], 0.5), fp_quantile(pooled_fp[i], 0.9)
        for f_, ls, lab in ((f50, '-', f'median gIVT rate ({100 * f50:.2g}%)'), (f90, '--', f'90th percentile ({100 * f90:.2g}%)')):
            if np.isfinite(f_):
                ax.plot(nn, null_threshold(nn, f_, 1e-3), color='white', lw=1.6, ls=ls)
                ax.plot(nn, null_threshold(nn, f_, 1e-3), color='#E00000', lw=0.9, ls=ls, label=lab)
        ax.plot(nn, null_threshold(nn, 0.0121, 8.55e-5), color='black', lw=0.9, ls=':', label='paper: fp 1.21%, P < 8.5e-5')
        ax.axhline(20, color='0.4', lw=0.8, ls='-.')
        ax.scatter([20], [20], color='black', s=18, zorder=5)
        ax.set_xscale('log')
        ax.set_xlim(20, 1e4)
        ax.set_ylim(0, 100)
        ax.set_xlabel('reads at the site (N$_{valid}$)')
        ax.set_ylabel('corrected occupancy (%)')
        ax.set_title(MOD_NAME[code], color=MOD_COLOUR[code], fontweight='bold')
        if k == 0:
            ax.legend(fontsize=6.5, loc='upper right')
        letter(ax, 'ABCD'[k])
    if last is not None:
        colourbar(fig, last, axes, label='site rows (all samples)')

    ax = fig.add_subplot(gs[1, 0])                                  # E null p histograms, m6A, per sample
    i = MOD_INDEX['a']
    for s in S:
        h = D[s]['plin'][i]
        if h.sum():
            ax.step(np.linspace(0, 1, 51)[:-1], 50 * h / h.sum(), where='post', color=c.colours[c.batch[s]], lw=0.7)
    ax.axhline(1, color='black', lw=0.8, ls='--')
    ax.set_yscale('log')
    ax.set_xlabel('randomized null p (Binomial(N, fp))')
    ax.set_ylabel('density')
    ax.set_title('m$^6$A, per sample' + (' (APPROXIMATE: no --errortable)' if approx else ''), fontsize=9,
                 color='#E00000' if approx else 'black')
    letter(ax, 'E')
    lam = {s: {m: null_stats(D[s]['plin'][MOD_INDEX[m]]) for m in ALL_MODS} for s in S}
    for k, (j, lab) in enumerate(((0, 'genomic-control λ'), (1, 'π$_0$ (Storey, λ = 0.5)'))):
        ax = fig.add_subplot(gs[1, 1 + k])                          # F lambda, G pi0
        for t, m in enumerate(ALL_MODS):
            strip_points(ax, t, [lam[s][m][j] for s in S], MOD_COLOUR[m], width=0.5, s=12, floor=None)
        ax.axhline(1, color='0.5', lw=0.6, ls='--')
        mod_axis(ax, ALL_MODS, fontsize=8)
        ax.set_ylabel(lab)
        if j == 0:
            ax.set_yscale('log')
        letter(ax, 'FG'[k])
    for s in S:
        c.note_sample(s, 'null_lambda_m6A', round(lam[s]['a'][0], 3) if np.isfinite(lam[s]['a'][0]) else np.nan)

    ax = fig.add_subplot(gs[1, 3])                                  # H the 20% rule against FDR, by coverage
    mids = np.sqrt(xe[:-1] * xe[1:])
    rows_h = []
    for code in MAIN_MODS:
        i = MOD_INDEX[code]
        pl = sum(D[s]['plog'][i] for s in S)                         # (NCB, 2, NPB)
        t, _ = bh_from_hist(pl.sum((0, 1)), fdr)
        sigbins = P_LOG_EDGES[1:] <= t
        sig_called = pl[:, 1][:, sigbins].sum(1)
        called = pl[:, 1].sum(1)
        sig_uncalled = pl[:, 0][:, sigbins].sum(1)
        with np.errstate(invalid='ignore', divide='ignore'):
            f1 = np.where(called >= 20, sig_called / called, np.nan)
            f2 = np.where(sig_called + sig_uncalled >= 20, sig_uncalled / (sig_called + sig_uncalled), np.nan)
        ax.plot(mids, 100 * f1, color=MOD_COLOUR[code], lw=1.2, marker='o', ms=3, label=MOD_NAME[code])
        ax.plot(mids, 100 * f2, color=MOD_COLOUR[code], lw=0.9, ls='--')
        rows_h += [(plain(code), int(NCOV_EDGES[b]), int(called[b]), int(sig_called[b]), int(sig_uncalled[b])) for b in range(NCB)]
    ax.set_xscale('log')
    ax.set_ylim(0, 100)
    ax.set_xlabel('reads at the site')
    ax.set_ylabel('%')
    ax.set_title(f'solid: 20% calls significant at FDR {fdr:g}\ndashed: significant rows below 20%', fontsize=8)
    ax.legend(fontsize=7)
    letter(ax, 'H')

    ax = fig.add_subplot(gs[2, 0:2])                                # I calls per sample: 20% rule against FDR
    rows_i = []
    for k, code in enumerate(MAIN_MODS + OME_MODS):
        i = MOD_INDEX[code]
        r20, both, extra_ = [], [], []
        for s in S:
            pl = D[s]['plog'][i]
            t, _ = bh_from_hist(pl.sum((0, 1)), fdr)
            sb = P_LOG_EDGES[1:] <= t
            r20.append(pl[:, 1].sum())
            both.append(pl[:, 1][:, sb].sum())
            extra_.append(pl[:, 0][:, sb].sum())
            rows_i.append((s, plain(code), int(pl[:, 1].sum()), int(pl[:, 1][:, sb].sum()), int(pl[:, 0][:, sb].sum()), t))
        paired_points(ax, k, r20, both, MOD_COLOUR[code])
        ax.text(k + 0.3, max(np.median(r20), 1), f'+{int(np.median(extra_)):,}', ha='left', va='center', fontsize=7, color='0.3')
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    mod_axis(ax, ALL_MODS)
    ax.set_xlim(-0.6, len(ALL_MODS) - 0.4)
    ax.set_ylabel('sites per sample')
    ax.set_title(f'open: ≥ 20 reads and ≥ 20%; filled: of those, significant at FDR {fdr:g}; '
                 f'grey: median significant sites below 20%', fontsize=9)
    letter(ax, 'I')
    ax = fig.add_subplot(gs[2, 2])                                  # J distribution of gIVT rates
    for code in ALL_MODS:
        h = pooled_fp[MOD_INDEX[code]]
        if h.sum():
            ax.step(FP_EDGES[1:], np.cumsum(h) / h.sum(), where='pre', color=MOD_COLOUR[code], lw=1.1,
                    ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_xscale('log')
    ax.set_xlim(1e-5, 1)
    ax.set_xlabel('gIVT false-positive rate of the row (fp)')
    ax.set_ylabel('cumulative fraction of rows')
    ax.legend(fontsize=6.5, ncol=2)
    letter(ax, 'J')
    ax = fig.add_subplot(gs[2, 3])                                  # K clipped and uncorrectable rows
    for k, (key, lab) in enumerate((('clipped', 'clipped to 0'), ('uncorr', 'not correctable'))):
        v = []
        for s in S:
            rows_ = D[s]['rows'].sum()
            v.append(100 * (D[s]['clipped'].sum() if key == 'clipped' else rows_ - D[s]['correctable'].sum()) / max(1, rows_))
        strip_points(ax, k, v, '0.3', width=0.5, s=14, floor=None)
    ax.set_xticks([0, 1], ['clipped to 0', 'not correctable'], fontsize=8)
    ax.set_ylabel('% of cov20 rows, per sample')
    letter(ax, 'K')
    fig_batch_legend(fig, c.colours)
    t1 = c.table('Fig12_fdr_by_coverage.tsv', ['mod', 'coverage_bin_from', 'calls_20pct', 'calls_20pct_significant',
                                               'significant_below_20pct'], rows_h)
    t2 = c.table('Fig12_calls_per_sample.tsv', ['sample', 'mod', 'calls_20pct', 'calls_20pct_significant',
                                                'significant_below_20pct', 'bh_p_threshold'], rows_i)
    c.legends['Fig12_ivt_calibration'] = f"""Does the gIVT correction behave like a null, and what does the 20% rule do at each depth? ({len(S)} samples)
Every cov20 row (position x modification with >= 20 valid reads) is tested against its own null: N_mod ~
Binomial(N_valid, fp), fp = the gIVT rate of the row's 9-mer at the sample's call threshold (recovered from
columns 11 and 19; {exact:,} rows that gIVT clipped to 0 got their exact rate from --errortable{'' if exact else ', none here: E-G are APPROXIMATE, those rows are left out of them'}).
A-D, POOLED funnels: corrected occupancy against reads; lines: the occupancy needed for P(X >= x) <= 1e-3 at the
median (solid) and 90th-percentile (dashed) gIVT rate of the modification; dotted: the paper's argument (highest
oligo false-positive rate 1.21%, which gives P < 8.5e-5 at 20 reads and 20%, black dot) carried across depth;
dash-dot: the 20% rule. At low depth 20% is close to the null; at high depth far smaller changes are significant.
E, PER SAMPLE (m6A): histogram of randomized null p-values (P(X > x) + U P(X = x), uniform under a correct
null whatever the discreteness): flat with a spike at 0 is a calibrated null; a tilt means the gIVT rates do not
transfer to that sample (overdispersion, threshold, chemistry). F, genomic-control lambda (median chi-square over
its null median; 1 = calibrated) and G, pi0 (fraction of null rows) per sample and modification.
H, POOLED by coverage: % of 20% calls significant at BH FDR {fdr:g} (solid) and % of significant rows below 20%
(dashed). I, PER SAMPLE: 20% calls (open), those significant at FDR {fdr:g} (filled), and the median number of
significant rows below 20% (grey numbers). J, the gIVT rates the rows were corrected with. K, rows clipped to 0
and rows gIVT could not correct. Tables: tables/{t1}, tables/{t2}."""
    notes.append(f'figure 12: null p-values {"approximate (no --errortable)" if approx else "exact"}; {clipped:,} clipped rows')
    return fig


# ---------------------------------------------------------------- figure 13: controls and plausibility
def load_tx_models(gtf, cache_dir):
    """Canonical protein-coding transcript models from a GENCODE GTF (tag Ensembl_canonical; otherwise the
    longest CDS per gene): per chrom, exons sorted by start with their transcript's strand, transcript-coordinate
    offset, transcript length and CDS bounds in transcript coordinates. Cached in <out>/.cache."""
    import pickle
    st = os.stat(gtf)
    cp = os.path.join(cache_dir, f'tx_models.{hashlib.md5((os.path.abspath(gtf) + str(st.st_size) + str(int(st.st_mtime))).encode()).hexdigest()}.pkl')
    if os.path.exists(cp):
        with open(cp, 'rb') as fh:
            return pickle.load(fh)
    rx_t = re.compile(r'transcript_id "([^"]+)"')
    rx_g = re.compile(r'gene_id "([^"]+)"')
    tx = {}
    canonical_seen = False
    with opener(gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 8)
            if len(f) < 9 or f[2] not in ('exon', 'CDS'):
                continue
            a = f[8]
            if 'protein_coding' not in a:
                continue
            can = 'Ensembl_canonical' in a
            canonical_seen |= can
            m = rx_t.search(a)
            if not m:
                continue
            t = tx.setdefault(m.group(1), {'chrom': f[0], 'strand': f[6], 'exons': [], 'cds': [], 'can': can,
                                           'gene': (rx_g.search(a) or m).group(1)})
            (t['exons'] if f[2] == 'exon' else t['cds']).append((int(f[3]) - 1, int(f[4])))
    if canonical_seen:
        keep = {k: v for k, v in tx.items() if v['can'] and v['cds']}
    else:
        best = {}
        for k, v in tx.items():
            if v['cds']:
                L = sum(e - s for s, e in v['cds'])
                if L > best.get(v['gene'], (0, None))[0]:
                    best[v['gene']] = (L, k)
        keep = {k: tx[k] for _, k in best.values()}
    per = {}
    for k, v in keep.items():
        ex = sorted(v['exons'])
        lens = np.array([e - s for s, e in ex])
        cum = np.concatenate([[0], np.cumsum(lens)[:-1]])
        L = int(lens.sum())

        def tpos(g):
            for (s_, e_), c0 in zip(ex, cum):
                if s_ <= g < e_:
                    p = c0 + (g - s_)
                    return p if v['strand'] == '+' else L - 1 - p
            return None
        lo, hi = min(s for s, _ in v['cds']), max(e for _, e in v['cds']) - 1
        a_, b_ = tpos(lo), tpos(hi)
        if a_ is None or b_ is None:
            continue
        cs, ce = (a_, b_ + 1) if v['strand'] == '+' else (b_, a_ + 1)
        if cs <= 0 or ce >= L:
            continue
        for (s_, e_), c0 in zip(ex, cum):
            per.setdefault(v['chrom'], []).append((s_, e_, v['strand'], c0, L, cs, ce))
    models = {}
    for ch, rows in per.items():
        rows.sort()
        a = np.array([(r[0], r[1], 1 if r[2] == '+' else -1, r[3], r[4], r[5], r[6]) for r in rows], np.int64)
        models[ch] = a
    with open(cp, 'wb') as fh:
        pickle.dump(models, fh)
    return models


def metagene_coord(models, chroms, starts, strands):
    """Metagene position (0-1 5'UTR, 1-2 CDS, 2-3 3'UTR) of each site on a canonical transcript, NaN otherwise."""
    out = np.full(len(starts), np.nan)
    chroms = np.asarray(chroms, dtype=object)
    starts = np.asarray(starts, np.int64)
    sd = np.where(np.asarray(strands, dtype=object) == '-', -1, 1)
    for ch in set(chroms.tolist()):
        if ch not in models:
            continue
        a = models[ch]
        sel = np.where(chroms == ch)[0]
        j0 = np.searchsorted(a[:, 0], starts[sel], side='right') - 1
        done = np.zeros(len(sel), bool)
        for back in range(4):
            j = j0 - back
            okj = (j >= 0) & ~done
            jj = np.clip(j, 0, None)
            hit = okj & (starts[sel] < a[jj, 1]) & (a[jj, 2] == sd[sel])
            if hit.any():
                r = a[jj[hit]]
                p = r[:, 3] + (starts[sel][hit] - r[:, 0])
                p = np.where(r[:, 2] == 1, p, r[:, 4] - 1 - p)
                L, cs, ce = r[:, 4], r[:, 5], r[:, 6]
                x = np.where(p < cs, p / cs, np.where(p < ce, 1 + (p - cs) / (ce - cs), 2 + (p - ce) / np.maximum(L - ce, 1)))
                out[sel[hit]] = x
                done[hit] = True
    return out


KNOWN_ALIAS = {'m6a': 'a', 'm5c': 'm', 'psi': '17802', 'pseu': '17802', 'ψ': '17802', 'ino': '17596',
               'inosine': '17596', 'am': '69426', 'cm': '19228', 'um': '19227', 'gm': '19229'}


def parse_known(specs):
    out = []
    for sp in specs or []:
        name, sep, path = sp.partition('=')
        if not sep:
            sys.exit(f'cohort_plots.py: --known-sites expects MOD=BED, got {sp}')
        code = name if name in MOD_NAME else KNOWN_ALIAS.get(name.lower())
        if code is None:
            sys.exit(f'cohort_plots.py: --known-sites: unknown modification {name}')
        out.append((code, os.path.basename(path).split('.')[0], path))
    return out



def known_panels(c, rows, ms, pct, mods):
    """Overlap and occupancy panels for the --known-sites sets, one pair per modification."""
    panels = []
    known = parse_known(c.args.known_sites)
    if not known or rows is None:
        return panels
    loaded = {}                                                 # path -> intervals, read once

    def iv_of(path):
        if path not in loaded:
            loaded[path] = bed_intervals(path)
        return loaded[path]
    by_mod = {}
    for code, name, path in known:
        by_mod.setdefault(code, []).append((name, path))
    styles = [('-', 'o'), ('--', 's'), (':', '^'), ('-.', 'D'), ((0, (5, 1, 1, 1)), 'v'), ((0, (1, 2)), 'P')]
    overlap_rows, occ_rows = [], []

    def known_overlap(code):
        def draw(ax):
            edges = [20, 30, 40, 50, 60, 70, 80, 90, 101]
            xx = np.arange(len(edges) - 1)
            sel = np.where(mods == code)[0]
            hits = []
            for k, (name, path) in enumerate(by_mod[code]):
                hit, _ = point_in_bed(iv_of(path), rows['chrom'][sel], rows['start'][sel], rows['strand'][sel])
                hits.append(hit >= 0)
            sets = list(by_mod[code]) + ([('any set', None)] if len(by_mod[code]) > 1 else [])
            if len(by_mod[code]) > 1:
                hits.append(np.any(hits, axis=0))
            for k, ((name, _), inset) in enumerate(zip(sets, hits)):
                med = []
                for lo, hi in zip(edges[:-1], edges[1:]):
                    per = []
                    for j, s in enumerate(ms):
                        b_ = (np.nan_to_num(pct[sel, j]) >= lo) & (np.nan_to_num(pct[sel, j]) < hi)
                        if b_.sum() >= 10:
                            per.append(100 * inset[b_].mean())
                            overlap_rows.append((name, plain(code), s, lo, int(b_.sum()), per[-1]))
                    med.append((np.median(per), np.min(per), np.max(per)) if per else (np.nan,) * 3)
                med = np.array(med)
                anyset = name == 'any set'
                ls, mk = ('-', None) if anyset else styles[k % len(styles)]
                col = 'black' if anyset else MOD_COLOUR[code]
                ax.plot(xx, med[:, 0], color=col, lw=1.6 if anyset else 1.0, ls=ls, marker=mk, ms=3, label=name)
                if anyset or len(sets) == 1:
                    ax.fill_between(xx, med[:, 1], med[:, 2], color=col, alpha=0.12, lw=0)
            ax.set_xticks(range(len(edges) - 1), [str(e) for e in edges[:-1]], fontsize=8)
            ax.set_ylim(0, 100)
            ax.set_xlabel(f'{MOD_NAME[code]} corrected occupancy (%, bin start)')
            ax.set_ylabel('calls in the known set (%)')
            ax.legend(fontsize=6, loc='best')
            ax.set_title(f'{MOD_NAME[code]}: median over samples' + ('; band: any set, range' if len(sets) > 1 else ''), fontsize=8)
        return draw

    def known_occupancy(code):
        def draw(ax):
            sel = np.where(mods == code)[0]
            best = None
            for name, path in by_mod[code]:
                hit, val = point_in_bed(iv_of(path), rows['chrom'][sel], rows['start'][sel], rows['strand'][sel])
                ok = (hit >= 0) & np.isfinite(val)
                if ok.sum() >= 20 and (best is None or ok.sum() > best[0]):
                    best = (int(ok.sum()), name, ok, val)
            if best is None:
                empty(ax, f'no {plain(code)} set with occupancy in column 5\noverlapping ≥ 20 sites')
                return
            n_, name, ok, val = best
            v = val[ok] * (100 if np.nanmax(val[ok]) <= 1.0 else 1)
            with np.errstate(invalid='ignore'):
                mean = np.nanmean(pct[sel][ok], 1)
            m = hist2d(ax, mean, v, f'DRS {MOD_NAME[code]} %, mean over samples', f'{name} %', cmap='magma')
            if m is not None:
                colourbar(ax.figure, m, ax, label='sites')
                okm = np.isfinite(mean)
                occ_rows.append((plain(code), name, int(okm.sum()), float(np.corrcoef(mean[okm], v[okm])[0, 1]),
                                 float(stats.spearmanr(mean[okm], v[okm])[0])))
        return draw
    for code in [m for m in ALL_MODS if m in by_mod]:
        panels.append((f'known_{code}', known_overlap(code), 1))
        if any(np.isfinite(v[2]).any() for _, path in by_mod[code] for v in iv_of(path).values()):
            panels.append((f'known_occ_{code}', known_occupancy(code), 1))
    c._known_tables = (overlap_rows, occ_rows)
    return panels


def fig13_controls(c, notes):
    panels = []
    D = c.stream
    S = c.samples
    try:
        rows, ms, pct, mods, nv = c.mods_matrix()
    except FileNotFoundError:
        rows, ms, pct, mods, nv = None, [], None, None, None
    mt = {}
    for s in S:
        if s in D and len(D[s]['mito']):
            for pos, md, minus, n, x, rawp, corr, fp in D[s]['mito']:
                if not minus:
                    mt.setdefault((int(pos), int(md)), {})[s] = (n, rawp, corr)

    def mito_sites(ax):
        xs = []
        for k, (pos, code, lab, rna, base) in enumerate(MT_RRNA_SITES):
            if code is None:
                continue
            xs.append((k, lab))
            vals = mt.get((pos, MOD_INDEX[code]), {})
            for s, (n, rawp, corr) in vals.items():
                ax.scatter(len(xs) - 1 + jitter(1, 0.3, hash(s) % 97)[0], corr, color=c.colours[c.batch[s]], s=16, zorder=3)
            if vals:
                ax.hlines(np.nanmedian([v[2] for v in vals.values()]), len(xs) - 1.3, len(xs) - 0.7, color='black', lw=1.2)
        ax.set_xticks(range(len(xs)), [f'{lab}\n{MOD_NAME[MT_RRNA_SITES[k][1]]}' for k, lab in xs], fontsize=7)
        ax.set_ylim(-2, 102)
        ax.set_ylabel('corrected occupancy (%)')
        ax.set_title('mt-rRNA sites, expected model, per sample', fontsize=9)
        c.table('Fig13_mt_rrna_sites.tsv', ['mt_position', 'site', 'rRNA', 'model', 'sample', 'n_valid', 'raw_pct', 'corrected_pct'],
                [(pos, lab.replace('$', '').replace('^', '').replace('_', ''), rna, plain(ALL_MODS[md]), s, n, rp, cr)
                 for pos, code, lab, rna, base in MT_RRNA_SITES for md in range(NM)
                 for s, (n, rp, cr) in mt.get((pos, md), {}).items()])

    def mito_specificity(ax):
        Mx = np.full((NM, len(MT_RRNA_SITES)), np.nan)
        for j, (pos, code, lab, rna, base) in enumerate(MT_RRNA_SITES):
            for i in range(NM):
                v = mt.get((pos, i), {})
                if v:
                    Mx[i, j] = np.nanmedian([x[1] for x in v.values()])
        cmap = HEAT_CMAP.copy()
        cmap.set_bad('0.92')
        ax.pcolormesh(np.arange(len(MT_RRNA_SITES) + 1) - 0.5, np.arange(NM + 1) - 0.5, np.ma.masked_invalid(Mx), cmap=cmap,
                      vmin=0, vmax=100, edgecolors='white', linewidth=0.5)
        for j, (pos, code, lab, rna, base) in enumerate(MT_RRNA_SITES):
            for i in range(NM):
                if np.isfinite(Mx[i, j]):
                    ax.text(j, i, f'{Mx[i, j]:.0f}', ha='center', va='center', fontsize=6, color='white' if Mx[i, j] > 55 else 'black')
            if code is not None:
                ax.add_patch(plt.Rectangle((j - 0.5, MOD_INDEX[code] - 0.5), 1, 1, fill=False, ec='black', lw=1.2))
        ax.set_xticks(range(len(MT_RRNA_SITES)), [x[2] for x in MT_RRNA_SITES], rotation=60, fontsize=6.5, ha='right')
        ax.set_yticks(range(NM), [MOD_NAME[m] for m in ALL_MODS], fontsize=7)
        ax.set_ylim(NM - 0.5, -0.5)
        ax.tick_params(length=0)
        ax.set_title('raw Dorado % of every model (median over samples); box = expected', fontsize=8)

    have_mt = bool(mt)
    if have_mt:
        panels += [('mito_sites', mito_sites, 1), ('mito_spec', mito_specificity, 1)]
    xt = None
    if D:
        xs_ = sum(D[s]['xt_sum'] for s in S if s in D)
        xn_ = sum(D[s]['xt_n'] for s in S if s in D)
        bg = sum(D[s]['bg_sum'] for s in S if s in D) / np.maximum(1, sum(D[s]['bg_n'] for s in S if s in D))
        xt = (xs_, xn_, bg)

    def crosstalk(anchor):
        def draw(ax):
            xs_, xn_, bg = xt
            a = MOD_INDEX[anchor]
            for t, code in enumerate(ALL_MODS):
                with np.errstate(invalid='ignore', divide='ignore'):
                    y = (xs_[a, t] / np.maximum(xn_[a, t], 1)) / max(bg[t], 1e-9)
                y[xn_[a, t] < 20] = np.nan
                if t == a:
                    y[10] = np.nan
                if np.isfinite(y).any():
                    ax.plot(XT_OFFSETS, y, color=MOD_COLOUR[code], lw=1.1, marker='o', ms=2.5,
                            ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
            ax.axhline(1, color='0.5', lw=0.6)
            ax.axvspan(-4.5, 4.5, color='0.93', zorder=0)
            ax.set_yscale('log')
            log_axis(ax.yaxis)
            ax.set_xlabel(f'offset from a ≥ 50% {MOD_NAME[anchor]} site (nt, 5′ to 3′)')
            ax.set_ylabel('raw occupancy / background')
            ax.legend(fontsize=6, ncol=2, loc='upper right')
        return draw
    if xt is not None:
        for anchor in ('a', '17802', 'm', '17596'):
            if xt[1][MOD_INDEX[anchor]].sum() > 0:
                panels.append((f'xt_{anchor}', crosstalk(anchor), 1))

    if c.fasta.ok and rows is not None:
        def drach(ax):
            sel = np.where(mods == 'a')[0]
            ctx = c.fasta.context(rows['chrom'][sel], rows['start'][sel], rows['strand'][sel], flank=2)
            dr = is_drach(ctx)
            edges = [-1, 5, 20, 30, 40, 50, 60, 70, 80, 90, 101]
            labs = ['<5', '5-20'] + [str(e) for e in edges[2:-1]]
            out_rows = []
            for s in ms:
                j = ms.index(s)
                p_ = pct[sel, j]
                meas = np.isfinite(p_) & (np.nan_to_num(nv[sel, j]) >= 20) if nv is not None else np.isfinite(p_)
                fr = []
                for lo, hi in zip(edges[:-1], edges[1:]):
                    b_ = meas & (p_ >= max(lo, 0)) & (p_ < hi)
                    fr.append(100 * dr[b_].mean() if b_.sum() >= 20 else np.nan)
                    out_rows.append((s, labs[len(fr) - 1], int(b_.sum()), fr[-1]))
                ax.plot(range(len(labs)), fr, color=c.colours[c.batch[s]], lw=0.8, marker='o', ms=2.5)
            ax.set_xticks(range(len(labs)), labs, fontsize=7, rotation=45)
            ax.set_ylim(0, 100)
            ax.set_xlabel('m$^6$A corrected occupancy (%, bin)')
            ax.set_ylabel('DRACH (%)')
            ax.set_title('cohort m$^6$A sites, per sample', fontsize=9)
            c.table('Fig13_drach_by_occupancy.tsv', ['sample', 'occupancy_bin', 'sites', 'drach_pct'], out_rows)
        panels.append(('drach', drach, 1))
    if c.args.gtf and os.path.exists(c.args.gtf) and rows is not None:
        def metagene(ax):
            models = load_tx_models(c.args.gtf, os.path.join(c.out, '.cache'))
            x = metagene_coord(models, rows['chrom'], rows['start'], rows['strand'])
            called_any = (np.nan_to_num(pct) >= 20).any(1)
            edges = np.linspace(0, 3, 61)
            for code in MAIN_MODS:
                s_ = (mods == code) & called_any & np.isfinite(x)
                if s_.sum() >= 20:
                    h, _ = np.histogram(x[s_], bins=edges, density=True)
                    ax.plot((edges[:-1] + edges[1:]) / 2, h, color=MOD_COLOUR[code], lw=1.3, label=f'{MOD_NAME[code]} ({s_.sum():,})')
            for b_ in (1, 2):
                ax.axvline(b_, color='0.5', lw=0.6, ls='--')
            ax.set_xticks([0.5, 1.5, 2.5], ["5′UTR", 'CDS', "3′UTR"])
            ax.set_ylabel('density of called sites')
            ax.set_title('canonical protein-coding transcripts', fontsize=9)
            ax.legend(fontsize=7)
        panels.append(('metagene', metagene, 1))
    if c.args.repeats and os.path.exists(c.args.repeats) and rows is not None:
        def alu(ax):
            iv = bed_intervals(c.args.repeats, name_filter=lambda n: n.startswith('Alu'))
            sel = np.where(mods == '17596')[0]
            hit, _ = point_in_bed({k: (v[0], v[1], v[2], np.array(['.'] * len(v[0]))) for k, v in iv.items()},
                                  rows['chrom'][sel], rows['start'][sel])
            inA = hit >= 0
            for j, s in enumerate(ms):
                p_ = pct[sel, j]
                called = np.nan_to_num(p_) >= 20
                bgm = np.isfinite(p_) & (np.nan_to_num(p_) < 10)
                if called.sum() >= 5 and bgm.sum() >= 5:
                    ax.plot([0, 1], [100 * inA[bgm].mean(), 100 * inA[called].mean()], color=c.colours[c.batch[s]], lw=0.8, marker='o', ms=4)
            ax.set_xticks([0, 1], ['measured, < 10%', 'called, ≥ 20%'])
            ax.set_xlim(-0.4, 1.4)
            ax.set_ylabel('inosine sites in Alu (%)')
        panels.append(('alu', alu, 1))
    if D and any(D[s]['snp_rows'].sum() for s in S if s in D):
        def snp(ax):
            hs = sum(D[s]['hist_snp'] for s in S if s in D)
            hn = sum(D[s]['hist_nosnp'] for s in S if s in D)
            mids = (PCT_EDGES[:-1] + PCT_EDGES[1:]) / 2
            for code in MAIN_MODS:
                i = MOD_INDEX[code]
                if hs[i].sum() >= 20:
                    ax.plot(mids, hs[i] / hs[i].sum(), color=MOD_COLOUR[code], lw=1.2, label=f'{MOD_NAME[code]}, SNV in 9-mer')
                if hn[i].sum():
                    ax.plot(mids, hn[i] / hn[i].sum(), color=MOD_COLOUR[code], lw=0.8, ls=':')
            ax.axvline(50, color='0.5', lw=0.6, ls='--')
            ax.set_xlim(20, 100)
            ax.set_xlabel('corrected occupancy of 20% calls (%)')
            ax.set_ylabel('fraction of calls')
            frac = [100 * D[s]['snp_called'][MOD_INDEX['a']] / max(1, D[s]['called20'][MOD_INDEX['a']]) for s in S if s in D]
            ax.set_title(f'donor SNVs within ±4 nt (dotted: none); m$^6$A calls with one: {np.median(frac):.2f}% (median)', fontsize=8)
            ax.legend(fontsize=6.5)
        panels.append(('snp', snp, 1))
    if not panels:
        notes.append('figure 13 skipped: no mt-rRNA rows, cross-talk data, --reference, --gtf, --repeats or variants')
        return None
    ncol = 4
    nrow = math.ceil(len(panels) / ncol)
    fig, gs = figure(nrow, ncol)
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for k, (name, draw, _) in enumerate(panels):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        draw(ax)
        letter(ax, letters[k])
    fig_batch_legend(fig, c.colours)
    desc = {
        'mito_sites': 'mitochondrial rRNA sites with a Dorado model (16S Ψ1397 = m.3067, Gm1145 = m.2815, Um1369 = m.3039, '
                      'Gm1370 = m.3040; 12S m5C841 = m.1488): corrected occupancy per sample (points), median (line). These are '
                      'near-constitutive in mature mitoribosomes, so they are a within-sample positive control and their spread '
                      'between samples is a calibration track (check chrM variants and heteroplasmy at them first)',
        'mito_spec': 'every model\'s raw occupancy (median over samples) at all known mt-rRNA modifications, including those no '
                     'model is trained for (m5U, m4C, m6,6A, m1A): off-diagonal values are cross-reactions. m4C839 sits two '
                     'bases from m5C841',
        'xt_a': 'cross-talk: raw occupancy of every model at -10..+10 nt (RNA 5′ to 3′) around sites with ≥ 50% corrected m6A, over '
                'its background (mean raw occupancy at all of that model\'s rows), pooled over samples; grey: the 9-mer window. '
                'ONT\'s synthetic data show false m5C calls of 13.6% within 2 nt of an m6A against 0.9% further away; this '
                'measures the echo in these data and so the exclusion window each pair of models needs',
        'xt_17802': 'the same around ≥ 50% Ψ', 'xt_m': 'the same around ≥ 50% m5C', 'xt_17596': 'the same around ≥ 50% inosine',
        'drach': 'fraction of cohort m6A sites in a DRACH context by the sample\'s corrected occupancy (one line per sample; "<5" = '
                 'measured but not called there): real m6A is mostly DRACH, so a low fraction at high occupancy flags false calls',
        'metagene': 'metagene of called sites (union over samples) on canonical protein-coding transcripts: m6A should rise toward '
                    'the stop codon and 3′UTR start (Meyer et al. 2012)',
        **{f'known_{m}': f'{MOD_NAME[m].replace("$", "").replace("^", "")} calls found in each --known-sites set, by the call\'s occupancy, median over samples '
                         '(black: in any set). The sets come from other cell lines (HEK293T, HeLa), so this measures how many calls '
                         'sit at sites known to be modifiable, a lower bound on specificity in blood, not their correctness' for m in ALL_MODS},
        **{f'known_occ_{m}': f'{MOD_NAME[m].replace("$", "").replace("^", "")} occupancy against the overlapping set with the most shared sites '
                             '(as the paper\'s Fig. 6E); between cell types, stoichiometry differs at conserved sites, so r is a floor'
           for m in ALL_MODS},
        'alu': 'inosine sites in Alu elements: called (≥ 20%) against measured but uncalled (< 10%) sites of the same sample',
        'snp': 'occupancy of 20% calls with a donor SNV (Clair3-RNA PASS, editing excluded) within ±4 nt (solid) and without '
               '(dotted), pooled: the reference 9-mer\'s gIVT rate is wrong for the alternative allele\'s reads, and SNVs in the '
               'motif core shift stoichiometry by a median 47 points (Nelson et al. 2026), so piles at 50% and 100% flag artefacts',
    }
    c.legends['Fig13_controls_plausibility'] = 'Positive controls, cross-talk and biological plausibility of the calls.\n' + \
        '\n'.join(f'{letters[k]}, {desc.get(name, name)}.' for k, (name, _, _) in enumerate(panels)) + \
        '\nPanels appear when their input is there: --reference (DRACH), --gtf (metagene), --repeats, ' \
        'variant tables (SNVs); mt-rRNA and cross-talk come from the stream pass.'
    return fig


# ---------------------------------------------------------------- figure 16: orthogonal site sets
def fig16_orthogonal(c, notes):
    if not c.args.known_sites:
        notes.append('figure 16 skipped: no --known-sites (orthogonal_to_bed.py makes them)')
        return None
    rows, ms, pct, mods, nv = c.mods_matrix()
    panels = known_panels(c, rows, ms, pct, mods)
    if not panels:
        notes.append('figure 16 skipped: no known-site panel could be drawn')
        return None
    ncol = 4
    fig, gs = figure(math.ceil(len(panels) / ncol), ncol)
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for k, (name, draw, _) in enumerate(panels):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        draw(ax)
        letter(ax, letters[k])
    fig_batch_legend(fig, c.colours)
    orows, crows = c._known_tables
    t1 = c.table('Fig16_known_site_overlap.tsv', ['set', 'mod', 'sample', 'occupancy_bin_from', 'calls', 'pct_in_set'], orows)
    t2 = c.table('Fig16_known_site_occupancy.tsv', ['mod', 'set', 'sites', 'pearson_r', 'spearman_rho'], crows)
    sets = '; '.join(f'{MOD_NAME[code].replace("$", "").replace("^", "")}: {name}' for code, name, _ in parse_known(c.args.known_sites))
    kinds = ', '.join(f'{letters[k]}' for k, (n_, _, _) in enumerate(panels) if n_.startswith('known_occ'))
    c.legends['Fig16_orthogonal_sites'] = f"""Calls against orthogonal site sets ({sets}).
For each modification with a set: the fraction of the sample's calls (>= 20 reads, >= 20% corrected) that are in
each set, by the call's occupancy, median over samples (line style = set; black = in any set, band = its range
over samples){'; then, where a set reports occupancy, DRS occupancy (mean over samples) against the set with the most shared sites (panels ' + kinds + ', as the paper\'s Fig. 6E)' if kinds else ''}.
The sets are from HEK293T (and other cell lines), not blood, and differ among themselves (the paper: GLORI1 against
GLORI2 Jaccard 0.64, two GLORI1 datasets 0.37), so overlap is a lower bound on specificity and r a floor on
agreement: m6A sites are often shared between cell types while their stoichiometry is cell-type specific. Rising
overlap with occupancy is the expected pattern for real calls; a flat, low line at high occupancy flags false
positives. Sets were made with orthogonal_to_bed.py, which places every site on the modified base by checking
the reference. Tables: tables/{t1}, tables/{t2}."""
    return fig


# ---------------------------------------------------------------- figure 14: measurement against biology
def seed_of(s):
    return int(hashlib.md5(str(s).encode()).hexdigest()[:8], 16)


def split_halves(N, M, rng, frac=0.5):
    """Binomial read split of every (site, sample): (n1, m1, n2, m2); NaN where not measured."""
    ok = np.isfinite(N) & np.isfinite(M) & (np.nan_to_num(N) > 0)
    Ni = np.where(ok, N, 0).astype(np.int64)
    Mi = np.clip(np.where(ok, M, 0).astype(np.int64), 0, Ni)
    n1 = rng.binomial(Ni, frac)
    m1 = np.zeros_like(Ni)
    go = (Ni > 0) & (n1 > 0)
    m1[go] = rng.hypergeometric(Mi[go], Ni[go] - Mi[go], n1[go])
    f = lambda a: np.where(ok, a, np.nan).astype(float)            # noqa: E731
    return f(n1), f(m1), f(Ni - n1), f(Mi - m1)


def corrected_pct(m, n, fp, min_n=10):
    with np.errstate(invalid='ignore', divide='ignore'):
        v = np.clip(100 * (m / n - np.nan_to_num(fp)), 0, 100)
    return np.where(np.nan_to_num(n) >= min_n, v, np.nan)


def fig14_measurement(c, notes, max_sites=30000, max_pairs=300):
    try:
        rows, ms, pct, mods, nv = c.mods_matrix()
    except FileNotFoundError:
        notes.append('figure 14 skipped: no cohort modification matrices')
        return None
    nvc, nmc, fpc, _ = c.mods_counts()
    if nvc is None:
        notes.append('figure 14 skipped: needs cohort.<ref>.mods.n_valid and n_mod matrices')
        return None
    rng = np.random.default_rng(c.args.seed)
    k = len(ms)
    fig, gs = figure(2, 4, widths=[PANEL[0] * 1.6, PANEL[0], PANEL[0], PANEL[0]])
    ax = fig.add_subplot(gs[0, 0])                                  # A split-half vs donor pairs
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    if len(pairs) > max_pairs:
        pairs = [pairs[x] for x in rng.choice(len(pairs), max_pairs, replace=False)]
    out_rows, codes_done = [], []
    for t, code in enumerate(ALL_MODS):
        sel = np.where(mods == code)[0]
        if len(sel) < 50:
            continue
        if len(sel) > max_sites:
            sel = np.sort(rng.choice(sel, max_sites, replace=False))
        N, M, F, P = nvc[sel], nmc[sel], fpc[sel], pct[sel]
        called = (np.nan_to_num(N) >= 20) & (np.nan_to_num(P) >= 20)
        nc = called.sum(1)
        n1, m1, n2, m2 = split_halves(N, M, rng)
        o1, o2 = corrected_pct(m1, n1, F), corrected_pct(m2, n2, F)
        within, same, cross = [], [], []
        for j in range(k):
            s_ = np.isfinite(o1[:, j]) & np.isfinite(o2[:, j]) & ((nc - called[:, j]) >= 1)
            if s_.sum() >= 10:
                within.append(icc21(np.column_stack([o1[s_, j], o2[s_, j]])))
                out_rows.append((plain(code), 'within sample (split half)', ms[j], ms[j], within[-1], int(s_.sum())))
        for i, j in pairs:
            s_ = np.isfinite(o1[:, i]) & np.isfinite(o1[:, j]) & ((nc - called[:, i] - called[:, j]) >= 1)
            if s_.sum() >= 10:
                v = icc21(np.column_stack([o1[s_, i], o1[s_, j]]))
                sb = c.batch[ms[i]] == c.batch[ms[j]]
                (same if sb else cross).append(v)
                out_rows.append((plain(code), 'donors, same batch' if sb else 'donors, other batch', ms[i], ms[j], v, int(s_.sum())))
        x0 = len(codes_done)
        for off, vals, face in ((-0.25, within, MOD_COLOUR[code]), (0.0, same, 'white'), (0.25, cross, '0.6')):
            vals = np.array(vals, float)
            vals = vals[np.isfinite(vals)]
            if len(vals):
                ax.scatter(x0 + off + jitter(len(vals), 0.18, t), vals, s=7, facecolor=face, edgecolor=MOD_COLOUR[code], lw=0.6, zorder=3)
                ax.hlines(np.median(vals), x0 + off - 0.1, x0 + off + 0.1, color='black', lw=1.2, zorder=4)
        codes_done.append(code)
    mod_axis(ax, codes_done)
    ax.set_ylabel('ICC(2,1) at half depth')
    ax.set_ylim(top=1.02)
    ax.set_title('filled: split halves of one sample; open: donor pairs, same batch; grey: other batch', fontsize=8)
    letter(ax, 'A')
    tA = c.table('Fig14_split_half_icc.tsv', ['mod', 'comparison', 'sample_a', 'sample_b', 'icc21', 'sites'], out_rows)

    need = max(3, math.ceil(0.75 * k))
    ratios, ratio_rows = {}, []
    for code in ALL_MODS:
        sel = mods == code
        N, P, F = nvc[sel], pct[sel], fpc[sel]
        meas = np.isfinite(P) & (np.nan_to_num(N) >= 20)
        ok = meas.sum(1) >= need
        if ok.sum() < 20:
            continue
        Pm = np.where(meas[ok], P[ok], np.nan)
        with np.errstate(invalid='ignore', divide='ignore'):
            obs = np.nanstd(Pm, 1, ddof=1)
            th = np.clip(np.nanmean(Pm, 1) / 100 + np.nanmean(np.where(meas[ok], F[ok], np.nan), 1), 1e-4, 1 - 1e-4)
            exp = 100 * np.sqrt(np.nanmean(np.where(meas[ok], th[:, None] * (1 - th[:, None]) / N[ok], np.nan), 1))
        ratios[code] = (obs, exp)
        r_ = obs / exp
        ratio_rows.append((plain(code), int(ok.sum()), float(np.nanmedian(r_)), float(np.nanmean(r_ >= 2))))
    ax = fig.add_subplot(gs[0, 1])                                  # B observed vs binomial SD, m6A
    if 'a' in ratios:
        obs, exp = ratios['a']
        okk = np.isfinite(obs) & np.isfinite(exp) & (obs > 0) & (exp > 0)
        be = np.logspace(-0.5, 2, 40)
        h, xe_, ye_ = np.histogram2d(exp[okk], obs[okk], bins=[be, be])
        m = ax.pcolormesh(xe_, ye_, np.ma.masked_equal(h.T, 0), cmap='viridis', norm=LogNorm(vmin=1), edgecolors='face', linewidth=0.1)
        ax.plot(be, be, color='black', ls='--', lw=0.8)
        ax.plot(be, 2 * be, color='black', ls=':', lw=0.8)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('binomial SD (points)')
        ax.set_ylabel('observed SD across samples')
        stat_box(ax, [f'{okk.sum():,} sites', f'{100 * np.mean(obs[okk] / exp[okk] >= 2):.0f}% ≥ 2x'])
        colourbar(fig, m, ax, label='sites')
        ax.set_title('m$^6$A', color=MOD_COLOUR['a'], fontweight='bold')
    else:
        empty(ax, f'too few m6A sites measured in ≥ {need} samples')
    letter(ax, 'B')
    ax = fig.add_subplot(gs[0, 2])                                  # C per modification
    cd = [r for r in ratio_rows]
    if cd:
        codes = [next(m for m in ALL_MODS if plain(m) == r[0]) for r in cd]
        ax.bar(range(len(cd)), [100 * r[3] for r in cd], color=[MOD_COLOUR[m] for m in codes], width=0.7)
        for i, r in enumerate(cd):
            ax.text(i, 100 * r[3] + 1, f'{r[2]:.1f}x', ha='center', fontsize=7)
        mod_axis(ax, codes, fontsize=8)
        ax.set_ylabel('sites with SD ≥ 2x binomial (%)')
        ax.set_title('numbers: median SD ratio', fontsize=8)
    letter(ax, 'C')
    tC = c.table('Fig14_overdispersion.tsv', ['mod', 'sites', 'median_sd_ratio', 'fraction_ge2'], ratio_rows)

    ax = fig.add_subplot(gs[0, 3])                                  # D rarefaction, m6A
    depth = c.depth() if os.path.exists(os.path.join(c.cohort, f'cohort.{c.ref}.read_qc.tsv')) else {}
    rare_rows = []
    if c.stream:
        i = MOD_INDEX['a']
        for s, d in c.stream.items():
            f_ = np.array(list(THIN_FRACS) + [1.0])
            y = np.array(list(d['thin'][i]) + [d['called20'][i]], float)
            x = f_ * depth.get(s, np.nan) / 1e6 if np.isfinite(depth.get(s, np.nan)) else f_
            ax.plot(x, y, color=c.colours[c.batch[s]], lw=0.9, marker='o', ms=2.5)
            for code in ALL_MODS:
                j = MOD_INDEX[code]
                rare_rows += [(s, plain(code), float(ff), int(v)) for ff, v in zip(f_, list(c.stream[s]['thin'][j]) + [c.stream[s]['called20'][j]])]
        ax.set_xlabel('non-globin reads (M), thinned' if depth else 'fraction of reads')
        ax.set_ylabel('m$^6$A sites (≥ 20 reads, ≥ 20%)')
        ax.set_title('binomial thinning of each sample', fontsize=9)
    else:
        empty(ax, 'needs the stream pass')
    letter(ax, 'D')
    ax = fig.add_subplot(gs[1, 0])                                  # E rarefaction, all mods (relative)
    if c.stream:
        f_ = np.array(list(THIN_FRACS) + [1.0])
        for code in ALL_MODS:
            j = MOD_INDEX[code]
            rel = []
            for s, d in c.stream.items():
                full = d['called20'][j]
                if full >= 20:
                    rel.append(np.array(list(d['thin'][j]) + [full], float) / full)
            if rel:
                ax.plot(f_, np.median(rel, 0), color=MOD_COLOUR[code], lw=1.2, marker='o', ms=3, ls='-' if code in MAIN_MODS else '--',
                        label=MOD_NAME[code])
        ax.plot([0, 1], [0, 1], color='0.6', lw=0.6, ls=':')
        ax.set_xlabel('fraction of reads kept')
        ax.set_ylabel('sites / sites at full depth')
        ax.legend(fontsize=6.5, ncol=2)
        ax.set_title('median over samples; above the diagonal = saturating', fontsize=8)
    letter(ax, 'E')
    c.table('Fig14_rarefaction.tsv', ['sample', 'mod', 'fraction', 'sites'], rare_rows)

    ax = fig.add_subplot(gs[1, 1:3])                                # F reliability of one sample
    rel_rows, vio, vcodes = [], [], []
    X1 = np.ones((k, 1))
    tau_med = {}
    for code in ALL_MODS:
        sel = np.where(mods == code)[0]
        if len(sel) < 50:
            continue
        if len(sel) > 20000:
            sel = np.sort(rng.choice(sel, 20000, replace=False))
        N = nvc[sel]
        mcorr = np.where(np.isfinite(pct[sel]) & (np.nan_to_num(N) >= 20), pct[sel] * N / 100, np.nan)
        res = re_logit(mcorr, np.where(np.isfinite(mcorr), N, np.nan), X1)
        mm = np.nan_to_num(mcorr)
        v = np.where(np.isfinite(mcorr), 1 / (mm + 0.5) + 1 / (np.nan_to_num(N) - mm + 0.5), np.nan)
        with np.errstate(invalid='ignore'):
            vbar = np.nanmedian(v, 1)
            share = res['tau2'] / (res['tau2'] + vbar)
        share = share[np.isfinite(share)]
        if len(share) >= 20:
            vio.append(share)
            vcodes.append(code)
            tau_med[code] = float(np.nanmedian(res['tau2']))
            rel_rows.append((plain(code), len(share), float(np.median(share)), float(np.mean(share > 0.5)), tau_med[code]))
    if vio:
        parts = ax.violinplot(vio, positions=range(len(vio)), showextrema=False, widths=0.8)
        for body, code in zip(parts['bodies'], vcodes):
            body.set_facecolor(MOD_COLOUR[code])
            body.set_alpha(0.6)
            body.set_edgecolor('none')
        for i, v in enumerate(vio):
            ax.hlines(np.median(v), i - 0.25, i + 0.25, color='black', lw=1.2)
        mod_axis(ax, vcodes)
        ax.set_ylim(0, 1)
        ax.set_ylabel('between-sample share of variance')
        ax.set_title('τ² / (τ² + sampling variance) per site, logit scale: the reliability of one sample at its depth', fontsize=8)
    else:
        empty(ax, 'too few sites')
    letter(ax, 'F')
    tF = c.table('Fig14_reliability.tsv', ['mod', 'sites', 'median_share', 'fraction_share_gt_half', 'median_tau2'], rel_rows)
    c.tau_median = tau_med
    ax = fig.add_subplot(gs[1, 3])                                  # G power for a two-group design
    cov = np.logspace(np.log10(20), 3, 50)
    zz = stats.norm.isf(1e-3 / 2) + stats.norm.isf(0.2)
    for code in MAIN_MODS:
        if code not in tau_med:
            continue
        for npg, ls in ((5, ':'), (10, '--'), (20, '-')):
            se = np.sqrt((tau_med[code] + 1 / (cov * 0.21)) * 2 / npg)
            ax.plot(cov, 100 * (expit(logit(0.3) + zz * se) - 0.3), color=MOD_COLOUR[code], ls=ls, lw=1,
                    label=f'{MOD_NAME[code]}' if npg == 20 else None)
    ax.set_xscale('log')
    log_axis(ax.xaxis)
    ax.set_xlabel('reads per sample at the site')
    ax.set_ylabel('detectable Δ from 30% (points)')
    ax.set_title('5 (dotted), 10 (dashed), 20 (solid) per group', fontsize=8)
    ax.legend(fontsize=7)
    letter(ax, 'G')
    fig_batch_legend(fig, c.colours)
    c.legends['Fig14_measurement_vs_biology'] = f"""How much of the between-sample variation in occupancy is measurement, and how much is the donors? ({k} samples)
A, ICC(2,1) at half depth: between the two halves of each sample's reads (filled; a binomial split of N_valid and
a hypergeometric split of N_mod, so only read sampling differs: the reliability ceiling), and between donor pairs
each thinned to half depth (open: same batch; grey: different batches). Sites are those called in at least one
other sample, so the selection does not depend on the pair. Donor pairs below the split-half points: real donor
(or batch) differences; donor pairs at the split-half level: the donors do not differ beyond sampling noise.
B, per m6A site measured in >= {need} of {k} samples: observed SD of corrected occupancy across samples against
the SD binomial sampling alone gives at each sample's depth (the same test as figure 9A for isoforms);
dashed: equal; dotted: twice. C, % of sites at >= 2x per modification (number: median ratio).
D, rarefaction: sites passing 20 reads and 20% when each sample's reads are binomially thinned to 10-75%
(computed on every cov20 row; rows with fewer than 4 modified calls cannot pass and are skipped), against the
thinned depth. Curves still rising at full depth mean site counts measure depth (cf. Nelson et al. 2026: no
plateau at 28 M reads). E, the same relative to full depth, median over samples.
F, per site, the between-sample variance tau^2 (random-effects logit, empirical-Bayes REML, no covariates) as a
share of the median sampling variance: the reliability of a single sample's value at its depth; near 0 the site's
differences between samples are almost all counting noise. G, smallest detectable change from 30% (alpha 1e-3,
80% power) with 5, 10 or 20 samples per group, from the median tau^2 per modification.
Tables: tables/{tA}, tables/{tC}, tables/{tF}, tables/Fig14_rarefaction.tsv."""
    return fig


# ---------------------------------------------------------------- figure 15: Bayesian measurement model
PI_GRID = np.linspace(0, 1, 101)


def _loglik(m, n, fp, sens):
    th = np.clip(sens * PI_GRID[None, :] + fp[:, None] * (1 - PI_GRID[None, :]), 1e-12, 1 - 1e-12)
    return m[:, None] * np.log(th) + (n - m)[:, None] * np.log1p(-th)


def eb_prior(m, n, fp, sens, rng, max_rows=200000, iters=60):
    """Nonparametric prior on the occupancy grid (EM over site-samples; g-modelling)."""
    if len(m) > max_rows:
        i = rng.choice(len(m), max_rows, replace=False)
        m, n, fp = m[i], n[i], fp[i]
    L = _loglik(m, n, fp, sens)
    L -= L.max(1, keepdims=True)
    lik = np.exp(L)
    w = np.full(len(PI_GRID), 1 / len(PI_GRID))
    kern = np.array([0.25, 0.5, 0.25])
    for _ in range(iters):                                          # EMS: EM plus a light smoothing step
        post = lik * w
        post /= post.sum(1, keepdims=True)
        w = np.convolve(post.mean(0), kern, mode='same') + 1e-10
        w /= w.sum()
    return w


def posterior(m, n, fp, sens, w, tau, chunk=50000):
    """Posterior mean, 2.5% and 97.5% quantiles and P(pi >= tau) on the grid."""
    out = np.full((len(m), 4), np.nan)
    lw = np.log(w)
    cdf_idx = PI_GRID >= tau
    for s0 in range(0, len(m), chunk):
        sl = slice(s0, s0 + chunk)
        L = _loglik(m[sl], n[sl], fp[sl], sens) + lw
        L -= L.max(1, keepdims=True)
        p = np.exp(L)
        p /= p.sum(1, keepdims=True)
        cum = np.cumsum(p, 1)
        out[sl, 0] = p @ PI_GRID
        out[sl, 1] = PI_GRID[np.argmax(cum >= 0.025, 1)]
        out[sl, 2] = PI_GRID[np.argmax(cum >= 0.975, 1)]
        out[sl, 3] = p[:, cdf_idx].sum(1)
    return out


def fig15_bayes(c, notes):
    try:
        rows, ms, pct, mods, nv = c.mods_matrix()
    except FileNotFoundError:
        notes.append('figure 15 skipped: no cohort modification matrices')
        return None
    nvc, nmc, fpc, clipped = c.mods_counts()
    if nvc is None:
        notes.append('figure 15 skipped: needs the n_valid and n_mod matrices')
        return None
    a = c.args
    rng = np.random.default_rng(a.seed)
    sens, tau, pcall = a.sensitivity, a.posterior_min, a.posterior_prob
    k = len(ms)
    priors, post_all, calls = {}, {}, {}
    long_rows = []
    for code in ALL_MODS:
        sel = np.where(mods == code)[0]
        if len(sel) < 20:
            continue
        N, M, F, P = nvc[sel], nmc[sel], fpc[sel], pct[sel]
        meas = np.isfinite(N) & np.isfinite(M) & (np.nan_to_num(N) >= 20) & np.isfinite(P)
        si, sj = np.where(meas)
        m_, n_, f_ = M[si, sj], N[si, sj], np.nan_to_num(F[si, sj])
        w = eb_prior(m_, n_, f_, sens, rng)
        po = posterior(m_, n_, f_, sens, w, tau)
        priors[code] = w
        post_all[code] = (sel, si, sj, po, P[si, sj], n_)
        rule20 = P[si, sj] >= 20
        rulep = po[:, 3] >= pcall
        calls[code] = [(int(rule20[sj == j].sum()), int(rulep[sj == j].sum()), int((rule20 & rulep)[sj == j].sum())) for j in range(k)]
        kk = np.where(rule20 | rulep)[0]
        r_ = sel[si[kk]]
        long_rows.append(pd.DataFrame({
            'chrom': rows['chrom'][r_], 'start': rows['start'][r_], 'strand': rows['strand'][r_], 'mod': plain(code),
            'sample': np.array(ms, dtype=object)[sj[kk]], 'n_valid': n_[kk].astype(np.int64), 'n_mod': m_[kk].astype(np.int64),
            'gIVT_fp': f_[kk], 'corrected_pct': P[si[kk], sj[kk]], 'posterior_mean_pct': 100 * po[kk, 0],
            'cri_lo_pct': 100 * po[kk, 1], 'cri_hi_pct': 100 * po[kk, 2], f'p_pi_ge_{int(100 * tau)}': po[kk, 3]}))
    if not priors:
        notes.append('figure 15 skipped: too few measured sites')
        return None
    fig, gs = figure(2, 4, widths=[PANEL[0], PANEL[0] * 1.15, PANEL[0] * 1.2, PANEL[0]])
    ax = fig.add_subplot(gs[0, 0])                                  # A priors
    for code, w in priors.items():
        ax.plot(100 * PI_GRID, w, color=MOD_COLOUR[code], lw=1.1, ls='-' if code in MAIN_MODS else '--', label=MOD_NAME[code])
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel('occupancy π (%)')
    ax.set_ylabel('prior mass (grid of 1%)')
    ax.legend(fontsize=6.5, ncol=2)
    ax.set_title('empirical-Bayes prior per modification', fontsize=9)
    letter(ax, 'A')
    ax = fig.add_subplot(gs[0, 1])                                  # B shrinkage, m6A
    if 'a' in post_all:
        sel, si, sj, po, naive, n_ = post_all['a']
        idx = rng.choice(len(po), min(20000, len(po)), replace=False)
        sc = ax.scatter(naive[idx], 100 * po[idx, 0], c=np.log10(n_[idx]), s=2, cmap='viridis', lw=0, rasterized=True)
        ax.plot([0, 100], [0, 100], color='black', ls='--', lw=0.8)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_xlabel('gIVT-corrected occupancy (%)')
        ax.set_ylabel('posterior mean (%)')
        cb = colourbar(fig, sc, ax, label='log$_{10}$ reads')
        ax.set_title('m$^6$A: shallow sites shrink toward the prior', fontsize=9)
    letter(ax, 'B')
    ax = fig.add_subplot(gs[0, 2])                                  # C calls per sample, two rules
    codes = [m for m in ALL_MODS if m in calls]
    for i, code in enumerate(codes):
        paired_points(ax, i, [x[0] for x in calls[code]], [x[1] for x in calls[code]], MOD_COLOUR[code])
    ax.set_yscale('log')
    log_axis(ax.yaxis)
    mod_axis(ax, codes, fontsize=8)
    ax.set_xlim(-0.6, len(codes) - 0.4)
    ax.set_ylabel('sites per sample')
    ax.set_title(f'open: ≥ 20 reads, ≥ 20%; filled: P(π ≥ {int(100 * tau)}%) ≥ {pcall:g}', fontsize=8)
    letter(ax, 'C')
    ax = fig.add_subplot(gs[0, 3])                                  # D against depth
    depth = c.depth() if os.path.exists(os.path.join(c.cohort, f'cohort.{c.ref}.read_qc.tsv')) else {}
    if 'a' in calls and depth:
        x = np.array([depth.get(s, np.nan) / 1e6 for s in ms])
        for j, (col, lab, mk) in enumerate((('0.55', '20% rule', 'o'), ('#398DC3', 'posterior', 's'))):
            y = np.array([v[j] for v in calls['a']], float)
            ax.scatter(x, y, color=col, marker=mk, s=20, label=lab, zorder=3)
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= 3:
                ax.plot(np.sort(x[ok]), np.poly1d(np.polyfit(x[ok], y[ok], 1))(np.sort(x[ok])), color=col, lw=0.8, ls='--')
        ax.set_xlabel('non-globin reads (M)')
        ax.set_ylabel('m$^6$A sites')
        ax.legend(fontsize=7)
    else:
        empty(ax, 'needs read_qc and m6A sites')
    letter(ax, 'D')
    ax = fig.add_subplot(gs[1, 0:2])                                # E forest plot, most variable m6A sites
    if 'a' in post_all:
        sel, si, sj, po, naive, n_ = post_all['a']
        counts = np.bincount(si, minlength=len(sel))
        var = np.zeros(len(sel))
        for r_ in np.where(counts >= max(3, math.ceil(0.75 * k)))[0]:
            var[r_] = np.var(po[si == r_, 0])
        top = np.argsort(var)[::-1][:8]
        top = [t for t in top if var[t] > 0]
        for yi, r_ in enumerate(top):
            for i_ in np.where(si == r_)[0]:
                s = ms[sj[i_]]
                y = yi + (sj[i_] - k / 2) / (k * 1.4)
                ax.plot([100 * po[i_, 1], 100 * po[i_, 2]], [y, y], color=c.colours[c.batch[s]], lw=0.7)
                mk = c.gmarker.get(c.group.get(s), 'o') if c.group else 'o'
                ax.scatter([100 * po[i_, 0]], [y], color=c.colours[c.batch[s]], s=10, marker=mk, zorder=3)
        ax.set_yticks(range(len(top)), [f'{rows["chrom"][sel[r_]]}:{rows["start"][sel[r_]] + 1}' for r_ in top], fontsize=7)
        ax.set_ylim(len(top) - 0.5, -0.5)
        ax.set_xlim(0, 100)
        ax.set_xlabel('m$^6$A occupancy, posterior mean and 95% interval (%)')
        ax.set_title('the 8 m$^6$A sites whose posterior means vary most between samples', fontsize=9)
    letter(ax, 'E')
    ax = fig.add_subplot(gs[1, 2])                                  # F agreement of the rules
    for i, code in enumerate(codes):
        r20 = sum(x[0] for x in calls[code])
        rp = sum(x[1] for x in calls[code])
        both = sum(x[2] for x in calls[code])
        jac = both / max(1, r20 + rp - both)
        ax.bar(i, jac, color=MOD_COLOUR[code], width=0.7)
        ax.text(i, jac + 0.02, f'{r20 - both:,}\n{rp - both:,}', ha='center', fontsize=6)
    mod_axis(ax, codes, fontsize=8)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel('Jaccard (20% rule, posterior)')
    ax.set_title('numbers: only-20%, only-posterior site-samples', fontsize=8)
    letter(ax, 'F')
    ax = fig.add_subplot(gs[1, 3])                                  # G interval width by depth
    be = np.array([20, 30, 50, 80, 130, 200, 320, 500, 800, 1300, 1e9])
    for code in MAIN_MODS:
        if code not in post_all:
            continue
        _, _, _, po, naive, n_ = post_all[code]
        width = 100 * (po[:, 2] - po[:, 1])
        b_ = np.searchsorted(be, n_, side='right') - 1
        med = [np.median(width[(b_ == i) & (naive >= 20)]) if ((b_ == i) & (naive >= 20)).sum() >= 10 else np.nan for i in range(len(be) - 1)]
        ax.plot(be[:-1], med, color=MOD_COLOUR[code], marker='o', ms=3, lw=1.1, label=MOD_NAME[code])
    ax.set_xscale('log')
    ax.set_xlabel('reads at the site')
    ax.set_ylabel('95% interval width (points)')
    ax.set_title('sites ≥ 20% corrected', fontsize=8)
    ax.legend(fontsize=7)
    letter(ax, 'G')
    fig_batch_legend(fig, c.colours)
    pd.concat(long_rows, ignore_index=True).to_csv(os.path.join(c.tables, 'Fig15_posterior_calls.tsv.gz'), sep='\t',
                                                   index=False, float_format='%.4g', compression='gzip')
    c.table('Fig15_calls_per_sample.tsv', ['sample', 'mod', 'calls_20pct', 'calls_posterior', 'calls_both'],
            [(ms[j], plain(code), *calls[code][j]) for code in calls for j in range(k)])
    for code in ('a',):
        if code in calls:
            for j, s in enumerate(ms):
                c.note_sample(s, 'posterior_calls_m6A', calls[code][j][1])
    nclip = int(np.nansum(clipped)) if clipped is not None else 0
    c.legends['Fig15_bayesian_occupancy'] = f"""A measurement model for occupancy, instead of subtracting the gIVT rate and cutting at 20% ({k} samples).
For each site and sample: N_mod ~ Binomial(N_valid, theta), theta = s pi + fp (1 - pi): a modified molecule is
called modified with sensitivity s = {sens:g} (--sensitivity; the paper's synthetic oligos gave recall 0.90-0.98) and
an unmodified one with the gIVT false-positive rate fp of its 9-mer; pi is the occupancy. The prior on pi is
estimated per modification from all of its site-samples (nonparametric, on a 1% grid, by smoothed EM (EMS): empirical Bayes,
i.e. partial pooling across sites; Statistical Rethinking ch. 15 builds the same measurement-error idea). Sites are
the cohort's final-stage sites (called >= 20% in at least one sample), so the prior describes those sites, not all
positions. {nclip:,} site-samples where gIVT clipped the corrected value to 0 use fp as a lower bound unless the
stream pass recovered it (--errortable).
A, the priors. B, m6A: posterior mean against the corrected value, coloured by depth: shallow, extreme values move
toward the prior, deep ones stay. C, per sample: sites passing the 20% rule (open) and with P(pi >= {tau:g}) >= {pcall:g}
(filled). D, both counts against depth. E, the eight m6A sites (measured in >= 75% of samples) whose posterior
means vary most between samples: each sample's posterior mean and 95% interval, coloured by batch. F, agreement of
the two call rules (Jaccard over site-samples; numbers: calls made by one rule only). G, 95% interval width against
depth for sites >= 20%: what a single sample's value is worth. Every call by either rule, with its posterior:
tables/Fig15_posterior_calls.tsv.gz."""
    return fig


# ---------------------------------------------------------------- Table 1, legends, main
def table1(c):
    """Per-sample summary (the paper's Table 1 role): depth, composition, sites, tails, QC and calibration."""
    sites = {}
    try:
        _, ms, pct, mods, _ = c.mods_matrix()
        for code in ALL_MODS:
            called = (np.nan_to_num(pct[mods == code]) >= 20).sum(0)
            for s, n in zip(ms, called):
                sites.setdefault(s, {})[f'sites_{plain(code)}'] = int(n)
    except FileNotFoundError:
        pass
    comp = c.composition()
    thr = c.thresholds()
    cols = (['mapped_reads', 'globin_pct', 'mito_pct', 'non_globin_non_mito_reads', 'genes_ge10_reads']
            + [f'sites_{plain(m)}' for m in ALL_MODS] + ['sites_all', 'posterior_calls_m6A']
            + ['polya_median_per_read_nuclear', 'polya_median_of_gene_medians_nuclear', 'polya_median_per_read_mito',
               'pct_genes_major_isoform_not_dominant', 'identity_median', 'complete_5p', 'read_len_median',
               'speed_median', 'qscore_median', 'fail_frac_m6A', 'fail_frac_Psi', 'fail_frac_m5C', 'fail_frac_Ino',
               'threshold_m6A', 'null_lambda_m6A'] + [f'composition_{n}' for n in comp['names']])
    vcols = [k for k in ('basecall_model', 'mod_models', 'dorado_version') if k in c.sheet.columns]
    idx = dict(zip(c.sheet['sample_id'], range(len(c.sheet))))
    rows = []
    for i, s in enumerate(c.samples):
        d = dict(c.summary.get(s, {}))
        d.update(sites.get(s, {}))
        if s in sites:
            d['sites_all'] = sum(sites[s].values())
        d['threshold_m6A'] = thr.get(s, {}).get('a', float('nan'))
        for j, n in enumerate(comp['names']):
            d[f'composition_{n}'] = float(comp['values'][i, j])
        rows.append([s, c.batch[s]] + [d.get(k, float('nan')) for k in cols] + [c.sheet[k].iloc[idx[s]] if s in idx else '' for k in vcols])
    return c.table('Table1_cohort_summary.tsv', ['sample', 'batch'] + cols + vcols, rows)


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


def existing(path):
    return path if path and os.path.exists(path) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True)
    ap.add_argument('--samplesheet', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ref', default='GRCh38')
    ap.add_argument('--batch-column', default='batch')
    ap.add_argument('--read-stats', help='directory of read_stats.sh outputs (default: <results>/read_stats)')
    ap.add_argument('--sample-qc', help='cohort_extract.py outputs (default: <results>/sample_qc)')
    ap.add_argument('--mask-contigs', default=SITE_FILTER_DEFAULT,
                    help=f'regex of contigs to drop, POSIX ERE (awk reads it too; default: {SITE_FILTER_DEFAULT}); "" for none')
    ap.add_argument('--mask-bed', help='BED of regions to drop from the modification figures')
    ap.add_argument('--group', help='samplesheet column to compare (two levels): adds figure 8')
    ap.add_argument('--levels', nargs=2, metavar=('A', 'B'))
    ap.add_argument('--fdr', type=float, default=0.1)
    ap.add_argument('--covariates', nargs='+', help='samplesheet columns for figure 10 (default: every non-id, non-path column)')
    ap.add_argument('--model-covariates', nargs='*',
                    help='covariates in the figure 8 models (default: the batch column and composition PC1; none: give the flag alone)')
    ap.add_argument('--composition', help='TSV of measured composition (sample + one column per cell type; CBC or CIBERSORTx)')
    ap.add_argument('--composition-columns', nargs='+', help='samplesheet columns holding measured composition')
    ap.add_argument('--reference', default=existing(CLUSTER_REFERENCE), help='uncompressed FASTA with .fai (DRACH; exact gIVT rates)')
    ap.add_argument('--gtf', default=existing(CLUSTER_GTF), help='annotation GTF for the metagene (GENCODE: Ensembl_canonical)')
    ap.add_argument('--errortable', default=existing(CLUSTER_ERRORTABLE), help="gIVT 9-mer table (the pipeline's givt_error_table)")
    ap.add_argument('--known-sites', nargs='+', metavar='MOD=BED', help='orthogonal site sets (figure 16; orthogonal_to_bed.py writes them), e.g. m6A=GLORI2022.m6A.bed')
    ap.add_argument('--repeats', help='RepeatMasker BED (name in column 4) for the inosine-in-Alu panel')
    ap.add_argument('--icc-selection', choices=['others', 'either'], default='others',
                    help='sites for pairwise ICC: called in another sample (default) or in either of the pair (v4)')
    ap.add_argument('--sensitivity', type=float, default=0.95, help='figure 15: P(called modified | modified)')
    ap.add_argument('--posterior-min', type=float, default=0.10, help='figure 15: call if P(occupancy >= this) >= --posterior-prob')
    ap.add_argument('--posterior-prob', type=float, default=0.95)
    ap.add_argument('--compare-batch', help='adds figure S1: this batch against all others')
    ap.add_argument('--only', nargs='+', choices=FIGURES, help='make only these figures')
    ap.add_argument('--no-stream', action='store_true', help='skip the cov20_fpcorr stream pass (figures 11-15 lose panels)')
    ap.add_argument('--chunk-rows', type=int, default=2_000_000, help='rows per chunk in the stream pass')
    ap.add_argument('--threads', type=int, default=min(8, os.cpu_count() or 1), help='files processed at once (figure 3, stream pass)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--style', default=LAB_STYLE, help=f'matplotlib style sheet (default: {LAB_STYLE}); "" for none')
    ap.add_argument('--font-dir', nargs='*', default=[], help='directories of .ttf/.otf fonts to register (Helvetica)')
    args = ap.parse_args()
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s %(message)s', datefmt='%H:%M:%S'))
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    log.propagate = False
    os.makedirs(args.out, exist_ok=True)
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    notes = apply_style(args.style, args.font_dir)
    mask = SiteMask(args.mask_contigs, args.mask_bed)
    c = Cohort(args, mask)
    rs_dir = args.read_stats or c.path('read_stats')
    qc = load_sample_qc(c, args.sample_qc or c.path('sample_qc'), notes)
    want = set(args.only or FIGURES)
    if not args.no_stream and want & {'Fig5', 'Fig8', 'Fig10', 'Fig11', 'Fig12', 'Fig13', 'Fig14', 'Fig15'}:
        stream_pass(c, args.threads, notes)
    figures = [('Fig1', 'Fig1_read_qc', lambda: fig1_reads(c, rs_dir, notes)),
               ('Fig2', 'Fig2_cohort_overview', lambda: fig2_overview(c, notes)),
               ('Fig3', 'Fig3_ivt_correction', lambda: fig3_ivt(c, args.out, args.threads, notes)),
               ('Fig4', 'Fig4_modification_calls', lambda: fig4_mods(c, notes)),
               ('Fig5', 'Fig5_between_donor_agreement', lambda: fig5_concordance(c, notes)),
               ('Fig6', 'Fig6_isoform_modifications', lambda: fig6_isoforms(c, notes)),
               ('Fig7', 'Fig7_polya', lambda: fig7_polya(c, notes))]
    if args.group:
        figures.append(('Fig8', 'Fig8_group_comparison', lambda: fig8_groups(c, args.fdr, notes)))
    figures += [('Fig9', 'Fig9_isoform_usage', lambda: fig9_usage(c, args.fdr, notes)),
                ('Fig10', 'Fig10_cohort_drivers', lambda: fig10_drivers(c, qc, notes)),
                ('Fig11', 'Fig11_integrity_signal', lambda: fig11_integrity(c, qc, notes)),
                ('Fig12', 'Fig12_ivt_calibration', lambda: fig12_calibration(c, notes)),
                ('Fig13', 'Fig13_controls_plausibility', lambda: fig13_controls(c, notes)),
                ('Fig14', 'Fig14_measurement_vs_biology', lambda: fig14_measurement(c, notes)),
                ('Fig15', 'Fig15_bayesian_occupancy', lambda: fig15_bayes(c, notes)),
                ('Fig16', 'Fig16_orthogonal_sites', lambda: fig16_orthogonal(c, notes))]
    if args.compare_batch:
        figures.append(('FigS1', 'FigS1_batch_comparison', lambda: figS1_batch(c, args.compare_batch, notes)))
    figures.append(('FigS2', 'FigS2_jaccard', lambda: figS2_jaccard(c, notes)))
    made, failed = [], []
    for key, name, make in figures:
        if args.only and key not in args.only:
            continue
        log.info('%s ...', name)
        try:
            fig = make()
        except FileNotFoundError as e:
            notes.append(f'{name} skipped: missing input {e.filename or e}')
            continue
        except SystemExit:
            raise
        except Exception as e:                                      # one figure's failure must not cost the others
            import traceback
            traceback.print_exc()
            failed.append(name)
            notes.append(f'{name} FAILED: {type(e).__name__}: {e}')
            plt.close('all')
            continue
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
        parts = [f'{plain(m)} {mask.removed[(rule, m)]:,}' for m in ALL_MODS if (rule, m) in mask.removed]
        if parts:
            notes.append(f'masked by {label}: cohort sites ' + ', '.join(parts))
    with open(os.path.join(args.out, 'notes.txt'), 'w') as fh:
        fh.write('\n'.join(notes) + '\n')
    for n in notes:
        print(f'note: {n}', file=sys.stderr)
    print(f'cohort_plots.py: {len(made)} figures in {args.out}: {" ".join(f"{m}.pdf" for m in made)}'
          + (f'; FAILED: {" ".join(failed)}' if failed else ''))
    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
