#!/usr/bin/env python3
"""Orthogonal RNA-modification site sets (supplementary xlsx / xlsb, GEO tsv / bed) -> BED6 for
cohort_plots.py --known-sites: chrom, 0-based start, end, set name, occupancy (0-100, '.' if none), strand.

    orthogonal_to_bed.py batch --dir /projects/Genometechlab/andrew/orthogonal --out <dir> [--inspect]
    orthogonal_to_bed.py one FILE --mod m6A --name GLORI2022 --out <dir> [--sheet REGEX] [--chrom-col C]
                         [--pos-col C] [--strand-col C] [--occ-col C ...] [--base 0|1] [--shift -1|0|1]
                         [--exclude FILE] [--inspect]

Columns are found from the header (the first row of the first 30 that names a chromosome and a position
column; title rows above it are skipped), or a site-id column such as chr1:12345:+. Occupancy is a column
named like ratio / level / stoichiometry / fraction (replicate columns averaged, a mean column preferred;
p-value, depth, count and knockout / mutant columns never used); override any choice with --*-col.

Coordinates are settled against the reference, not assumed: every combination of 0- or 1-based positions
and a 1-nt shift (in RNA 5'->3' direction) is scored by how often the reference base, on the site's strand,
is the modification's base (A for m6A and inosine, T for pseudouridine, C for m5C), and the best is used.
That also places sites that a method reports on a neighbouring base (DART-seq's C-to-U edit 3' of the m6A,
miCLIP crosslink mutations and truncations). A best match below --min-match (0.85) is reported as a likely
build or transcript-coordinate problem (lift it over first) and the set is not written unless --force.
Missing strands are inferred from the base; strand words (Forward/Reverse, plus/minus) are understood. A column
of transcript ids where chromosomes should be is ignored in favour of a genomic site id (PRAISE chr_site,
chr14_73219165-73219168: every position of a range, descending = minus strand). NJU-seq-style 2'-O-Me tables
(--mod Nm) are checked against their own modified-base column (--base-col) and split into Am, Cm, Gm and Um.
Replicate sheets can be combined, keeping sites in >= --min-sheets of them (SLIC-seq: 2 of 3 HEK293T replicates).
Occupancy never comes from control / input / untreated / knockout columns or from integer BED scores (counts).
--exclude drops sites also in a control set (DART-seq APOBEC1-YTHmut).

batch applies the presets below to the files in --dir, writes one BED per set, report.tsv, and prints the
--known-sites arguments for cohort_plots.py. .xlsb needs pyxlsb (pip install --user pyxlsb).
"""
import argparse
import glob
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

CLUSTER_REFERENCE = '/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
MOD_BASE = {'m6A': 'A', 'Ino': 'A', 'Psi': 'T', 'm5C': 'C', 'Nm': 'N', 'Am': 'A', 'Cm': 'C', 'Gm': 'G', 'Um': 'T'}
NM_SPLIT = {'A': 'Am', 'C': 'Cm', 'G': 'Gm', 'T': 'Um'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# filename regex -> settings; first match wins. sheet: regex for the sheet to use (largest match if several).
PRESETS = [
    (r'^GLORI2022_.*\.xlsx$', dict(mod='m6A', name='GLORI2022', sheet=r'(?i)hek|293')),
    # Sun et al. 2025 (Nat Methods): the 10 ng sheets are the paper's GLORI2.0 (101,613) and GLORI1.0 (76,452) sets
    (r'^GLORI2025_MOESM3.*\.xlsb$', dict(mod='m6A', name='GLORI2_10ng', sheet=r'^10ng')),
    (r'^GLORI2025_MOESM5.*\.xlsb$', dict(mod='m6A', name='GLORI1_10ng', sheet=r'^10ng')),
    (r'^m6ASACseq_.*\.bed(\.gz)?$', dict(mod='m6A', name='m6ASACseq')),
    (r'^DARTseq_.*APOBEC1YTH_processed\.hg38\.bed(\.gz)?$',
     dict(mod='m6A', name='DARTseq', exclude=r'^DARTseq_.*APOBEC1YTHmut_processed\.hg38\.bed(\.gz)?$')),
    (r'^miCLIP_.*abcam.*CIMS.*\.hg38\.bed(\.gz)?$', dict(mod='m6A', name='miCLIP_CIMS')),
    (r'^miCLIP_.*sysy.*CITS.*\.hg38\.bed(\.gz)?$', dict(mod='m6A', name='miCLIP_CITS')),
    (r'^BIDseq_.*\.xlsx$', dict(mod='Psi', name='BIDseq', sheet=r'(?i)hek|293|wt')),
    (r'^PRAISE_.*\.xlsx$', dict(mod='Psi', name='PRAISE', sheet=r'(?i)dataset 2$')),
    (r'^UBSseq_.*HEK293T-WT.*\.tsv(\.gz)?$', dict(mod='m5C', name='UBSseq')),
    (r'^SLICseq_.*\.xlsx$', dict(mod='Ino', name='SLICseq', sheet=r'^HEK293T-rep', min_sheets=2)),
    (r'^NJUseq_.*\.xlsx$', dict(mod='Nm', name='NJUseq', sheet=r'(?i)hek|293', base_col='Nm')),
]
SKIP_NOTES = [
    (r'^DARTseq_.*YTHmut', 'APOBEC1-YTHmut control: used to exclude background edits from DARTseq'),
    (r'^(DARTseq|miCLIP)_(?!.*\.hg38\.).*\.bed(\.txt)?(\.gz)?$', 'hg19 original: the .hg38 version is used'),
    (r'^UBSseq_.*HeLa', 'HeLa, not HEK293T: convert with "one" if wanted (shNSUN2 vs shControl defines NSUN2 sites)'),
    (r'(MD5SUMS|\.chain(\.gz)?|\.log|liftover_report\.tsv|inspect\.txt)$', 'not a site table'),
]

CHROM = {'chr', 'chrom', 'chromosome', 'seqnames', 'seqname', 'chromname', 'chrname', 'contig', 'seqid', 'chromosomename'}
POS = {'pos', 'position', 'site', 'sites', 'coordinate', 'genomicposition', 'gpos', 'loc', 'location', 'start', 'chromstart',
       'pos1', 'sitepos', 'siteposition', 'genomepos', 'genomicsite', 'modsite', 'modificationsite', 'genomicloci', 'locus'}
END = {'end', 'chromend', 'stop', 'pos2'}
STRAND = {'strand', 'str', 'genestrand', 'sitestrand'}
OCC_WORDS = ('ratio', 'level', 'stoich', 'fraction', 'occupancy', 'editrate', 'editingrate', 'editinglevel', 'methylation',
             'modrate', 'modlevel', 'percent', 'frac', 'rate')
OCC_NEVER = ('pval', 'pvalue', 'fdr', 'qval', 'padj', 'depth', 'coverage', 'count', 'reads', 'number', 'total', 'cov',
             'score', 'log', 'fold')
CONDITION_NEVER = re.compile(r'(?i)(^|[_\W])(ko|kd|knock\w*|mut\w*|inhib\w*|stm\w*|treat\w*|sh\w+|si\w+|dko|del|'
                             r'ctrl|control|input|mock|untreated|ivt|neg\w*)([_\W]|$)')
SITE_ID = re.compile(r'^(chr[0-9A-Za-z_]+?|[0-9]{1,2}|[XYM]T?)[:_|-](\d+)(?:-(\d+))?(?:[:_|]([+-]))?\s*$')
CHROM_LIKE = re.compile(r'^(chr)?([0-9]{1,2}|X|Y|M|MT|Un)([_.]\w+)?$')
PLUS = {'+', 'plus', 'forward', 'fwd', 'f', '1', '+1', 'sense', 'pos', 'positive', 'w', 'watson'}
MINUS = {'-', 'minus', 'reverse', 'rev', 'r', '-1', 'antisense', 'neg', 'negative', 'c', 'crick'}


def norm(x):
    return re.sub(r'[^a-z0-9]', '', str(x).lower())


# ---------------------------------------------------------------- reference
class Ref:
    def __init__(self, path):
        if not path or not os.path.exists(path) or not os.path.exists(path + '.fai'):
            sys.exit(f'orthogonal_to_bed.py: needs an uncompressed reference FASTA with .fai (--reference); got {path}')
        self.fai = {}
        with open(path + '.fai') as fh:
            for line in fh:
                f = line.split('\t')
                self.fai[f[0]] = tuple(int(x) for x in f[1:5])
        self.mm = np.memmap(path, dtype=np.uint8, mode='r')

    def bases(self, chroms, pos0):
        """Upper-case reference base (forward strand) at 0-based positions; 'N' off contig."""
        out = np.full(len(pos0), ord('N'), np.uint8)
        chroms = np.asarray(chroms, dtype=object)
        pos0 = np.asarray(pos0, np.int64)
        for c in set(chroms.tolist()):
            if c not in self.fai:
                continue
            L, off, lb, lw = self.fai[c]
            sel = np.where(chroms == c)[0]
            p = pos0[sel]
            ok = (p >= 0) & (p < L)
            out[sel[ok]] = self.mm[off + (p[ok] // lb) * lw + p[ok] % lb]
        low = (out >= 97) & (out <= 122)
        out[low] -= 32
        return np.array([chr(b) for b in out], dtype=object)

    def fix_chrom(self, chroms):
        """Match the reference's contig naming (1 -> chr1, MT -> chrM)."""
        names = set(self.fai)
        out = []
        for c in chroms:
            c = str(c).strip()
            if c in names:
                out.append(c)
            elif f'chr{c}' in names:
                out.append(f'chr{c}')
            elif c in ('MT', 'chrMT', 'M') and 'chrM' in names:
                out.append('chrM')
            elif c.startswith('chr') and c[3:] in names:
                out.append(c[3:])
            else:
                out.append(c)
        return np.array(out, dtype=object)


# ---------------------------------------------------------------- reading
def read_sheets(path):
    """{sheet: DataFrame with header=None, all text}."""
    low = path.lower()
    if low.endswith(('.xlsx', '.xlsm', '.xls')):
        return pd.read_excel(path, sheet_name=None, header=None, dtype=str)
    if low.endswith('.xlsb'):
        try:
            return pd.read_excel(path, sheet_name=None, header=None, dtype=str, engine='pyxlsb')
        except ImportError:
            sys.exit(f'orthogonal_to_bed.py: {os.path.basename(path)} is .xlsb: pip install --user pyxlsb (or save it as .xlsx)')
    sep = ',' if re.search(r'\.csv(\.gz)?$', low) else '\t'
    op = gzip.open if low.endswith('.gz') else open
    with op(path, 'rt') as fh:
        lines = [ln.rstrip('\n') for ln in fh if not ln.startswith(('track', 'browser'))]
    rows = [ln.split(sep) for ln in lines if ln.strip()]
    width = max(len(r) for r in rows) if rows else 0
    return {'-': pd.DataFrame([r + [None] * (width - len(r)) for r in rows], dtype=str)}


def looks_bed(df):
    """No header and BED-shaped: chrom-like, integer start, integer end >= start."""
    if df.shape[1] < 3 or not len(df):
        return False
    r = df.iloc[0]
    try:
        a, b = int(r[1]), int(r[2])
    except (TypeError, ValueError):
        return False
    return b >= a and not str(r[0]).lstrip('#').lower().startswith(('chrom', 'chr\t'))


def find_header(df):
    best, bi = 0, None
    for i in range(min(30, len(df))):
        cells = [norm(x) for x in df.iloc[i].tolist() if x is not None and str(x) != 'nan']
        sc = sum(1 for x in cells if x in CHROM | POS | END | STRAND) + sum(1 for x in cells if any(w in x for w in OCC_WORDS))
        nxt = df.iloc[i + 1].tolist() if i + 1 < len(df) else []
        has_core = any(x in CHROM for x in cells) or any(SITE_ID.match(str(y)) for y in nxt if y is not None)
        if sc > best and (has_core or sc >= 3):
            best, bi = sc, i
    return bi


def numeric(s):
    return pd.to_numeric(pd.Series(s).astype(str).str.replace('%', '', regex=False).str.strip(), errors='coerce')


def detect(df, over):
    """Column mapping from a header row (or BED layout), with overrides. Returns (data, mapping, notes)."""
    notes = []
    if looks_bed(df) and not any(over.get(k) for k in ('chrom', 'pos')):
        data = df.copy()
        n = data.shape[1]
        data.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand'][:n] + [f'c{i}' for i in range(6, n)]
        occ = []
        if n >= 5:
            v = numeric(data['score'])
            frac = ((v % 1) != 0).mean() if v.notna().any() else 0
            if v.notna().mean() > 0.9 and v.max() <= 100 and v.min() >= 0 and v.nunique() > 2 and frac > 0.5:
                occ = ['score']
        m = dict(chrom='chrom', pos='start', end='end', strand='strand' if n >= 6 else None, occ=occ, siteid=None, base=0)
        return data, m, ['BED layout (0-based starts)']
    hi = find_header(df)
    if hi is None:
        hi = 0
        notes.append('no recognisable header: first row used')
    names = [str(x) if x is not None else f'col{i}' for i, x in enumerate(df.iloc[hi].tolist())]
    seen = {}
    for i, n_ in enumerate(names):                                   # unique names
        if n_ in seen:
            seen[n_] += 1
            names[i] = f'{n_}.{seen[n_]}'
        else:
            seen[n_] = 0
    data = df.iloc[hi + 1:].copy()
    data.columns = names
    data = data.dropna(how='all')
    if hi:
        notes.append(f'header on row {hi + 1} ({hi} title row(s) skipped)')
    nm = {c: norm(c) for c in names}
    pick = lambda pool: next((c for c in names if nm[c] in pool), None)   # noqa: E731
    m = dict(chrom=over.get('chrom') or pick(CHROM), pos=over.get('pos') or pick(POS - {'start'}) or pick({'start'}),
             end=over.get('end') or pick(END), strand=over.get('strand') or pick(STRAND), siteid=None, base=1)
    if m['pos'] and nm.get(m['pos']) in ('start', 'chromstart') and m['end']:
        e = numeric(data[m['end']]) - numeric(data[m['pos']])
        if (e == 1).mean() > 0.9:
            m['base'] = 0
            notes.append('start/end one apart: BED-style 0-based starts')
    if m['chrom'] and not over.get('chrom'):
        v = data[m['chrom']].dropna().astype(str).head(500)
        if len(v) and v.str.match(CHROM_LIKE).mean() < 0.5:
            notes.append(f'"{m["chrom"]}" holds {v.iloc[0]!r}-like values, not chromosomes: ignored')
            m['chrom'], m['pos'] = None, None
    if not m['chrom'] or not m['pos']:
        for c in names:
            v = data[c].dropna().astype(str).head(200)
            if len(v) and v.str.match(SITE_ID).mean() > 0.8:
                m['siteid'] = c
                notes.append(f'positions from site ids in "{c}"')
                break
    if over.get('occ'):
        occ = list(over['occ'])
    else:
        cand = []
        for c in names:
            k = nm[c]
            if not any(w in k for w in OCC_WORDS) or any(w in k for w in OCC_NEVER) or CONDITION_NEVER.search(c):
                continue
            v = numeric(data[c])
            if v.notna().mean() > 0.8 and v.min() >= 0 and v.max() <= 100 and v.nunique() > 2:
                cand.append(c)
        pref = [c for c in cand if re.search(r'(?i)mean|avg|ave(rage)?\b|merged|combined|pooled', c)]
        occ = pref[:1] or cand
        if len(occ) > 1:
            notes.append(f'occupancy = mean of {", ".join(occ)}')
    m['occ'] = occ
    return data, m, notes


def strand_of(values):
    out = []
    for x in values:
        k = str(x).strip().lower()
        out.append('+' if k in PLUS else '-' if k in MINUS else None)
    return np.array(out, dtype=object)


def sites_from(data, m, base_col=None):
    """(chrom, raw position, strand or None, occupancy or NaN, expected base or None, dropped rows).
    A site-id range (chr14_73219165-73219168) becomes every position in it; descending ranges are minus strand."""
    occ = np.full(len(data), np.nan)
    if m['occ']:
        with np.errstate(all='ignore'), __import__('warnings').catch_warnings():
            __import__('warnings').simplefilter('ignore')
            occ = np.nanmean(np.column_stack([numeric(data[c]).to_numpy(float) for c in m['occ']]), 1)
    want = None
    if base_col and base_col in data.columns:
        want = data[base_col].astype(str).str.strip().str.upper().str.replace('U', 'T').to_numpy(dtype=object)
    if m['siteid']:
        mt = data[m['siteid']].astype(str).str.extract(SITE_ID.pattern)
        a, b = numeric(mt[1]).to_numpy(float), numeric(mt[2]).to_numpy(float)
        st = strand_of(data[m['strand']]) if m['strand'] else strand_of(mt[3])
        rng_strand = np.where(np.isfinite(b) & (b < a), '-', np.where(np.isfinite(b) & (b > a), '+', None))
        st = np.where(st == None, rng_strand, st)                     # noqa: E711
        ch = mt[0].to_numpy(dtype=object)
        rows = []
        for i in range(len(data)):
            if ch[i] is None or ch[i] != ch[i] or not np.isfinite(a[i]):
                continue
            lo, hi = (a[i], b[i]) if np.isfinite(b[i]) else (a[i], a[i])
            lo, hi = int(min(lo, hi)), int(max(lo, hi))
            if hi - lo > 50:                                        # not a site
                continue
            for p in range(lo, hi + 1):
                rows.append((ch[i], p, st[i], occ[i], None if want is None else want[i]))
        bad = int(np.sum(~np.isfinite(a)))
        if not rows:
            return (np.array([], dtype=object), np.array([], np.int64), np.array([], dtype=object), np.array([]), None, len(data))
        r = list(zip(*rows))
        return (np.array(r[0], dtype=object), np.array(r[1], np.int64), np.array(r[2], dtype=object), np.array(r[3], float),
                None if want is None else np.array(r[4], dtype=object), bad)
    chrom, pos = data[m['chrom']], numeric(data[m['pos']])
    st = strand_of(data[m['strand']]) if m['strand'] else np.array([None] * len(data), dtype=object)
    ok = chrom.notna().to_numpy() & np.isfinite(pos.to_numpy(float)) & (chrom.astype(str).str.strip() != '.').to_numpy()
    return (chrom.astype(str).str.strip().to_numpy(dtype=object)[ok], pos.to_numpy(float)[ok].astype(np.int64), st[ok], occ[ok],
            None if want is None else want[ok], int((~ok).sum()))


# ---------------------------------------------------------------- coordinates
def settle(ref, chrom, pos, strand, mod, base_hint, force_base=None, force_shift=None, sample=20000, want_site=None):
    """Choose 0/1-based and a 1-nt shift by agreement with the modification's base (or with a per-site base
    column, want_site). Returns (pos0, strand, ref base, base, shift, match, table of all match rates)."""
    want = MOD_BASE[mod]
    if want_site is not None:
        want = 'site'
    idx = np.arange(len(pos))
    if len(idx) > sample:
        idx = np.random.default_rng(0).choice(len(pos), sample, replace=False)
    bases = [force_base] if force_base is not None else [0, 1]
    shifts = [force_shift] if force_shift is not None else [0, -1, 1]
    table = {}
    for b in bases:
        for d in shifts:
            st = strand[idx]
            fwd = np.where(st == '-', -d, d)
            p0 = pos[idx] - b + np.where(st == None, 0, fwd)          # noqa: E711
            rb = ref.bases(chrom[idx], p0)
            if want == 'N':
                table[(b, d)] = float(np.mean(rb != 'N'))
                continue
            if want == 'site':
                w = want_site[idx]
                rc_ = np.array([COMP.get(x, 'N') for x in rb], dtype=object)
                table[(b, d)] = float(np.mean(((st == '+') & (rb == w)) | ((st == '-') & (rc_ == w))))
                continue
            plus = (st == '+') & (rb == want)
            minus = (st == '-') & (rb == COMP[want])
            unk = (st == None) & ((rb == want) | (rb == COMP[want]))  # noqa: E711
            table[(b, d)] = float(np.mean(plus | minus | unk))
    best = max(table.values())
    order = sorted(table, key=lambda k: (0 if table[k] >= best - 0.02 else 1, k[1] != 0, k[0] != base_hint, -float(table[k])))
    b, d = order[0]
    st = strand.copy()
    p0 = pos - b + np.where(st == '-', -d, np.where(st == '+', d, 0))
    rb = ref.bases(chrom, p0)
    if want not in ('N', 'site'):
        miss = st == None                                            # noqa: E711
        st[miss & (rb == want)] = '+'
        st[miss & (rb == COMP[want])] = '-'
    return p0, st, rb, b, d, table[(b, d)], table


def convert(path, mod, name, ref, over=None, sheet=None, exclude=None, min_match=0.85, force=False, inspect=False,
            base_col=None, min_sheets=None):
    over = over or {}
    sheets = read_sheets(path)
    tried = []
    for sn, df in sheets.items():
        if df is None or not len(df):
            continue
        data, m, notes = detect(df, over)
        usable = bool((m['chrom'] and m['pos']) or m['siteid'])
        tried.append((sn, len(data), usable, m, notes, data))
    pool = [t for t in tried if t[2] and (not sheet or re.search(sheet, str(t[0])) or len(tried) == 1)]
    if not pool:
        pool = [t for t in tried if t[2]]
    info = [f'  sheet {t[0]!r}: {t[1]:,} rows, ' + ('usable' if t[2] else 'no chrom/position columns') for t in tried]
    if not pool:
        return None, info + ['  -> nothing usable; give --chrom-col / --pos-col (see --inspect)']
    use = [t for t in pool if not sheet or re.search(sheet, str(t[0]))] if min_sheets else []
    if min_sheets and len(use) >= min_sheets:                       # replicate sheets: sites in >= min_sheets of them
        info.append(f'  combining {len(use)} sheets ({", ".join(str(t[0]) for t in use)}): sites in >= {min_sheets}')
    else:
        if min_sheets:
            info.append(f'  fewer than {min_sheets} sheets match {sheet!r}: the largest sheet is used alone')
            min_sheets = None
        use = [max(pool, key=lambda t: t[1])]
    parts = []
    for k, (sn, n, _, m, notes, data) in enumerate(use):
        info.append(f'  using sheet {sn!r}: chrom={m["chrom"] or "-"} pos={m["pos"] or "-"} strand={m["strand"] or "-"} '
                    f'site-id={m["siteid"] or "-"} occupancy={",".join(m["occ"]) or "-"}' + ('; ' + '; '.join(notes) if notes else ''))
        if inspect and k == 0:
            info.append('  columns: ' + ' | '.join(map(str, data.columns)))
            for _, r in data.head(3).iterrows():
                info.append('    ' + ' | '.join(str(x) for x in r.tolist()))
        bc = base_col or over.get('base_col')
        parts.append(sites_from(data, m, bc) + (k, m['base']))
    chrom = np.concatenate([p_[0] for p_ in parts])
    pos = np.concatenate([p_[1] for p_ in parts])
    strand = np.concatenate([p_[2] for p_ in parts])
    occ = np.concatenate([p_[3] for p_ in parts])
    want_site = None if parts[0][4] is None else np.concatenate([p_[4] for p_ in parts])
    bad = sum(p_[5] for p_ in parts)
    sheet_of = np.concatenate([np.full(len(p_[1]), p_[6]) for p_ in parts])
    if want_site is not None:
        info.append(f'  per-site base from column "{base_col or over.get("base_col")}": the convention is checked against it')
    chrom = ref.fix_chrom(chrom)
    p0, st, rb, b, d, match, table = settle(ref, chrom, pos, strand, mod, parts[0][7], over.get('base'), over.get('shift'),
                                            want_site=want_site)
    rates = ', '.join(f'{"0" if k[0] == 0 else "1"}-based shift {k[1]:+d}: {v:.2f}' for k, v in sorted(table.items()))
    info.append(f'  {len(pos):,} sites ({bad:,} rows without a position); reference-base agreement: {rates}')
    checked = mod != 'Nm' or want_site is not None
    info.append(f'  -> {b}-based, shift {d:+d} (' + (f'{100 * match:.1f}% on the expected base' if checked else
                'Nm: any base is possible, so the convention cannot be checked: the file\'s own convention is assumed; '
                'give --base if you know it') + ')')
    if checked and match < min_match and not force:
        info.append(f'  -> NOT WRITTEN: below {min_match:.0%} on the expected base: another build, transcript coordinates, '
                    'or the wrong modification? (lift over, or --force)')
        return None, info
    keep = np.array([s_ in ('+', '-') for s_ in st]) & (rb != 'N')
    rna = np.array([x if s_ == '+' else COMP.get(x, 'N') for x, s_ in zip(rb, st)], dtype=object)
    if want_site is not None:
        keep &= rna == want_site
    elif mod != 'Nm':
        keep &= rna == MOD_BASE[mod]
    occ_ = occ * (100 if np.nanmax(np.append(occ, -1)) <= 1.0 else 1) if np.isfinite(occ).any() else occ
    df = pd.DataFrame({'chrom': chrom[keep], 'start': p0[keep], 'strand': st[keep], 'occ': occ_[keep], 'base': rna[keep],
                       'sheet': sheet_of[keep]})
    df = df.groupby(['chrom', 'start', 'strand', 'base'], as_index=False).agg(occ=('occ', 'mean'), sheets=('sheet', 'nunique'))
    if min_sheets:
        before = len(df)
        df = df[df['sheets'] >= min_sheets]
        info.append(f'  {len(df):,} of {before:,} sites in >= {min_sheets} sheets')
    if exclude:
        ex, xinfo = convert(exclude, mod, 'exclude', ref, over=None, min_match=0.5, force=True)
        if ex is not None:
            k = set(zip(ex['chrom'], ex['start'], ex['strand']))
            before = len(df)
            df = df[[t not in k for t in zip(df['chrom'], df['start'], df['strand'])]]
            info.append(f'  excluded {before - len(df):,} sites also in {os.path.basename(exclude)}')
    info.append(f'  -> {len(df):,} unique sites kept ({int(keep.size - keep.sum()):,} off the expected base or strand dropped)'
                + (f'; occupancy median {np.nanmedian(df["occ"]):.1f}%' if df['occ'].notna().any() else '; no occupancy'))
    df['name'] = name
    return df, info


def write_bed(df, path):
    df = df.sort_values(['chrom', 'start', 'strand'])
    with open(path, 'w') as out:
        for c, s, n, o, st in zip(df['chrom'], df['start'], df['name'], df['occ'], df['strand']):
            out.write(f'{c}\t{s}\t{s + 1}\t{n}\t{"." if o != o else f"{o:.2f}"}\t{st}\n')


def emit(df, mod, name, out):
    """Write one BED (or four for Nm); return [(mod, path, sites)]."""
    os.makedirs(out, exist_ok=True)
    res = []
    if mod == 'Nm':
        for b, sub in df.groupby('base'):
            if b in NM_SPLIT:
                p = os.path.join(out, f'{name}.{NM_SPLIT[b]}.bed')
                write_bed(sub, p)
                res.append((NM_SPLIT[b], p, len(sub)))
    else:
        p = os.path.join(out, f'{name}.{mod}.bed')
        write_bed(df, p)
        res.append((mod, p, len(df)))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('batch')
    b.add_argument('--dir', required=True)
    o = sub.add_parser('one')
    o.add_argument('file')
    o.add_argument('--mod', required=True, choices=sorted(MOD_BASE))
    o.add_argument('--name', required=True)
    o.add_argument('--sheet', help='regex for the sheet to use')
    for k in ('chrom', 'pos', 'end', 'strand'):
        o.add_argument(f'--{k}-col')
    o.add_argument('--occ-col', nargs='+')
    o.add_argument('--base', type=int, choices=[0, 1])
    o.add_argument('--shift', type=int, choices=[-1, 0, 1])
    o.add_argument('--exclude', help='a control set whose sites are removed (same format rules)')
    o.add_argument('--base-col', help='column holding each site\'s modified base (checks the convention, e.g. NJU-seq "Nm")')
    o.add_argument('--min-sheets', type=int, help='combine every sheet matching --sheet; keep sites in at least this many')
    for p in (b, o):
        p.add_argument('--out', required=True)
        p.add_argument('--reference', default=CLUSTER_REFERENCE)
        p.add_argument('--min-match', type=float, default=0.85)
        p.add_argument('--force', action='store_true', help='write a set even below --min-match')
        p.add_argument('--inspect', action='store_true', help='show sheets, columns and first rows; write nothing')
    a = ap.parse_args()
    ref = Ref(a.reference)
    written, report = [], []
    if a.cmd == 'one':
        over = dict(chrom=a.chrom_col, pos=a.pos_col, end=a.end_col, strand=a.strand_col, occ=a.occ_col, base=a.base, shift=a.shift)
        df, info = convert(a.file, a.mod, a.name, ref, over, a.sheet, a.exclude, a.min_match, a.force, a.inspect,
                           a.base_col, a.min_sheets)
        print(f'{os.path.basename(a.file)} -> {a.name} ({a.mod})\n' + '\n'.join(info))
        if df is not None and not a.inspect:
            written += emit(df, a.mod, a.name, a.out)
    else:
        files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(a.dir, '*')) if os.path.isfile(f))
        for f in files:
            preset = next((p for rx, p in PRESETS if re.search(rx, f)), None)
            if preset is None:
                why = next((w for rx, w in SKIP_NOTES if re.search(rx, f)), 'no preset (convert with "one")')
                print(f'{f}: skipped: {why}')
                report.append((f, '.', '.', 'skipped', why))
                continue
            ex = None
            if preset.get('exclude'):
                ex = next((os.path.join(a.dir, g) for g in files if re.search(preset['exclude'], g)), None)
            try:
                df, info = convert(os.path.join(a.dir, f), preset['mod'], preset['name'], ref, {}, preset.get('sheet'), ex,
                                   a.min_match, a.force, a.inspect, preset.get('base_col'), preset.get('min_sheets'))
            except SystemExit as e:
                df, info = None, [f'  -> {e}']
            except Exception as e:                                  # one unreadable file must not stop the batch
                df, info = None, [f'  -> FAILED: {type(e).__name__}: {e}']
            print(f'{f} -> {preset["name"]} ({preset["mod"]})\n' + '\n'.join(info))
            if df is not None and not a.inspect:
                res = emit(df, preset['mod'], preset['name'], a.out)
                written += res
                report += [(f, m_, p_, n_, info[-2].strip()) for m_, p_, n_ in res]
            else:
                report.append((f, preset['mod'], '.', 'not written', info[-1].strip()))
        if not a.inspect:
            os.makedirs(a.out, exist_ok=True)
            pd.DataFrame(report, columns=['file', 'mod', 'bed', 'sites', 'note']).to_csv(os.path.join(a.out, 'report.tsv'), sep='\t', index=False)
    if written:
        print('\n--known-sites ' + ' '.join(f'{m_}={p_}' for m_, p_, n_ in written if n_))


if __name__ == '__main__':
    main()
