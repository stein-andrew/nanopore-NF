#!/usr/bin/env python3
"""cmp_ortho.py: score basecalling setups against orthogonal HEK293T modification maps.

    python cmp_ortho.py --sample 08_01_25_RNA004_HEK293 --ortho-dir /projects/Genometechlab/andrew/orthogonal \\
        --fasta GRCh38_no_alt.fa --mask-bed clusters.bed --out ortho/08_01_nofilter \\
        --run A_d212_v6=pileup/A_d212_v6/08_01.bed --run B_d130_v52=... --run C_d212_v52=...

Truth sets (file names as downloaded into --ortho-dir; filters follow Ghohabi Esfahani et al. 2026):
  m6A    GLORI1.0 (Sun 2025, MOESM5, 10 ng sheet), GLORI2.0 (Sun 2025, MOESM3, 10 ng), GLORI1.0_2022 (Liu 2022),
         and the consensus GLORI1.0∩2.0
  Psi    BID-seq (GSE179798), PRAISE (dataset 2, NM_/XM_ only; position ranges kept as candidates), BID∩PRAISE
  m5C    UBS-seq (GSE225614, rRNA/tRNA removed)
  Inosine SLIC-seq (GSE169710 Data S2, HEK293T replicates, UTR/exonic/intergenic): union, and >=2 of 3 replicates
  2'-O-Me NJU-seq (Cell Rep Methods 2024, HEK293T sheet), split by modified base into 2OmeA/C/G/U
  m6A    also m6A-SAC-seq (GSE162356, HEK293, calibrated stoichiometry), miCLIP CIMS and CITS (GSE63753, HEK293),
         DART-seq (GSE125780, HEK293T; A 5' of each edited C, YTHmut background removed), and the consensus
         [miCLIP and DART-seq are hg19 in GEO: run liftover_truthsets.py first; it writes the *.hg38.bed.gz files read here]
         m6A_2+methods: sites found by at least two of GLORI (any version), m6A-SAC-seq, miCLIP, DART-seq
  Psi    also PRAISE_PUS-dependent (dataset 3, sites lost in PUS knockouts)
  m5C    also UBS-seq_NSUN2-dependent: HEK293T sites lost or halved on NSUN2 knockdown in HeLa (cross-cell-line)
Offsets in the reference check are transcript-relative (+1 = one base 3' on either strand).

Before scoring, every set is checked against the reference: at each listed site the reference base must be the
modified base (A for m6A, inosine and Am, C for m5C and Cm, T for Psi and Um, G for Gm), or its complement on the
minus strand. Offsets 0, +1 and -1 are tried; the best offset is applied when at least 90% of sites match, sites
that still mismatch are dropped, and a set below 75% is excluded as being on another build or convention. Sets
without a strand (PRAISE) take the strand implied by the reference base.

Scoring (per modification and truth set, for each run):
  universe     reference positions with >= --min-cov valid reads for that modification in EVERY run
  recall(c)    truth sites in the universe whose occupancy is >= c (for a PRAISE range: the highest candidate)
  precision(c) universe sites at >= c that are truth sites. Truth sets are incomplete, so this is a lower bound
               and is meant for comparing runs, not as an absolute FDR.
  c runs 1..100%; 20% is the lab's operating point. AUPRC is the area under the precision-recall curve.
  Matched recall: the lowest recall any run reaches at 20%; each run's precision at that recall.
  Stoichiometry: occupancy against the set's own level (GLORI, BID-seq, PRAISE, UBS-seq): r, RMSE, mean difference.
  Uncertainty: 95% CIs from a paired block bootstrap (--bootstrap, 100 kb blocks), the same resampled blocks for every run,
  for each run and for each difference between runs (tables/ci.tsv; share_boot_le_0 = bootstrap share of differences <= 0).
  Paper style: Jaccard on (chromosome, position) of each run's valid sites (>= 20 reads, >= 20%) against the full
  set as listed, as in the paper's figure 6, for continuity.
"""
import argparse
import io
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MODS = ['m6A', 'Inosine', '2OmeA', 'm5C', '2OmeC', 'Psi', '2OmeU', '2OmeG']
CODE = {'m6A': 'a', 'Inosine': '17596', '2OmeA': '69426', 'm5C': 'm', '2OmeC': '19228', 'Psi': '17802',
        '2OmeU': '19227', '2OmeG': '19229'}
BASE = {'m6A': 'A', 'Inosine': 'A', '2OmeA': 'A', 'm5C': 'C', '2OmeC': 'C', 'Psi': 'T', '2OmeU': 'T', '2OmeG': 'G'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
SHORT = {'A_d212_v6': 'v6 / 2.1.2', 'B_d130_v52': 'v5.2 / 1.3.0', 'C_d212_v52': 'v5.2 / 2.1.2'}
COLORS = {'A_d212_v6': '#1b9e77', 'B_d130_v52': '#7570b3', 'C_d212_v52': '#d95f02'}
FILES = {'glori25_v2': 'GLORI2025_MOESM3_ESM.xlsb', 'glori25_v1': 'GLORI2025_MOESM5_ESM.xlsb',
         'glori22': 'GLORI2022_41587_2022_1487_MOESM3_ESM.xlsx', 'bid': 'BIDseq_GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx',
         'praise': 'PRAISE_41589_2023_1304_MOESM3_ESM.xlsx', 'ubs': 'UBSseq_GSE225614_HEK293T-WT_sites.tsv.gz',
         'slic': 'SLICseq_GSE169710_Data_S2.xlsx', 'nju': 'NJUseq_mmc3.xlsx',
         'sac': 'm6ASACseq_GSE162356_HEK293_polyA.bed.gz',
         'miclip_cims': 'miCLIP_GSE63753_hek293.abcam.CIMS.m6A.9536.hg38.bed.gz',     # GEO files are hg19:
         'miclip_cits': 'miCLIP_GSE63753_hek293.sysy.CITS.m6A.12051.hg38.bed.gz',     # lifted by liftover_truthsets.py
         'ubs_hela_ctrl': 'UBSseq_GSE225614_HeLa-shControl_sites.tsv.gz',
         'ubs_hela_kd': 'UBSseq_GSE225614_HeLa-shNSUN2_sites.tsv.gz',
         'dart': 'DARTseq_GSM4024125_HEK293T_APOBEC1YTH_processed.hg38.bed.gz',
         'dart_mut': 'DARTseq_GSM4024129_HEK293T_APOBEC1YTHmut_processed.hg38.bed.gz'}
M6A_METHODS = {'GLORI': ['GLORI1.0', 'GLORI2.0', 'GLORI1.0_2022'], 'SAC': ['m6A-SAC-seq'],
               'miCLIP': ['miCLIP-CIMS', 'miCLIP-CITS'], 'DART': ['DART-seq']}                 # independent chemistries for the consensus
CUTOFFS = np.arange(1, 101)
TRAPZ = getattr(np, 'trapezoid', None) or np.trapz            # numpy 2 renamed trapz


# ------------------------------------------------------------------ reading
def read_excel(path, **kw):
    """pyxlsb for .xlsb when installed; otherwise openpyxl from memory (also reads .xlsx content in any file)."""
    if path.endswith('.xlsb'):
        try:
            return pd.read_excel(path, engine='pyxlsb', **kw)
        except ImportError:
            pass
    with open(path, 'rb') as fh:
        return pd.read_excel(io.BytesIO(fh.read()), engine='openpyxl', **kw)


def norm_chrom(c):
    """'1' -> 'chr1', 'MT'/'M' -> 'chrM'; names already starting with chr are left as they are."""
    c = str(c)
    if c in ('MT', 'M', 'chrMT'):
        return 'chrM'
    return 'chr' + c if re.fullmatch(r'[0-9]+|X|Y', c) else c


def frame(set_name, mod, chrom, pos, strand, level, site_id=None):
    n = len(pos)
    col = lambda x: np.full(n, x) if np.ndim(x) == 0 else np.asarray(x)
    d = pd.DataFrame({'set': set_name, 'mod': mod, 'chrom': [norm_chrom(c) for c in col(chrom)],
                      'pos': pd.to_numeric(pd.Series(col(pos)), errors='coerce').to_numpy(),
                      'strand': col(strand).astype(str),
                      'level': pd.to_numeric(pd.Series(col(level)), errors='coerce').to_numpy()})
    d['site_id'] = np.arange(len(d)) if site_id is None else site_id
    d = d.dropna(subset=['pos'])
    d['pos'] = d['pos'].astype(np.int64)
    return d


def praise_frame(d, name, lev):
    """PRAISE chr_site: single positions (strand from the reference base), contiguous U-run ranges (all positions)
    and splice-junction-crossing ranges (the two endpoints only); ranges are written 5'->3' along the transcript."""
    rows, kinds = [], {'single': 0, 'contiguous range': 0, 'junction-crossing range (endpoints only)': 0, 'unparsed': 0}
    for i, (cs, site) in enumerate(zip(d['chr_site'].astype(str), d['site'].astype(str))):
        chrom, _, span = cs.rpartition('_')
        m = re.fullmatch(r'(\d+)(?:-(\d+))?', span)
        if not chrom or not m:
            kinds['unparsed'] += 1
            continue
        a = int(m.group(1))
        if m.group(2) is None:
            rows.append((chrom, a, '.', i, lev.iloc[i])); kinds['single'] += 1
            continue
        b = int(m.group(2)); strand = '+' if b > a else '-'
        t = re.fullmatch(r'(\d+)-(\d+)', site.strip())
        t_len = abs(int(t.group(2)) - int(t.group(1))) + 1 if t else None
        if t_len is not None and abs(b - a) + 1 == t_len:
            rows += [(chrom, q, strand, i, lev.iloc[i]) for q in range(min(a, b), max(a, b) + 1)]
            kinds['contiguous range'] += 1
        else:
            rows += [(chrom, a, strand, i, lev.iloc[i]), (chrom, b, strand, i, lev.iloc[i])]
            kinds['junction-crossing range (endpoints only)'] += 1
    r = pd.DataFrame(rows, columns=['chrom', 'pos', 'strand', 'site_id', 'level'])
    return (frame(name, 'Psi', r['chrom'], r['pos'], r['strand'], r['level'], site_id=r['site_id'].to_numpy()),
            f'{name}: {len(d):,} mRNA sites -> {len(r):,} candidate positions; ' + ', '.join(f'{k} {v:,}' for k, v in kinds.items()))


def load_sets(odir, glori_sheet):
    sets, notes = [], []
    p = lambda k: os.path.join(odir, FILES[k])
    for key, name in [('glori25_v1', 'GLORI1.0'), ('glori25_v2', 'GLORI2.0')]:
        if os.path.exists(p(key)):
            d = read_excel(p(key), sheet_name=glori_sheet)
            lev = d[[c for c in d.columns if str(c).startswith('m6A_level')]].mean(axis=1)
            sets.append(frame(name, 'm6A', d['Chr'], d['Site'], d['Strand'], lev))
    if os.path.exists(p('glori22')):
        d = read_excel(p('glori22'))
        lev = d[[c for c in d.columns if str(c).startswith('m6A_level')]].mean(axis=1) * 100
        sets.append(frame('GLORI1.0_2022', 'm6A', d['Chr'], d['Sites'], d['Strand'], lev))
    if os.path.exists(p('bid')):
        d = read_excel(p('bid'), header=3)
        sets.append(frame('BID-seq', 'Psi', d['chr'], d['pos'], d['strand'], d['Frac_Ave %']))
    if os.path.exists(p('praise')):
        for sheet, name, levcols in [('Supplementary Dataset 2', 'PRAISE', ['rep1_deletion_ratio', 'rep2_deletion_ratio']),
                                     ('Supplementary Dataset 3', 'PRAISE_PUS-dependent', ['wt_treated_deletion_ratio'])]:
            try:
                d = read_excel(p('praise'), sheet_name=sheet, header=2)
            except ValueError:
                continue
            d = d[d['chr_name'].astype(str).str.match(r'^(NM_|XM_)')].reset_index(drop=True)
            fr, note = praise_frame(d, name, d[levcols].mean(axis=1) * 100)
            sets.append(fr); notes.append(note)
    if os.path.exists(p('ubs')):
        d = pd.read_csv(p('ubs'), sep='\t')
        d = d[(d['chromosome'] != '.') & ~d['gene_type'].astype(str).str.contains('rRNA|tRNA', case=False)]
        sets.append(frame('UBS-seq', 'm5C', d['chromosome'], d['position'], d['strand'], d['ratio'] * 100))
    if os.path.exists(p('ubs')) and os.path.exists(p('ubs_hela_ctrl')) and os.path.exists(p('ubs_hela_kd')):
        hek = pd.read_csv(p('ubs'), sep='\t')
        hek = hek[(hek['chromosome'] != '.') & ~hek['gene_type'].astype(str).str.contains('rRNA|tRNA', case=False)]
        ctrl, kd = (pd.read_csv(p(k), sep='\t') for k in ('ubs_hela_ctrl', 'ubs_hela_kd'))
        key = lambda c, q, st: c.astype(str).map(norm_chrom) + ':' + q.astype(str) + ':' + st.astype(str)
        kd_ratio = pd.Series(kd['ratio'].to_numpy(), index=key(kd['ref'], kd['pos'], kd['strand']))
        kd_ratio = kd_ratio[~kd_ratio.index.duplicated()]
        ck = key(ctrl['ref'], ctrl['pos'], ctrl['strand'])
        r_kd = kd_ratio.reindex(ck).to_numpy()
        dep = set(ck[np.isnan(r_kd) | (r_kd <= 0.5 * ctrl['ratio'].to_numpy())])      # lost or halved on NSUN2 knockdown
        hk = key(hek['chromosome'], hek['position'], hek['strand'])
        d = hek[hk.isin(dep).to_numpy()]
        sets.append(frame('UBS-seq_NSUN2-dependent', 'm5C', d['chromosome'], d['position'], d['strand'], d['ratio'] * 100))
        notes.append(f'UBS-seq NSUN2-dependent: {len(dep):,} HeLa sites lost or halved on shNSUN2; {len(d):,} of them in the HEK293T set '
                     '(knockdown data are HeLa: a cross-cell-line subset)')
    if os.path.exists(p('sac')):
        # columns: chrom, start, end, name, mean replicate fraction (%), strand, calibrated stoichiometry (%)
        d = pd.read_csv(p('sac'), sep='\t', header=None, usecols=[0, 2, 5, 6], names=['chrom', 'end', 'strand', 'level'])
        sets.append(frame('m6A-SAC-seq', 'm6A', d['chrom'], d['end'], d['strand'], d['level']))
    if os.path.exists(p('dart')):
        # edited C sites (C-to-U next to m6A); YTHmut edits are binding-independent background and are removed.
        # The m6A is the A immediately 5' of the edited C on the transcript: genomic -1 on '+', +1 on '-'.
        cols = ['chrom', 'start', 'end', 'c4', 'c5', 'strand']
        d = pd.read_csv(p('dart'), sep='\t', header=None, names=cols, usecols=range(6))
        if os.path.exists(p('dart_mut')):
            mut = pd.read_csv(p('dart_mut'), sep='\t', header=None, names=cols, usecols=range(6))
            bg = set(mut['chrom'].astype(str) + ':' + mut['end'].astype(str) + ':' + mut['strand'])
            n0 = len(d)
            d = d[~(d['chrom'].astype(str) + ':' + d['end'].astype(str) + ':' + d['strand']).isin(bg)]
            notes.append(f'DART-seq: {n0:,} APOBEC1-YTH edit sites, {n0 - len(d):,} also edited by APOBEC1-YTHmut removed')
        a_pos = np.where(d['strand'] == '-', d['end'] + 1, d['end'] - 1)
        sets.append(frame('DART-seq', 'm6A', d['chrom'], a_pos, d['strand'], np.nan))
    for key_, name in (('miclip_cims', 'miCLIP-CIMS'), ('miclip_cits', 'miCLIP-CITS')):
        if os.path.exists(p(key_)):
            d = pd.read_csv(p(key_), sep='\t', header=None, usecols=[0, 2, 5], names=['chrom', 'end', 'strand'])
            sets.append(frame(name, 'm6A', d['chrom'], d['end'], d['strand'], np.nan))
    if os.path.exists(p('slic')):
        sheets = read_excel(p('slic'), sheet_name=None)
        reps = []
        for name, d in sheets.items():
            if not name.startswith('HEK293T-rep'):
                continue
            tokens = d['Location'].astype(str).str.split(';')
            keep = tokens.apply(lambda t: any(x in ('UTR3', 'UTR5', 'exonic', 'intergenic') for x in t))
            reps.append(d.loc[keep, ['Chromosome', 'position', 'strand']].drop_duplicates().assign(rep=name))
        if reps:
            u = pd.concat(reps).groupby(['Chromosome', 'position', 'strand'], as_index=False)['rep'].nunique()
            sets.append(frame('SLIC-seq', 'Inosine', u['Chromosome'], u['position'], u['strand'], np.nan))
            v = u[u['rep'] >= 2]
            sets.append(frame('SLIC-seq_2of3', 'Inosine', v['Chromosome'], v['position'], v['strand'], np.nan))
            notes.append(f'SLIC-seq: {len(reps)} HEK293T replicate sheets, {len(u):,} sites in any, {len(v):,} in >=2')
    if os.path.exists(p('nju')):
        d = read_excel(p('nju'), sheet_name='HEK293T', header=1)
        chrom = 'chr' + d['Chr'].astype(str).str.replace(r'^chr', '', regex=True).replace({'MT': 'M'})
        strand = d['Strand'].map({'Forward': '+', 'Reverse': '-'}).fillna('.')
        nm = d['Nm'].astype(str).str.upper().str.replace('U', 'T')
        for b, mod in [('A', '2OmeA'), ('C', '2OmeC'), ('G', '2OmeG'), ('T', '2OmeU')]:
            k = nm == b
            sets.append(frame('NJU-seq', mod, chrom[k], d.loc[k, 'Position'], strand[k], np.nan))
    for s in sets:
        s['site_id'] = s['set'] + ':' + s['mod'] + ':' + s['site_id'].astype(str)
    return pd.concat(sets, ignore_index=True), notes


# ------------------------------------------------------------------ reference check
class Fasta:
    def __init__(self, path):
        self.fh = open(path, 'rb')
        self.idx = {}
        with open(path + '.fai') as fai:
            for line in fai:
                n, length, off, lb, lw = line.split('\t')[:5]
                self.idx[n] = (int(length), int(off), int(lb), int(lw))

    def bases(self, chrom, pos):
        """reference base at 1-based positions; '' where the contig is absent or the position is out of range."""
        out = np.full(len(pos), '', dtype=object)
        for c in pd.unique(chrom):
            if c not in self.idx:
                continue
            length, off, lb, lw = self.idx[c]
            sel = np.flatnonzero(chrom == c)
            for i in sel:
                p0 = pos[i] - 1
                if 0 <= p0 < length:
                    self.fh.seek(off + (p0 // lb) * lw + p0 % lb)
                    out[i] = self.fh.read(1).decode().upper()
        return out


OFFSETS = {'0': (0, False), 'genomic+1': (1, False), 'genomic-1': (-1, False),
           'transcript+1': (1, True), 'transcript-1': (-1, True)}   # transcript: +1 = one base 3' on either strand


def check_against_reference(truth, fasta):
    """Per set, the share of listed positions whose reference base is the modified base (or its complement on '-'),
    as given and under genomic (0/1-based) and transcript-relative (CLIP, deamination) one-base offsets. The best
    offset is applied when >= 90% match (75-90%: warn); sets below 75% are excluded."""
    qc, kept = [], []
    for (set_name, mod), d in truth.groupby(['set', 'mod'], sort=False):
        d = d.copy()
        exp_plus = BASE[mod]
        chrom, pos, st = d['chrom'].to_numpy(), d['pos'].to_numpy(), d['strand'].to_numpy()
        on_ref = np.array([c in fasta.idx for c in chrom])
        res, frac = {}, {}
        for lab, (off, tx) in OFFSETS.items():
            shift = np.where(st == '-', -off, off) if tx else np.full(len(d), off)
            b = fasta.bases(chrom, pos + shift)
            plus, minus = b == exp_plus, b == COMP[exp_plus]
            ok = np.where(st == '+', plus, np.where(st == '-', minus, plus | minus))
            res[lab] = (ok, plus, minus, shift)
            frac[lab] = (ok & on_ref).sum() / max(int(on_ref.sum()), 1)
        best = '0' if frac['0'] >= 0.9 else max(frac, key=frac.get)        # prefer the coordinates as given
        status = 'ok' if frac[best] >= 0.9 else ('warn' if frac[best] >= 0.75 else 'excluded')
        ok, plus, minus, shift = res[best]
        d['pos'] = d['pos'] + shift
        d.loc[d['strand'] == '.', 'strand'] = np.where(plus, '+', '-')[(d['strand'] == '.').to_numpy()]
        keep = ok & on_ref
        row = {'set': set_name, 'mod': mod, 'positions_listed': len(d), 'sites_listed': d['site_id'].nunique(),
               'on_reference_contigs': int(on_ref.sum())}
        row.update({('match_offset_0' if k == '0' else f'match_{k}'): v for k, v in frac.items()})
        row.update({'offset_applied': best, 'status': status,
                    'positions_kept': int(keep.sum()) if status != 'excluded' else 0,
                    'sites_kept': d.loc[keep, 'site_id'].nunique() if status != 'excluded' else 0})
        qc.append(row)
        if status != 'excluded':
            kept.append(d[keep])
    return pd.concat(kept, ignore_index=True) if kept else truth.iloc[0:0], pd.DataFrame(qc)


def add_consensus(truth):
    out = [truth]
    key = lambda d: d['chrom'] + ':' + d['pos'].astype(str) + ':' + d['strand']
    for a, b, name, mod in [('GLORI1.0', 'GLORI2.0', 'GLORI1.0∩2.0', 'm6A'), ('BID-seq', 'PRAISE', 'BID∩PRAISE', 'Psi')]:
        da, db = truth[(truth['set'] == a) & (truth['mod'] == mod)], truth[(truth['set'] == b) & (truth['mod'] == mod)]
        if da.empty or db.empty:
            continue
        hit_b = db.assign(k=key(db))
        d = da.assign(k=key(da)).merge(hit_b[['k', 'level']], on='k', suffixes=('', '_b'))
        d['level'] = d[['level', 'level_b']].mean(axis=1)
        d = d.drop(columns=['k', 'level_b']).drop_duplicates('site_id').assign(set=name)
        d['site_id'] = name + ':' + d['site_id']
        out.append(d)
    m = truth[truth['mod'] == 'm6A'].copy()
    if not m.empty:
        m['k'] = key(m)
        support = {meth: set(m.loc[m['set'].isin(members), 'k']) for meth, members in M6A_METHODS.items()}
        support = {k_: v for k_, v in support.items() if v}
        if len(support) >= 2:
            count = pd.Series(0, index=pd.Index(sorted(set().union(*support.values()))))
            for v in support.values():
                count[list(v)] += 1
            keep = set(count.index[count >= 2])
            d = m[m['k'].isin(keep)].sort_values('level', na_position='last').drop_duplicates('k').copy()
            lev = m[m['k'].isin(keep)].groupby('k')['level'].mean()
            d['level'] = lev.reindex(d['k']).to_numpy()
            d = d.drop(columns='k').assign(set='m6A_2+methods')
            d['site_id'] = 'm6A_2+methods:' + d['site_id']
            out.append(d)
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ nanopore pileups
def read_bed(path, codes, min_cov, contig_re, mask):
    parts = []
    for ch in pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 3, 5, 9, 10], chunksize=4_000_000,
                          names=['chrom', 'start', 'code', 'strand', 'nvalid', 'pct'],
                          dtype={'chrom': str, 'start': np.int64, 'code': str, 'strand': str, 'nvalid': np.int64,
                                 'pct': np.float64}):
        ch = ch[(ch['nvalid'] >= min_cov) & ch['code'].isin(codes)]
        if contig_re:
            ch = ch[~ch['chrom'].str.contains(contig_re, regex=True)]
        parts.append(ch)
    d = pd.concat(parts, ignore_index=True)
    d['pos'] = d['start'] + 1
    if mask is not None and len(d):
        drop = np.zeros(len(d), dtype=bool)
        for c, iv in mask.groupby('chrom'):
            sel = np.flatnonzero(d['chrom'].to_numpy() == c)
            if not len(sel):
                continue
            s, e = iv['start'].to_numpy(), iv['end'].to_numpy()
            o = np.argsort(s); s, e = s[o], np.maximum.accumulate(e[o])
            z = d['start'].to_numpy()[sel]
            j = np.searchsorted(s, z, side='right') - 1
            drop[sel] = (j >= 0) & (z < e[np.clip(j, 0, None)])
        d = d[~drop]
    return d


def keys(chrom, pos, strand, chrom_ids):
    c = pd.Series(chrom).map(chrom_ids).fillna(-1).to_numpy(np.int64)
    return np.where(c < 0, -1, c * (1 << 34) + pos.astype(np.int64) * 2 + (np.asarray(strand) == '-'))


# ------------------------------------------------------------------ scoring
def key_blocks(k, block_size):
    """genomic block of each position key (chromosome id, position // block_size)."""
    return (k >> 34) * (1 << 20) + ((k & ((1 << 34) - 1)) >> 1) // block_size


def ge_counts(block_idx, values, nb):
    """nb x 100 matrix: per block, how many values are >= c for c = 1..100 (values are percentages)."""
    b = np.clip(np.floor(values).astype(np.int64), 0, 100)
    h = np.bincount(block_idx * 101 + b, minlength=nb * 101).reshape(nb, 101)
    return np.cumsum(h[:, ::-1], axis=1)[:, ::-1][:, 1:101].astype(np.float32)


def auprc_rows(rec, prec):
    out = np.full(len(rec), np.nan)
    for i in range(len(rec)):
        ok = np.isfinite(rec[i]) & np.isfinite(prec[i])
        if ok.sum() > 1:
            o = np.argsort(rec[i][ok])
            out[i] = float(TRAPZ(prec[i][ok][o], rec[i][ok][o]))
    return out


def run_pairs(runs):
    pref = [(x, y) for x, y in [('A_d212_v6', 'C_d212_v52'), ('C_d212_v52', 'B_d130_v52'), ('A_d212_v6', 'B_d130_v52')]
            if x in runs and y in runs]
    return pref or [(runs[i], runs[j]) for i in range(len(runs)) for j in range(i + 1, len(runs))]


def score(truth, runs, tables, chrom_ids, min_cov, n_boot=0, block_size=100_000, seed=1):
    curves, summary, stoich, ranked, cis = [], [], [], {}, []
    for (set_name, mod), t in truth.groupby(['set', 'mod'], sort=False):
        code = CODE[mod]
        per_run = {}
        for r in runs:
            d = tables[r][tables[r]['code'] == code]
            k = keys(d['chrom'], d['pos'].to_numpy(), d['strand'], chrom_ids)
            o = np.argsort(k)
            per_run[r] = (k[o], d['pct'].to_numpy()[o])
        universe = per_run[runs[0]][0]
        for r in runs[1:]:
            universe = np.intersect1d(universe, per_run[r][0], assume_unique=True)
        pct = {r: per_run[r][1][np.searchsorted(per_run[r][0], universe)] for r in runs}
        tk = keys(t['chrom'], t['pos'].to_numpy(), t['strand'], chrom_ids)
        in_u = np.isin(tk, universe)
        is_truth = np.isin(universe, tk[in_u])
        site_ids = t['site_id'].to_numpy()
        covered_sites = pd.unique(site_ids[in_u])
        res, site_pct_r = {}, {}
        for r in runs:
            v = pd.Series(pct[r][np.searchsorted(universe, tk[in_u])]).groupby(site_ids[in_u]).max()
            site_pct = v.reindex(covered_sites).to_numpy()
            site_pct_r[r] = site_pct
            rec = np.array([(site_pct >= c).mean() if len(site_pct) else np.nan for c in CUTOFFS])
            called = np.array([(pct[r] >= c).sum() for c in CUTOFFS])
            tp = np.array([((pct[r] >= c) & is_truth).sum() for c in CUTOFFS])
            prec = np.where(called > 0, tp / np.maximum(called, 1), np.nan)
            ok = np.isfinite(rec) & np.isfinite(prec)
            o = np.argsort(rec[ok])
            auprc = float(TRAPZ(prec[ok][o], rec[ok][o])) if ok.sum() > 1 else np.nan
            res[r] = (rec, prec, called, auprc)
            for c, a, b, n in zip(CUTOFFS, rec, prec, called):
                curves.append({'set': set_name, 'mod': mod, 'run': r, 'cutoff_pct': c, 'recall': a, 'precision': b,
                               'n_called': int(n)})
            lev = t.drop_duplicates('site_id').set_index('site_id')['level'].reindex(covered_sites).to_numpy()
            fin = np.isfinite(lev) & np.isfinite(site_pct)
            if fin.sum() >= 10 and np.std(site_pct[fin]) > 0 and np.std(lev[fin]) > 0:
                diff = site_pct[fin] - lev[fin]
                stoich.append({'set': set_name, 'mod': mod, 'run': r, 'n_sites': int(fin.sum()),
                               'pearson_r': np.corrcoef(site_pct[fin], lev[fin])[0, 1],
                               'spearman_rho': pd.Series(site_pct[fin]).corr(pd.Series(lev[fin]), method='spearman'),
                               'rmse': float(np.sqrt(np.mean(diff ** 2))), 'mean_diff_nanopore_minus_set': diff.mean(),
                               'loa_low': diff.mean() - 1.96 * diff.std(ddof=1), 'loa_high': diff.mean() + 1.96 * diff.std(ddof=1)})
        i20 = int(np.flatnonzero(CUTOFFS == 20)[0])
        matched = np.nanmin([res[r][0][i20] for r in runs]) if len(covered_sites) else np.nan
        for r in runs:
            rec, prec, called, auprc = res[r]
            at = np.flatnonzero(rec >= matched) if np.isfinite(matched) else []
            p_m = prec[at[-1]] if len(at) else np.nan
            p20, r20 = prec[i20], rec[i20]
            summary.append({'set': set_name, 'mod': mod, 'run': r, 'sites_after_check': t['site_id'].nunique(),
                            'sites_covered_all_runs': len(covered_sites), 'universe_positions': len(universe),
                            'recall_at_20': r20, 'n_called_at_20': int(called[i20]), 'precision_at_20': p20,
                            'f1_at_20': 2 * p20 * r20 / (p20 + r20) if (p20 + r20) > 0 else np.nan,
                            'auprc': auprc, 'matched_recall': matched, 'precision_at_matched_recall': p_m})
        ranked[(set_name, mod)] = (res, len(covered_sites))
        if n_boot > 0 and len(covered_sites) >= 10:
            cis += bootstrap(set_name, mod, runs, universe, pct, is_truth, tk[in_u], site_ids[in_u], covered_sites,
                             site_pct_r, res, i20, n_boot, block_size, seed)
    return pd.DataFrame(curves), pd.DataFrame(summary), pd.DataFrame(stoich), ranked, pd.DataFrame(cis)


def bootstrap(set_name, mod, runs, universe, pct, is_truth, tk_u, sid_u, covered_sites, site_pct_r, res, i20,
              n_boot, block_size, seed):
    """Paired block bootstrap: the same resampled genomic blocks for every run, so differences get paired CIs."""
    bu_raw = key_blocks(universe, block_size)
    site_block = pd.Series(key_blocks(tk_u, block_size)).groupby(sid_u).min().reindex(covered_sites).to_numpy()
    ub, inv = np.unique(np.concatenate([bu_raw, site_block]), return_inverse=True)
    nb = len(ub); bu, bs = inv[:len(universe)], inv[len(universe):]
    W = np.random.default_rng(seed).multinomial(nb, np.full(nb, 1 / nb), size=n_boot).astype(np.float32)
    nsite = W @ np.bincount(bs, minlength=nb).astype(np.float32)
    stats = {}
    for r in runs:
        called = W @ ge_counts(bu, pct[r], nb)
        tp = W @ ge_counts(bu[is_truth], pct[r][is_truth], nb)
        hits = W @ ge_counts(bs, site_pct_r[r], nb)
        with np.errstate(divide='ignore', invalid='ignore'):
            rec = hits / nsite[:, None]
            prec = np.where(called > 0, tp / np.maximum(called, 1e-9), np.nan)
        stats[r] = {'recall_at_20': rec[:, i20], 'precision_at_20': prec[:, i20], 'auprc': auprc_rows(rec, prec)}
    point = {r: {'recall_at_20': res[r][0][i20], 'precision_at_20': res[r][1][i20], 'auprc': res[r][3]} for r in runs}
    rows = []
    for m in ['recall_at_20', 'precision_at_20', 'auprc']:
        for r in runs:
            lo, hi = np.nanpercentile(stats[r][m], [2.5, 97.5])
            rows.append({'set': set_name, 'mod': mod, 'metric': m, 'comparison': r, 'estimate': point[r][m],
                         'ci_low': lo, 'ci_high': hi, 'share_boot_le_0': np.nan, 'blocks': nb, 'n_boot': n_boot})
        for x, y in run_pairs(runs):
            d = stats[x][m] - stats[y][m]
            lo, hi = np.nanpercentile(d, [2.5, 97.5])
            rows.append({'set': set_name, 'mod': mod, 'metric': m, 'comparison': f'{x} - {y}',
                         'estimate': point[x][m] - point[y][m], 'ci_low': lo, 'ci_high': hi,
                         'share_boot_le_0': float(np.mean(d[np.isfinite(d)] <= 0)) if np.isfinite(d).any() else np.nan,
                         'blocks': nb, 'n_boot': n_boot})
    return rows


def paper_style(truth_listed, runs, tables, min_pct):
    rows = []
    for (set_name, mod), t in truth_listed.groupby(['set', 'mod'], sort=False):
        ts = set(zip(t['chrom'], t['pos']))
        for r in runs:
            d = tables[r][(tables[r]['code'] == CODE[mod]) & (tables[r]['pct'] >= min_pct)]
            ds = set(zip(d['chrom'], d['pos']))
            both = len(ts & ds)
            rows.append({'set': set_name, 'mod': mod, 'run': r, 'nanopore_sites': len(ds), 'set_sites': len(ts),
                         'shared': both, 'jaccard': both / max(len(ts | ds), 1)})
    return pd.DataFrame(rows)


def plot_curves(ranked, runs, sample, path):
    keys_ = [k for k in ranked if ranked[k][1] > 0]
    if not keys_:
        return
    n = len(keys_); ncol = min(4, n); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow), squeeze=False, layout='constrained')
    i20 = int(np.flatnonzero(CUTOFFS == 20)[0])
    for ax, k in zip(axes.flat, keys_):
        res, ncov = ranked[k]
        for r in runs:
            rec, prec, _, auprc = res[r]
            ax.plot(rec, prec, color=COLORS.get(r), lw=1.6, label=f'{SHORT.get(r, r)} (AUPRC {auprc:.2f})')
            ax.plot(rec[i20], prec[i20], 'o', color=COLORS.get(r), ms=6, mec='black', mew=0.6)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('recall'); ax.set_ylabel('precision (lower bound)')
        ax.set_title(f'{k[1]} vs {k[0]} (n={ncov:,})', fontsize=10)
        ax.legend(fontsize=7, frameon=False, loc='upper right')
    for ax in list(axes.flat)[len(keys_):]:
        ax.axis('off')
    fig.suptitle(f'{sample}: precision-recall against orthogonal sets (circles: 20% occupancy)', fontsize=12)
    fig.savefig(path); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--run', action='append', required=True, metavar='NAME=BED')
    ap.add_argument('--ortho-dir', required=True)
    ap.add_argument('--fasta', required=True, help='reference FASTA with .fai (the one the pileups used)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-cov', type=int, default=20)
    ap.add_argument('--min-pct', type=float, default=20)
    ap.add_argument('--glori-sheet', default='10ng_mRNA_input')
    ap.add_argument('--bootstrap', type=int, default=500, help='paired block-bootstrap replicates for CIs (0: none)')
    ap.add_argument('--block-size', type=int, default=100_000, help='genomic block size for the bootstrap (bp)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--style')
    a = ap.parse_args()
    if a.style:
        plt.style.use(a.style)
    runs, beds = [], {}
    for item in a.run:
        name, _, path = item.partition('=')
        if not os.path.exists(path):
            sys.exit(f'cmp_ortho.py: --run {item}: no such file')
        runs.append(name); beds[name] = path
    for f in (a.fasta, a.fasta + '.fai'):
        if not os.path.exists(f):
            sys.exit(f'cmp_ortho.py: {f} not found')
    os.makedirs(f'{a.out}/tables', exist_ok=True); os.makedirs(f'{a.out}/figures', exist_ok=True)

    listed, notes = load_sets(a.ortho_dir, a.glori_sheet)
    for n in notes:
        print(n, file=sys.stderr)
    truth, qc = check_against_reference(listed, Fasta(a.fasta))
    truth = add_consensus(truth)
    qc.to_csv(f'{a.out}/tables/orthogonal_qc.tsv', sep='\t', index=False, float_format='%.4f')
    print(qc[['set', 'mod', 'sites_listed', 'match_offset_0', 'match_genomic+1', 'match_transcript+1', 'match_transcript-1',
              'offset_applied', 'status', 'sites_kept']].to_string(index=False), file=sys.stderr)

    mask = pd.read_csv(a.mask_bed, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                       comment='#') if a.mask_bed else None
    codes = sorted({CODE[m] for m in truth['mod'].unique()})
    tables = {}
    for r in runs:
        tables[r] = read_bed(beds[r], codes, a.min_cov, a.contig_regex, mask)
        print(f'{r}: {len(tables[r]):,} rows with >= {a.min_cov} valid reads', file=sys.stderr)
    chroms = sorted(set().union(*[set(tables[r]['chrom'].unique()) for r in runs]) | set(truth['chrom'].unique()))
    chrom_ids = {c: i for i, c in enumerate(chroms)}

    curves, summary, stoich, ranked, cis = score(truth, runs, tables, chrom_ids, a.min_cov, a.bootstrap, a.block_size, a.seed)
    cis.to_csv(f'{a.out}/tables/ci.tsv', sep='\t', index=False, float_format='%.4g')
    curves.to_csv(f'{a.out}/tables/pr_curves.tsv', sep='\t', index=False, float_format='%.5g')
    summary.to_csv(f'{a.out}/tables/summary.tsv', sep='\t', index=False, float_format='%.4g')
    stoich.to_csv(f'{a.out}/tables/stoichiometry.tsv', sep='\t', index=False, float_format='%.4g')
    paper_style(listed[listed['set'].isin(qc.loc[qc['status'] != 'excluded', 'set'])], runs, tables,
                a.min_pct).to_csv(f'{a.out}/tables/paper_style_jaccard.tsv', sep='\t', index=False, float_format='%.4g')
    plot_curves(ranked, runs, a.sample, f'{a.out}/figures/pr_curves.pdf')
    with open(f'{a.out}/legends.md', 'w') as fh:
        fh.write(f'# {a.sample}: orthogonal comparison\n\nRuns: ' + ', '.join(f'{r} = {beds[r]}' for r in runs) +
                 f'\n\nMinimum valid reads {a.min_cov} in every run; operating point {a.min_pct:g}%; mask '
                 f'{a.mask_bed or "none"}; contigs dropped matching `{a.contig_regex}`; GLORI sheet {a.glori_sheet}.\n\n'
                 'Definitions and filters: see the header of cmp_ortho.py. Precision is a lower bound because every '
                 'truth set is incomplete; compare runs with it rather than reading it as an FDR.\n\n' +
                 '\n'.join(f'- {n}' for n in notes) + '\n')
    print(f'wrote {a.out}/tables, figures/pr_curves.pdf and legends.md', file=sys.stderr)


if __name__ == '__main__':
    main()
