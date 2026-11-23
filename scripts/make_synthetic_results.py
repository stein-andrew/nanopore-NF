#!/usr/bin/env python3
"""A small synthetic RNA_analysis.nf results tree, in the pipeline's file formats, for smoke-testing
cohort_plots.py and cohort_extract.py end to end (every figure, every optional input).

    make_synthetic_results.py <outdir> [--samples 8] [--seed 7]

Writes <outdir>/results/..., <outdir>/analysis_samplesheet.csv, <outdir>/ref/{genome.fa,.fai,genes.gtf,
givt_9mer.tsv.gz}, <outdir>/known_m6A.bed, <outdir>/rmsk.bed and BAM-like SAM text + IsoQuant read
assignments for cohort_extract.py (<outdir>/extract/). Nothing here is biology: sizes are tiny and the
generative model is only rich enough to exercise every code path (batch, group, composition, isoform
differences, donor SNVs, a duplicated donor, a sex mismatch, clipped and uncorrectable gIVT rows).
"""
import argparse
import gzip
import hashlib
import os

import numpy as np
from scipy.special import expit, logit

MODS = {'a': 'A', '17596': 'A', '69426': 'A', 'm': 'C', '19228': 'C', '17802': 'T', '19227': 'T', '19229': 'G'}
NSITES = {'a': 700, 'm': 250, '17802': 180, '17596': 180, '69426': 80, '19228': 80, '19227': 80, '19229': 80}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
MARKERS = {'neutrophil': ['FCGR3B', 'CSF3R', 'CXCR2', 'CXCR1', 'FPR1', 'AQP9', 'CMTM2', 'PROK2'],
           'monocyte': ['CD14', 'CSF1R', 'FCN1', 'VCAN', 'CD68'],
           'T cell': ['CD3E', 'CD3D', 'CD3G', 'CD2', 'LCK', 'CD28'],
           'B cell': ['MS4A1', 'CD79A', 'CD79B', 'CD19', 'BANK1', 'PAX5'],
           'NK cell': ['GNLY', 'NKG7', 'KLRD1', 'PRF1', 'KLRF1'],
           'platelet': ['PPBP', 'PF4', 'ITGA2B', 'GP9', 'GP1BA', 'TUBB1']}
SPECIAL = ['XIST', 'RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY', 'HBB', 'HBA1', 'HBA2', 'DDX17', 'RPS24']
MT_SITES = [(1076, 'T'), (1486, 'C'), (1488, 'C'), (1592, 'A'), (1593, 'A'), (2617, 'A'), (2815, 'G'), (3039, 'T'),
            (3040, 'G'), (3067, 'T')]
MT_TRUE = {(1488, 'm'): 0.92, (2815, '19229'): 0.95, (3039, '19227'): 0.9, (3040, '19229'): 0.93, (3067, '17802'): 0.94}


def rc(s):
    return ''.join(COMP[b] for b in reversed(s))


def fp_of(kmer, mod):
    h = int(hashlib.md5(f'{kmer}{mod}'.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
    return float(10 ** (-3.6 + 2.2 * h))


def wgz(path, header, rows, sep='\t'):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'wt') as fh:
        if header:
            fh.write(sep.join(header) + '\n')
        for r in rows:
            fh.write(sep.join('' if (isinstance(v, float) and v != v) else (f'{v:.2f}' if isinstance(v, float) else str(v)) for v in r) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('outdir')
    ap.add_argument('--samples', type=int, default=8)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    O = a.outdir
    R = os.path.join(O, 'results')
    ref = 'GRCh38'
    for d in ('ref', 'extract', f'results/cohort/{ref}', f'results/isoforms/{ref}/high_confidence', 'results/read_stats'):
        os.makedirs(os.path.join(O, d), exist_ok=True)
    S = [f'S{i + 1:02d}' for i in range(a.samples)]
    n = len(S)
    batch = {s: ('B1' if i < n // 2 else 'B2') for i, s in enumerate(S)}
    group = {s: ('case' if i % 2 else 'ctrl') for i, s in enumerate(S)}
    donor = {s: i for i, s in enumerate(S)}
    donor[S[-1]] = donor[S[1]]                                       # last sample: same donor as the second
    sex_true = {s: ('F' if donor[s] % 2 == 0 else 'M') for s in S}
    sex_rec = dict(sex_true)
    sex_rec[S[4]] = 'M' if sex_true[S[4]] == 'F' else 'F'            # one recorded wrongly
    depth = {s: rng.uniform(0.6, 1.8) for s in S}
    neut = {s: rng.uniform(0.4, 0.78) for s in S}
    globin = {s: rng.uniform(0.2, 0.7) for s in S}
    thr_a = {s: round(rng.uniform(0.80, 0.88), 2) for s in S}

    # ---------------------------------------------------------------- genome and genes
    L = {'chr1': 320000, 'chr2': 320000}
    seq = {c: list(rng.choice(list('ACGT'), L[c])) for c in L}
    seq['chrM'] = list(rng.choice(list('ACGT'), 16569))
    for p, b in MT_SITES:
        seq['chrM'][p - 1] = b
    seq['chrUn_test'] = list(rng.choice(list('ACGT'), 6000))
    names = [g for v in MARKERS.values() for g in v] + SPECIAL
    names += [f'GENE{i}' for i in range(400 - len(names))]
    genes = []                                                      # (gid, name, chrom, strand, tx list)
    for i, nm in enumerate(names):
        ch = 'chr1' if i < 200 else 'chr2'
        g0 = 2000 + (i % 200) * 1500
        st = '+' if i % 2 == 0 else '-'
        ex = [(g0, g0 + 300), (g0 + 450, g0 + 750), (g0 + 900, g0 + 1200)]
        exb = [(g0, g0 + 300), (g0 + 900, g0 + 1200)]                  # exon-skipping isoform
        genes.append((f'ENSG{i:05d}', nm, ch, st, [(f'ENST{i:05d}a', ex, 'full-splice_match'),
                                                    (f'ENST{i:05d}b', exb, 'novel_in_catalog')]))
    with open(os.path.join(O, 'ref', 'genes.gtf'), 'w') as fh:
        for gid, nm, ch, st, txs in genes:
            fh.write(f'{ch}\tsyn\tgene\t{txs[0][1][0][0] + 1}\t{txs[0][1][-1][1]}\t.\t{st}\t.\tgene_id "{gid}"; gene_type "protein_coding"; gene_name "{nm}";\n')
            for k, (tid, ex, cat) in enumerate(txs):
                tag = ' tag "Ensembl_canonical";' if k == 0 else ''
                attr = f'gene_id "{gid}"; transcript_id "{tid}"; gene_type "protein_coding"; transcript_type "protein_coding"; gene_name "{nm}";{tag}'
                fh.write(f'{ch}\tsyn\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{st}\t.\t{attr}\n')
                for s_, e_ in ex:
                    fh.write(f'{ch}\tsyn\texon\t{s_ + 1}\t{e_}\t.\t{st}\t.\t{attr}\n')
                cds = [(ex[0][0] + 120, ex[0][1]), ex[1], (ex[-1][0], ex[-1][0] + 150)] if k == 0 else [(ex[0][0] + 120, ex[0][1]), (ex[1][0], ex[1][0] + 150)]
                for s_, e_ in cds:
                    fh.write(f'{ch}\tsyn\tCDS\t{s_ + 1}\t{e_}\t.\t{st}\t0\t{attr}\n')

    # ---------------------------------------------------------------- cohort sites (in exons of the canonical isoform)
    slots = []
    for gi, (gid, nm, ch, st, txs) in enumerate(genes):
        for s_, e_ in txs[0][1]:
            for p in range(s_ + 15, e_ - 15, 18):
                slots.append((gi, ch, p, st))
    order = rng.permutation(len(slots))
    sites, used = [], 0
    for code, k in NSITES.items():
        for _ in range(k):
            gi, ch, p, st = slots[order[used]]
            used += 1
            base = MODS[code]
            if st == '+':
                seq[ch][p] = base
            else:
                seq[ch][p] = COMP[base]
            if code == 'a' and rng.random() < 0.75:                    # DRACH for most m6A
                motif = rng.choice(['GGACT', 'AGACA', 'GAACT', 'TGACC'])
                if st == '+':
                    seq[ch][p - 2:p + 3] = list(motif)
                else:
                    seq[ch][p - 2:p + 3] = list(rc(motif))
            pi = rng.beta(1.2, 7) if rng.random() < 0.6 else rng.beta(4, 3)
            sites.append(dict(code=code, gi=gi, chrom=ch, pos=p, strand=st, pi=pi, expr=rng.lognormal(4.2, 0.6),
                              grp=(rng.random() < 0.06), iso=(rng.random() < 0.12)))
    for (p, code), pi in MT_TRUE.items():
        sites.append(dict(code=code, gi=None, chrom='chrM', pos=p - 1, strand='+', pi=pi, expr=4000.0, grp=False, iso=False))
    for k in range(30):                                             # masked contig
        p = 100 + 150 * k
        seq['chrUn_test'][p] = 'A'
        sites.append(dict(code='a', gi=None, chrom='chrUn_test', pos=p, strand='+', pi=0.5, expr=100.0, grp=False, iso=False))
    with open(os.path.join(O, 'ref', 'genome.fa'), 'w') as fh, open(os.path.join(O, 'ref', 'genome.fa.fai'), 'w') as fai:
        off = 0
        for ch in ('chr1', 'chr2', 'chrM', 'chrUn_test'):
            h = f'>{ch}\n'
            fh.write(h)
            off += len(h)
            s_ = ''.join(seq[ch])
            lines = [s_[i:i + 60] for i in range(0, len(s_), 60)]
            fai.write(f'{ch}\t{len(s_)}\t{off}\t60\t61\n')
            for ln in lines:
                fh.write(ln + '\n')
            off += sum(len(x) + 1 for x in lines)
    G = {ch: ''.join(v) for ch, v in seq.items()}

    def kmer(ch, p, st):
        k = G[ch][p - 4:p + 5]
        return rc(k) if st == '-' else k

    # background rows: every position within +-12 nt of a site whose base fits a model on that strand
    rows_bg = set()
    for s in sites:
        for d in range(-12, 13):
            p = s['pos'] + d
            if p < 5 or p > len(G[s['chrom']]) - 6:
                continue
            b = G[s['chrom']][p] if s['strand'] == '+' else COMP[G[s['chrom']][p]]
            for code, mb in MODS.items():
                if mb == b and (d != 0 or code != s['code']):
                    rows_bg.add((s['chrom'], p, s['strand'], code))
    for p, b in MT_SITES:                                           # every model at the mt-rRNA sites
        for code, mb in MODS.items():
            if mb == b and (p - 1, code) not in [(x['pos'], x['code']) for x in sites if x['chrom'] == 'chrM']:
                rows_bg.add(('chrM', p - 1, '+', code))
    site_key = {(s['chrom'], s['pos'], s['strand'], s['code']): s for s in sites}
    anchors = [s for s in sites if s['code'] == 'a' and s['pi'] > 0.5]
    xt_bump = {(s['chrom'], s['pos'] + d * (1 if s['strand'] == '+' else -1), s['strand']): 0.08 for s in anchors for d in (-2, -1, 1, 2)}
    all_rows = sorted(set(site_key) | rows_bg, key=lambda r: (r[0], r[1], r[3], r[2]))
    kmers = {r: kmer(r[0], r[1], r[2]) for r in all_rows}
    with gzip.open(os.path.join(O, 'ref', 'givt_9mer.tsv.gz'), 'wt') as fh:
        cols = ['false_positive_0.50', 'false_positive_0.60', 'false_positive_0.70', 'false_positive_0.80',
                'false_positive_0.90', 'false_positive_0.95']
        fh.write('kmer\tmod\t' + '\t'.join(cols) + '\n')
        seen = set()
        for r in all_rows:
            km, code = kmers[r], r[3]
            if len(km) == 9 and 'N' not in km and (km, code) not in seen and km[4] == MODS[code]:
                seen.add((km, code))
                f_ = fp_of(km, code)
                fh.write(f'{km}\t{code}\t' + '\t'.join(repr(f_) for _ in cols) + '\n')

    # ---------------------------------------------------------------- per-sample cov20_fpcorr.bed
    snv_near = [s for s in sites if s['code'] == 'a'][:40]
    long_vals = {}                                                  # (site key, sample) -> (n, m, pct, corr)
    for j, s in enumerate(S):
        d = os.path.join(R, 'mods', ref, s)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f'{s}.givt_flags.txt'), 'w') as fh:
            fh.write(f'a,{thr_a[s]} 17596,{thr_a[s]} 69426,{thr_a[s]} m,0.85 19228,0.85 19229,0.95 17802,0.9 19227,0.9\n')
        lines = []
        for r in all_rows:
            ch, p, st, code = r
            st_ = site_key.get(r)
            km = kmers[r]
            fp = fp_of(km, code) if len(km) == 9 else 0.01
            if st_ is not None:
                pi = expit(logit(np.clip(st_['pi'], 1e-3, 1 - 1e-3)) + rng.normal(0, 0.35) + (0.9 if st_['grp'] and group[s] == 'case' else 0)
                           + (0.4 * (neut[s] - 0.6) if st_['gi'] is not None and st_['gi'] % 7 == 0 else 0))
                mean = st_['expr'] * depth[s]
            else:
                pi = 0.0
                near = site_key.get((ch, p - (1 if st == '+' else -1), st, 'a'))
                mean = (near or {'expr': 60})['expr'] * depth[s] * rng.uniform(0.7, 1.1) if True else 60
            nv = int(rng.negative_binomial(5, 5 / (5 + mean)))
            if nv < 20:
                continue
            th = 0.95 * pi + fp * (1 - pi) + (xt_bump.get((ch, p, st), 0) if code in ('m', '17802') else 0)
            nm = int(rng.binomial(nv, min(th, 1)))
            pct = round(100 * nm / nv, 2)
            corr = '' if rng.random() < 0.015 else repr(max(0.0, pct - 100 * fp))
            nfail = int(rng.binomial(int(nv / 0.9), 0.1))
            lines.append((ch, p, code, st, nv, pct, nm, nfail, corr))
            if st_ is not None:
                long_vals[(r, s)] = (nv, nm, pct, corr)
        with open(os.path.join(d, f'{s}.cov20_fpcorr.bed'), 'w') as fh:
            for ch, p, code, st, nv, pct, nm, nfail, corr in lines:
                fh.write(f'{ch}\t{p}\t{p + 1}\t{code}\t{nv}\t{st}\t{p}\t{p + 1}\t255,0,0\t{nv}\t{pct:.2f}\t{nm}\t{nv - nm}\t0\t0\t{nfail}\t0\t0\t{corr}\n')
    # ---------------------------------------------------------------- cohort mods matrices (final-stage sites)
    final = sorted({r for (r, s), (nv, nm, pct, corr) in long_vals.items() if corr != '' and float(corr) >= 20},
                   key=lambda r: (r[0], r[1], r[3], r[2]))
    cm = os.path.join(R, 'cohort', ref)
    for metric in ('n_valid', 'n_mod', 'pct_fpcorr'):
        rows = []
        for r in final:
            row = [r[0], r[1], r[1] + 1, r[3], r[2]]
            for s in S:
                v = long_vals.get((r, s))
                if v is None:
                    row.append('')
                else:
                    row.append(v[0] if metric == 'n_valid' else v[1] if metric == 'n_mod' else (f'{float(v[3]):.2f}' if v[3] != '' else ''))
            rows.append(row)
        wgz(os.path.join(cm, f'cohort.{ref}.mods.{metric}.tsv.gz'), ['chrom', 'start', 'end', 'mod', 'strand'] + S, rows)

    # ---------------------------------------------------------------- expression
    gene_counts = np.zeros((len(genes), n))
    base = rng.lognormal(3.5, 1.2, len(genes))
    for gi, (gid, nm, ch, st, txs) in enumerate(genes):
        for j, s in enumerate(S):
            mu = base[gi] * depth[s] * 30
            if nm in MARKERS['neutrophil']:
                mu *= 8 * neut[s]
            elif any(nm in v for k, v in MARKERS.items() if k != 'neutrophil'):
                mu *= 4 * (1 - neut[s])
            if nm == 'XIST':
                mu = (3000 if sex_true[s] == 'F' else 5) * depth[s]
            if nm in ('RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY'):
                mu = (2000 if sex_true[s] == 'M' else 3) * depth[s]
            if nm in ('HBB', 'HBA1', 'HBA2'):
                mu = 2e5 * globin[s] * depth[s]
            if group[s] == 'case' and gi % 50 == 3:
                mu *= 3
            gene_counts[gi, j] = rng.negative_binomial(10, 10 / (10 + mu))
    mt_counts = rng.negative_binomial(10, 10 / (10 + 5e4), size=(5, n))
    allc = np.vstack([gene_counts, mt_counts])
    gnames = [g[1] for g in genes] + [f'MT-CO{i}' for i in range(1, 6)]
    gids = [g[0] for g in genes] + [f'ENSGMT{i}' for i in range(5)]
    cpm = allc / allc.sum(0) * 1e6
    wgz(os.path.join(cm, f'cohort.{ref}.gene_counts.tsv.gz'), ['gene_id', 'gene_name'] + S,
        [[g, nm] + [int(x) for x in allc[i]] for i, (g, nm) in enumerate(zip(gids, gnames))])
    wgz(os.path.join(cm, f'cohort.{ref}.gene_cpm.tsv.gz'), ['gene_id', 'gene_name'] + S,
        [[g, nm] + [f'{x:.4f}' for x in cpm[i]] for i, (g, nm) in enumerate(zip(gids, gnames))])
    mapped = allc.sum(0) * 1.15
    with open(os.path.join(cm, f'cohort.{ref}.read_qc.tsv'), 'w') as fh:
        fh.write('sample\tmapped_reads\tglobin_reads\tmito_reads\tnon_globin_reads\n')
        for j, s in enumerate(S):
            gl = allc[[gnames.index(x) for x in ('HBB', 'HBA1', 'HBA2')], j].sum()
            mt = mt_counts[:, j].sum()
            fh.write(f'{s}\t{int(mapped[j])}\t{int(gl)}\t{int(mt)}\t{int(mapped[j] - gl)}\n')

    # ---------------------------------------------------------------- per-sample gene_mods / feature_mods
    sites_by_key = {(s_['chrom'], s_['pos'], s_['strand'], s_['code']): s_ for s_ in sites}
    for j, s in enumerate(S):
        d = os.path.join(R, 'mods', ref, s)
        gm, fm = [], []
        usage = np.clip(0.7 + 0.25 * (neut[s] - 0.6) * 3, 0.3, 0.95)
        for r in final:
            v = long_vals.get((r, s))
            st_ = sites_by_key[r]
            if v is None or st_['gi'] is None:
                continue
            gid, nm_, ch, strand, txs = genes[st_['gi']]
            nv, nmod, pct, corr = v
            corr = float(corr) if corr != '' else float('nan')
            fp = (pct - corr) / 100 if corr == corr and corr > 0 else 0.002
            if corr == corr and corr >= 20:
                gm.append([r[0], r[1], r[1] + 1, '.', nv, r[2], r[3], gid, gid, nm_, '.', '.', nv, nmod, pct, corr, corr])
            fm.append(['gene', gid, gid, nm_, '.', '.', r[0], r[1], r[2], r[3], nv, nmod, pct, corr, corr, int(corr >= 20)])
            na = int(rng.binomial(nv, usage))
            for tid, cat, nn, extra in ((txs[0][0], txs[0][2], na, 0.0), (txs[1][0], txs[1][2], nv - na, 0.3 if st_['iso'] else 0.0)):
                if nn < 1:
                    continue
                pi_k = np.clip(nmod / nv - fp + extra, 0, 1)
                mk = int(rng.binomial(nn, min(1, 0.95 * pi_k + fp)))
                pk = 100 * mk / nn
                fm.append(['isoform', tid, gid, nm_, tid.replace('ENST', 'REFT'), cat, r[0], r[1], r[2], r[3], nn, mk, pk,
                           max(0.0, pk - 100 * fp), corr, int(nn >= 20 and pk - 100 * fp >= 20)])
        wgz(os.path.join(d, f'{s}.gene_mods.bed'), ['#chrom', 'start', 'end', 'name', 'score', 'strand', 'mod_code', 'feature_id',
                                                    'gene_id', 'ref_gene_name', 'ref_tx_id', 'structural_category', 'n_valid',
                                                    'n_mod', 'percent_modified', 'percent_fpcorr', 'sample_percent_fpcorr'], gm)
        wgz(os.path.join(d, f'{s}.feature_mods.tsv.gz'), ['level', 'feature', 'gene_id', 'ref_gene_name', 'ref_tx_id',
                                                          'structural_category', 'chrom', 'pos', 'strand', 'mod', 'n_valid', 'n_mod',
                                                          'percent', 'percent_fpcorr', 'sample_pct_fpcorr', 'final'], fm)

    # ---------------------------------------------------------------- isoform counts (high confidence)
    hc = os.path.join(R, 'isoforms', ref, 'high_confidence')
    rows, cls = [], []
    for gi, (gid, nm, ch, st, txs) in enumerate(genes):
        tot = gene_counts[gi]
        u = np.array([np.clip(0.75 + (0.6 * (neut[s] - 0.6) if gi % 9 == 0 else 0) + rng.normal(0, 0.02), 0.05, 0.98) for s in S])
        ca = rng.binomial(tot.astype(int), u)
        rows.append([txs[0][0], gid] + list(ca))
        rows.append([txs[1][0], gid] + list(tot.astype(int) - ca))
        cls += [[txs[0][0], txs[0][0].replace('ENST', 'REFT'), txs[0][2]], [txs[1][0], '.', txs[1][2]]]
    wgz(os.path.join(hc, f'cohort.{ref}.counts.tsv'), ['transcript_id', 'gene_id'] + S, rows)
    wgz(os.path.join(hc, f'cohort.{ref}.high_confidence.classification.tsv'), ['isoform_id', 'ref_tx_id', 'structural_category'], cls)

    # ---------------------------------------------------------------- poly(A)
    gene_med, gene_nt = {}, {}
    for j, s in enumerate(S):
        d = os.path.join(R, 'polya', ref, s)
        os.makedirs(d, exist_ok=True)
        reads, per_iso, per_gene = [], {}, {}
        shift = rng.normal(0, 4)
        for k in range(16000):
            if rng.random() < 0.15:
                reads.append([f'{s}r{k}', 'chrM', 'mito', int(rng.gamma(6, 8)), 'used', '.', '.'])
                continue
            gi = int(rng.integers(0, len(genes)))
            gid, nm, ch, st, txs = genes[gi]
            iso = txs[0][0] if rng.random() < 0.7 else txs[1][0]
            m_ = 70 + (gi % 13) * 8 + shift + (40 if (gi % 5 == 0 and iso.endswith('b')) else 0)
            if nm == 'DDX17':
                m_ = 160
            if nm == 'RPS24':
                m_ = 45
            if nm in ('HBB', 'HBA1'):
                m_ = 60
            pt = int(rng.gamma(9, m_ / 9))
            status = 'used' if rng.random() > 0.03 else 'missing'
            gene = gid if rng.random() > 0.05 else '.'
            reads.append([f'{s}r{k}', ch, 'nuclear', pt if status == 'used' else '', status, iso if gene != '.' else '.', gene])
            if status == 'used' and gene != '.':
                per_iso.setdefault((iso, gid, nm, ch, st, txs[0][1][0][0], txs[0][1][-1][1]), []).append(pt)
                per_gene.setdefault((gid, nm, ch, st, txs[0][1][0][0], txs[0][1][-1][1]), []).append(pt)
        wgz(os.path.join(d, f'{s}.polya_reads.tsv.gz'), ['read_id', 'contig', 'compartment', 'pt', 'status', 'isoform', 'gene'], reads)
        wgz(os.path.join(d, f'{s}.polya_features.tsv.gz'), ['level', 'feature_id', 'gene_id', 'ref_gene_name', 'chrom', 'start', 'end',
                                                            'strand', 'n_tails', 'median_pt'],
            [['gene', k[0], k[0], k[1], k[2], k[4], k[5], k[3], len(v), float(np.median(v))] for k, v in per_gene.items()])
        wgz(os.path.join(d, f'{s}.isoform_polya.bed'), ['#chrom', 'start', 'end', 'name', 'score', 'strand', 'feature_id', 'gene_id',
                                                        'ref_gene_name', 'n_tails', 'median_pt', 'ref_tx_id', 'structural_category'],
            [[k[3], k[5], k[6], k[0], len(v), k[4], k[0], k[1], k[2], len(v), float(np.median(v)),
              k[0].replace('ENST', 'REFT') if k[0].endswith('a') else '.', 'full-splice_match' if k[0].endswith('a') else 'novel_in_catalog']
             for k, v in per_iso.items()])
        for k, v in per_gene.items():
            gene_med.setdefault(k, {})[s] = float(np.median(v))
            gene_nt.setdefault(k, {})[s] = len(v)
    keys = sorted(gene_med)
    for metric, src in (('median_pt', gene_med), ('n_tails', gene_nt)):
        wgz(os.path.join(cm, f'cohort.{ref}.polya.gene.{metric}.tsv.gz'),
            ['feature_id', 'gene_id', 'ref_gene_name', 'chrom', 'start', 'end', 'strand'] + S,
            [[k[0], k[0], k[1], k[2], k[4], k[5], k[3]] + [src[k].get(s, 0 if metric == 'n_tails' else float('nan')) for s in S] for k in keys])

    # ---------------------------------------------------------------- read_stats
    for j, s in enumerate(S):
        idn = np.clip(rng.beta(60, 1.2, 3000) - (0.004 if batch[s] == 'B2' else 0), 0.8, 1)
        wgz(os.path.join(R, 'read_stats', f'{s}.read_identity.tsv.gz'), ['identity', 'aligned_len', 'ins_rate', 'del_rate'],
            [[f'{x:.5f}', int(rng.lognormal(6.8, 0.6)), f'{rng.gamma(2, 0.002):.5f}', f'{rng.gamma(2, 0.003):.5f}'] for x in idn])
        sub = []
        for rb in 'ACGT':
            for qb in 'ACGT':
                sub.append([rb, qb, int(1e6 if rb == qb else rng.integers(500, 5000))])
        wgz(os.path.join(R, 'read_stats', f'{s}.substitutions.tsv'), ['ref_base', 'read_base', 'count'], sub)

    # ---------------------------------------------------------------- variants (Clair3-RNA + known_novel)
    common = [('chr1' if k < 300 else 'chr2', int(rng.integers(2000, 300000))) for k in range(600)]
    near = [(st_['chrom'], st_['pos'] + int(rng.choice([-3, -2, 2, 3]))) for st_ in snv_near]
    geno = {dn: {v: rng.choice(['0/1', '1/1', '0/0'], p=[0.4, 0.15, 0.45]) for v in common + near} for dn in set(donor.values())}
    for s in S:
        d = os.path.join(R, 'variants', ref, s)
        os.makedirs(d, exist_ok=True)
        rows = []
        for ch, p in sorted(common + near):
            gt = geno[donor[s]][(ch, p)]
            if gt == '0/0' or rng.random() < 0.2:
                continue
            refb = G[ch][p]
            alt = 'G' if refb != 'G' else 'A'
            rows.append([ch, p + 1, refb, alt, 'PASS', 20, gt, 20, 40, 0.5 if gt == '0/1' else 1.0, 'known'])
        rows.append(['chr1', 5000, 'A', 'G', 'RNAEditing', 15, '0/1', 15, 50, 0.2, 'editing'])
        wgz(os.path.join(d, f'{s}.variants.tsv.gz'), ['chrom', 'pos', 'ref', 'alt', 'filter', 'qual', 'gt', 'gq', 'dp', 'vaf', 'class'], rows)

    # ---------------------------------------------------------------- known sites, repeats, samplesheet
    with open(os.path.join(O, 'known_m6A.bed'), 'w') as fh:
        for st_ in sites:
            if st_['code'] == 'a' and st_['chrom'] != 'chrUn_test' and rng.random() < 0.6:
                fh.write(f"{st_['chrom']}\t{st_['pos']}\t{st_['pos'] + 1}\tglori\t{min(1, max(0, st_['pi'] + rng.normal(0, 0.08))):.3f}\t{st_['strand']}\n")
    with open(os.path.join(O, 'rmsk.bed'), 'w') as fh:
        for st_ in sites:
            if st_['code'] == '17596' and rng.random() < 0.5:
                fh.write(f"{st_['chrom']}\t{st_['pos'] - 100}\t{st_['pos'] + 100}\tAluSx\t0\t+\n")
        fh.write('chr1\t100\t400\tL1PA2\t0\t+\n')
    with open(os.path.join(O, 'analysis_samplesheet.csv'), 'w') as fh:
        fh.write('sample_id,reference,data_type,bam,bai,basecall_model,mod_models,dorado_version,n_runs,batch,group,sex,age,rin\n')
        for i, s in enumerate(S):
            model = 'rna004_130bps_sup@v5.2.0' if i != n - 2 else 'rna004_130bps_sup@v5.1.0'
            fh.write(f'{s},{ref},RNA,bam/{ref}/{s}.{ref}.bam,bam/{ref}/{s}.{ref}.bam.bai,{model},'
                     f'inosine_m6A_2OmeA@v1 pseU_2OmeU@v1 m5C_2OmeC@v1 2OmeG@v1,1.3.0,{1 + i % 2},{batch[s]},{group[s]},'
                     f'{sex_rec[s]},{int(rng.integers(25, 80))},{rng.uniform(6.5, 9.5):.1f}\n')

    # ---------------------------------------------------------------- inputs for cohort_extract.py
    ex = os.path.join(O, 'extract')
    with open(os.path.join(ex, 'catalogue.gtf'), 'w') as fh:
        for gid, nm, ch, st, txs in genes:
            fh.write(f'{ch}\tisomatch\tgene\t{txs[0][1][0][0] + 1}\t{txs[0][1][-1][1]}\t.\t{st}\t.\tgene_id "{gid}";\n')
            for tid, exs, cat in txs:
                fh.write(f'{ch}\tisomatch\ttranscript\t{exs[0][0] + 1}\t{exs[-1][1]}\t.\t{st}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n')
                for s_, e_ in exs:
                    fh.write(f'{ch}\tisomatch\texon\t{s_ + 1}\t{e_}\t.\t{st}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n')
    for j, s in enumerate(S):
        trunc = 0.25 + 0.4 * rng.random()                            # sample-specific 5' truncation
        with gzip.open(os.path.join(ex, f'{s}.read_assignments.tsv.gz'), 'wt') as fa, open(os.path.join(ex, f'{s}.sam'), 'w') as fs:
            fa.write('# IsoQuant read assignments\n#read_id\tchr\tstrand\tisoform_id\tgene_id\tassignment_type\tassignment_events\texons\tadditional_info\n')
            fs.write(f'@HD\tVN:1.6\tSO:coordinate\n@RG\tID:{s}\tSM:{s}\tDS:basecall_model=rna004_130bps_sup@v5.2.0 modbase_models=x\n')
            for k in range(1500):
                gid, nm, ch, st, txs = genes[int(rng.integers(0, len(genes)))]
                tid, exs, cat = txs[0] if rng.random() < 0.7 else txs[1]
                L_ = sum(e - s_ for s_, e in exs)
                cov = L_ if rng.random() > trunc else int(rng.uniform(0.2, 0.9) * L_)
                # genomic blocks covering the 3'-most `cov` transcript bases
                blocks, need = [], cov
                for s_, e_ in (exs[::-1] if st == '+' else exs):
                    if need <= 0:
                        break
                    w = e_ - s_
                    take = min(w, need)
                    blocks.append((e_ - take, e_) if st == '+' else (s_, s_ + take))
                    need -= take
                blocks.sort()
                exs_txt = ','.join(f'{b0 + 1}-{b1}' for b0, b1 in blocks)
                atype = 'unique' if rng.random() > 0.15 else 'ambiguous'
                klass = 'full_splice_match' if cov == L_ else 'incomplete_splice_match'
                fa.write(f'{s}r{k}\t{ch}\t{st}\t{tid}\t{gid}\t{atype}\t.\t{exs_txt}\tgene_assignment={atype}; Classification={klass};\n')
                n_bases = sum(b1 - b0 for b0, b1 in blocks)
                stride, speed = 5, rng.normal(128 if batch[s] == 'B1' else 118, 8)
                nmoves = n_bases
                nsig = int(n_bases / speed * 4000)
                ncol = max(nmoves, nsig // stride)
                mv = np.zeros(ncol, int)
                mv[np.sort(rng.choice(ncol, nmoves, replace=False))] = 1
                mv[0] = 1
                flag = 0 if st == '+' else 16
                cigar = 'N'.join(f'{b1 - b0}M' if i_ == 0 else f'{b0 - blocks[i_ - 1][1]}N{b1 - b0}M' for i_, (b0, b1) in enumerate(blocks)) \
                    if len(blocks) == 1 else ''.join(f'{b1 - b0}M' if i_ == 0 else f'{b0 - blocks[i_ - 1][1]}N{b1 - b0}M' for i_, (b0, b1) in enumerate(blocks))
                seq_ = ''.join(G[ch][b0:b1] for b0, b1 in blocks)
                q = f'{rng.normal(18 if batch[s] == "B1" else 16, 2):.2f}'
                fs.write(f'{s}r{k}\t{flag}\t{ch}\t{blocks[0][0] + 1}\t60\t{cigar}\t*\t0\t0\t{seq_}\t*\tqs:f:{q}\t'
                         f'mv:B:c,{stride},' + ','.join(map(str, mv)) + f'\tpt:i:{int(rng.gamma(9, 10))}\tRG:Z:{s}\n')
    print(f'synthetic results in {R}: {n} samples, {len(final)} final-stage sites, {len(all_rows)} cov20 positions per sample (before depth)')


if __name__ == '__main__':
    main()
