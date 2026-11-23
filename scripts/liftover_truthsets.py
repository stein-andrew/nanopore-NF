#!/usr/bin/env python3
"""liftover_truthsets.py: lift the hg19 truth sets (miCLIP GSE63753, DART-seq GSE125780) to hg38 and verify them.

    python liftover_truthsets.py --chain hg19ToHg38.over.chain.gz --fasta GRCh38.fa --ortho-dir /projects/.../orthogonal

Writes <name>.hg38.bed.gz next to each input (same columns; start/end/strand lifted) and liftover_report.tsv.
Single bases are lifted with the UCSC chain directly: a base maps if it lies in an aligned block; bases covered by more
than one chain (duplicated in hg38) are dropped, as are bases in gaps. Inverted chains flip the strand.
Every site is then checked on hg38, before and after lifting:
  miCLIP   reference base A on the transcript strand (offsets 0 and +-1 reported) and the 5-mer listed in the file
  DART-seq reference base C at the edited site, and A immediately 5' of it on the transcript
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmp_ortho as co  # noqa: E402

COMP = str.maketrans('ACGTN', 'TGCAN')
# file, kind, column holding the listed 5-mer (None: parse the name's suffix; '-' for none)
INPUTS = [('miCLIP_GSE63753_hek293.abcam.CIMS.m6A.9536.bed.txt.gz', 'miclip', 6),
          ('miCLIP_GSE63753_hek293.sysy.CITS.m6A.12051.bed.gz', 'miclip', None),
          ('DARTseq_GSM4024125_HEK293T_APOBEC1YTH_processed.bed.gz', 'dart', '-'),
          ('DARTseq_GSM4024129_HEK293T_APOBEC1YTHmut_processed.bed.gz', 'dart', '-')]


def read_chains(path):
    """[(score, tName, qName, qStrand, qSize, t_starts, q_starts, sizes)] sorted by score, best first."""
    chains, cur, t, q = [], None, 0, 0
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == 'chain':
                if cur:
                    chains.append(cur)
                score, tName, _tSize, _tStrand, tStart, _tEnd, qName, qSize, qStrand, qStart = f[1:11]
                cur = [float(score), tName, qName, qStrand, int(qSize), [], [], []]
                t, q = int(tStart), int(qStart)
            else:
                size = int(f[0])
                cur[5].append(t); cur[6].append(q); cur[7].append(size)
                if len(f) == 3:
                    t += size + int(f[1]); q += size + int(f[2])
    if cur:
        chains.append(cur)
    out = [(s, tn, qn, qs, qz, np.array(ts), np.array(qq), np.array(sz)) for s, tn, qn, qs, qz, ts, qq, sz in chains]
    return sorted(out, key=lambda c: -c[0])


def lift(chains, chrom, start0, strand):
    """Lift 0-based single-base positions. Returns new chrom, start0, strand ('' / -1 where unmapped) and hit counts."""
    n = len(chrom)
    nc, np0, ns = np.full(n, '', dtype=object), np.full(n, -1, dtype=np.int64), np.full(n, '', dtype=object)
    hits = np.zeros(n, dtype=int)
    by_chrom = {c: np.flatnonzero(chrom == c) for c in pd.unique(chrom)}
    flip = {'+': '-', '-': '+', '.': '.'}
    for score, tn, qn, qs, qz, ts, qq, sz in chains:
        idx = by_chrom.get(tn)
        if idx is None or ts[0] > start0[idx].max() or ts[-1] + sz[-1] <= start0[idx].min():
            continue
        p = start0[idx]
        j = np.searchsorted(ts, p, side='right') - 1
        ok = (j >= 0) & (p < ts[np.maximum(j, 0)] + sz[np.maximum(j, 0)])
        if not ok.any():
            continue
        sel, jj = idx[ok], j[ok]
        qpos = qq[jj] + (start0[sel] - ts[jj])
        hits[sel] += 1
        first = sel[hits[sel] == 1]                      # best-scoring chain wins; later hits only mark ambiguity
        qp = qpos[hits[sel] == 1]
        if qs == '-':
            qp = qz - 1 - qp
        nc[first], np0[first] = qn, qp
        ns[first] = [flip[s] if qs == '-' else s for s in strand[first]]
    amb = hits > 1
    nc[amb], np0[amb], ns[amb] = '', -1, ''
    return nc, np0, ns, hits


def five_mer(fasta, chrom, pos1, strand):
    """5-mer centred on each 1-based position, read on the transcript strand."""
    cols = [fasta.bases(chrom, pos1 + d) for d in (-2, -1, 0, 1, 2)]
    s = np.array([''.join(b) for b in zip(*cols)], dtype=object)
    return np.array([x.translate(COMP)[::-1] if st == '-' else x for x, st in zip(s, strand)], dtype=object)


def check(fasta, kind, chrom, pos1, strand, motif):
    """Shares of sites passing the base checks (and the listed motif, for miCLIP)."""
    out = {}
    tx = lambda off: pos1 + np.where(strand == '-', -off, off)
    base = lambda off: fasta.bases(chrom, tx(off))
    plus = strand == '+'
    if kind == 'miclip':
        for off in (0, 1, -1):
            b = base(off)
            out[f'A_at_transcript_{off:+d}' if off else 'A_at_site'] = float(np.mean(np.where(plus, b == 'A', b == 'T')))
        if motif is not None:
            out['listed_5mer_matches'] = float(np.mean(five_mer(fasta, chrom, pos1, strand) == motif))
    else:
        b0, b5 = base(0), base(-1)
        out['C_at_edited_site'] = float(np.mean(np.where(plus, b0 == 'C', b0 == 'G')))
        out['A_5prime_of_C'] = float(np.mean(np.where(plus, b5 == 'A', b5 == 'T')))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--chain', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--ortho-dir', required=True)
    a = ap.parse_args()
    chains = read_chains(a.chain)
    fasta = co.Fasta(a.fasta)
    print(f'{len(chains):,} chains read from {a.chain}', file=sys.stderr)
    rows = []
    for name, kind, mcol in INPUTS:
        path = os.path.join(a.ortho_dir, name)
        if not os.path.exists(path):
            print(f'skip (not found): {path}', file=sys.stderr)
            continue
        d = pd.read_csv(path, sep='\t', header=None, dtype={0: str})
        chrom, start0, strand = d[0].astype(str).to_numpy(), d[1].to_numpy(np.int64), d[5].astype(str).to_numpy()
        motif = None
        if kind == 'miclip':
            motif = (d[mcol].astype(str).str.upper() if mcol is not None else d[3].astype(str).str.rsplit('_', n=1).str[-1].str.upper()).to_numpy()
            motif = np.array([m.replace('U', 'T') for m in motif], dtype=object)
        before = check(fasta, kind, chrom, start0 + 1, strand, motif)
        nc, np0, ns, hits = lift(chains, chrom, start0, strand)
        ok = np0 >= 0
        out = d[ok].copy()
        out[0], out[1], out[2], out[5] = nc[ok], np0[ok], np0[ok] + 1, ns[ok]
        dst = path.replace('.bed.txt.gz', '.hg38.bed.gz').replace('_processed.bed.gz', '_processed.hg38.bed.gz')
        if dst == path:
            dst = path.replace('.bed.gz', '.hg38.bed.gz')
        out.to_csv(dst, sep='\t', header=False, index=False, compression='gzip')
        after = check(fasta, kind, out[0].to_numpy(), out[1].to_numpy(np.int64) + 1, out[5].astype(str).to_numpy(),
                      motif[ok] if motif is not None else None)
        row = {'file': name, 'sites': len(d), 'lifted': int(ok.sum()), 'unmapped': int((hits == 0).sum()),
               'multi_chain_dropped': int((hits > 1).sum()), 'strand_flipped': int((ns[ok] != strand[ok]).sum()),
               **{f'{k} (hg19 coords)': v for k, v in before.items()}, **{f'{k} (lifted)': v for k, v in after.items()},
               'written': os.path.basename(dst)}
        rows.append(row)
    rep = pd.DataFrame(rows)
    rep.to_csv(os.path.join(a.ortho_dir, 'liftover_report.tsv'), sep='\t', index=False, float_format='%.4f')
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(rep.T.to_string(header=False))


if __name__ == '__main__':
    main()
