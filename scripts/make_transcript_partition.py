#!/usr/bin/env python3
"""make_transcript_partition.py: a strand-aware transcript-region partition of the genome from a GENCODE GTF, one
label per position and strand, as one labelled BED (for cmp_mappability.py --partition).

    python make_transcript_partition.py --gtf gencode.v47.annotation.gtf --fai GRCh38.fa.fai --out regions.bed.gz

Transcripts: those tagged --tag (default 'basic', GENCODE's representative set; most retained-intron and fragment
models are not in it); --tag '' takes every transcript. Per transcript (exons, CDS and stop_codon rows):
  coding span       CDS rows plus stop_codon rows (GENCODE's CDS excludes the stop codon)
  5' UTR            exonic positions upstream of the coding span (in the transcript's direction)
  CDS               exonic positions in the coding span
  stop codon +/-W   exonic positions within W nt of the stop codon along the spliced transcript (first stop
                    nucleotide - W to last + W, clipped to the transcript); coding transcripts without a stop_codon
                    row (CDS end not found) get none
  3' UTR rest       exonic positions downstream of the coding span, outside the stop window
  non-coding exon   exons of transcripts with no CDS row
  intron            within a transcript's first-to-last exon span, in none of its exons
  intergenic        no selected transcript on that strand (every contig of --fai, both strands)
Where transcripts on the same strand disagree, the label of highest priority wins:
  stop codon +/-W > CDS > 5' UTR > 3' UTR rest > non-coding exon > intron > intergenic
Output: BED6 (chrom, start, end, label, 0, strand), sorted by contig (--fai order), start, strand; the whole of every
contig on both strands, no overlaps within a strand. Header lines: '#labels<TAB>...' (display order, 5' to 3', then
non-coding, intron, intergenic) and '#priority<TAB>...'. A summary (bases per label) goes to stderr and <out>.summary.tsv.
"""
import argparse
import gzip
import re
import sys

import numpy as np
import pandas as pd

def labels(window):
    stop = f'stop codon +/-{window} nt'
    display = ["5' UTR", 'CDS', stop, "3' UTR rest", 'non-coding exon', 'intron', 'intergenic']
    priority = [stop, 'CDS', "5' UTR", "3' UTR rest", 'non-coding exon', 'intron', 'intergenic']
    return display, priority


def read_gtf(path, tag):
    cols = ['chrom', 'feature', 'start', 'end', 'strand', 'attr']
    keep = {'transcript', 'exon', 'CDS', 'stop_codon'}
    parts = []
    for ch in pd.read_csv(path, sep='\t', header=None, comment='#', usecols=[0, 2, 3, 4, 6, 8], names=cols,
                          dtype={'chrom': str, 'feature': str, 'start': np.int64, 'end': np.int64, 'strand': str,
                                 'attr': str}, chunksize=1_000_000):
        ch = ch[ch['feature'].isin(keep)].copy()
        ch['tid'] = ch['attr'].str.extract(r'transcript_id "([^"]+)"', expand=False)
        if tag:
            ch['tagged'] = ch['attr'].str.contains(f'tag "{re.escape(tag)}"', regex=True)
        ch = ch.drop(columns='attr')
        parts.append(ch)
    g = pd.concat(parts, ignore_index=True)
    if g['tid'].isna().any():
        sys.exit(f'{path}: {int(g["tid"].isna().sum())} rows without transcript_id')
    tx = g[g['feature'] == 'transcript']
    chosen = set(tx.loc[tx['tagged'], 'tid']) if tag else set(tx['tid'])
    n_all = tx['tid'].nunique()
    g = g[g['tid'].isin(chosen) & (g['feature'] != 'transcript')].drop(columns=[c for c in ['tagged'] if c in g])
    return g, n_all, len(chosen)


def tx_coord(g, ex):
    """transcript coordinate (0-based, 5' to 3') of genomic position g in exons ex (list of (start, end), 1-based
    inclusive, in transcript order) and the strand; None if g is in no exon."""
    o = 0
    for s, e, strand in ex:
        if s <= g <= e:
            return o + (g - s if strand == '+' else e - g)
        o += e - s + 1
    return None


def genomic_of_window(w0, w1, ex):
    """genomic intervals (1-based inclusive) of transcript coordinates w0..w1 over exons ex (transcript order)."""
    out, o = [], 0
    for s, e, strand in ex:
        n = e - s + 1
        a, b = max(w0, o), min(w1, o + n - 1)
        if a <= b:
            out.append((s + (a - o), s + (b - o)) if strand == '+' else (e - (b - o), e - (a - o)))
        o += n
    return out


def build(g, window):
    """intervals per class: DataFrame chrom, strand, start, end (1-based inclusive), cls (priority rank, higher wins)."""
    display, priority = labels(window)
    rank = {lab: len(priority) - 1 - i for i, lab in enumerate(priority)}   # intergenic 0 ... stop window 6
    ex = g[g['feature'] == 'exon']
    cod = g[g['feature'].isin(['CDS', 'stop_codon'])]
    span = ex.groupby('tid').agg(chrom=('chrom', 'first'), strand=('strand', 'first'), start=('start', 'min'),
                                 end=('end', 'max'))
    out = [span.assign(cls=rank['intron'])[['chrom', 'strand', 'start', 'end', 'cls']]]
    cspan = cod.groupby('tid').agg(lo=('start', 'min'), hi=('end', 'max'))
    nc = ex[~ex['tid'].isin(cspan.index)]
    out.append(nc.assign(cls=rank['non-coding exon'])[['chrom', 'strand', 'start', 'end', 'cls']])
    ce = ex[ex['tid'].isin(cspan.index)].join(cspan, on='tid')
    plus = ce['strand'] == '+'
    up = pd.DataFrame({'chrom': ce['chrom'], 'strand': ce['strand'],                       # upstream of the span
                       'start': np.where(plus, ce['start'], np.maximum(ce['start'], ce['hi'] + 1)),
                       'end': np.where(plus, np.minimum(ce['end'], ce['lo'] - 1), ce['end'])})
    down = pd.DataFrame({'chrom': ce['chrom'], 'strand': ce['strand'],
                         'start': np.where(plus, np.maximum(ce['start'], ce['hi'] + 1), ce['start']),
                         'end': np.where(plus, ce['end'], np.minimum(ce['end'], ce['lo'] - 1))})
    mid = pd.DataFrame({'chrom': ce['chrom'], 'strand': ce['strand'], 'start': np.maximum(ce['start'], ce['lo']),
                        'end': np.minimum(ce['end'], ce['hi'])})
    out += [down.assign(cls=rank["3' UTR rest"]), up.assign(cls=rank["5' UTR"]), mid.assign(cls=rank['CDS'])]
    stops = g[g['feature'] == 'stop_codon']
    exs = {t: d for t, d in ex.groupby('tid')}
    win, no_window, off_exon = [], 0, 0
    for t, sc in stops.groupby('tid'):
        d = exs[t]
        strand = d['strand'].iloc[0]
        d = d.sort_values('start', ascending=(strand == '+'))
        el = list(zip(d['start'], d['end'], [strand] * len(d)))
        tc = [tx_coord(p, el) for p in np.concatenate([sc['start'].to_numpy(), sc['end'].to_numpy()])]
        if any(c is None for c in tc):
            off_exon += 1
            continue
        for s, e in genomic_of_window(min(tc) - window, max(tc) + window, el):
            win.append((d['chrom'].iloc[0], strand, s, e))
    no_window = len(cspan) - stops['tid'].nunique()
    out.append(pd.DataFrame(win, columns=['chrom', 'strand', 'start', 'end']).assign(cls=rank[priority[0]]))
    iv = pd.concat(out, ignore_index=True)
    iv = iv[iv['start'] <= iv['end']]
    info = {'transcripts': int(span.shape[0]), 'coding': int(len(cspan)), 'non-coding': int(nc['tid'].nunique()),
            'coding without stop_codon row (no window)': int(no_window),
            'stop codon outside the exons (no window)': int(off_exon)}
    return iv, rank, info


def paint(iv, fai, rank):
    """one label per position and strand: paint classes in increasing priority over intergenic, then run-length."""
    name = {v: k for k, v in rank.items()}
    by = {k: d for k, d in iv.groupby(['chrom', 'strand'])}
    rows = []
    for chrom, length in fai:
        segs = []
        for strand in '+-':
            arr = np.zeros(length, dtype=np.uint8)
            d = by.get((chrom, strand))
            if d is not None:
                for c in sorted(d['cls'].unique()):
                    x = d[d['cls'] == c]
                    s, e = np.clip(x['start'].to_numpy() - 1, 0, length), np.clip(x['end'].to_numpy(), 0, length)
                    for a, b in zip(s, e):
                        arr[a:b] = c
            cut = np.flatnonzero(np.diff(arr)) + 1
            st = np.concatenate([[0], cut]); en = np.concatenate([cut, [length]])
            segs.append(pd.DataFrame({'chrom': chrom, 'start': st, 'end': en, 'label': [name[v] for v in arr[st]],
                                      'score': 0, 'strand': strand}))
        rows.append(pd.concat(segs).sort_values(['start', 'strand'], kind='stable'))
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--fai', required=True, help='contigs and lengths; the partition covers every contig listed')
    ap.add_argument('--tag', default='basic', help="keep transcripts with this GTF tag ('' = all transcripts)")
    ap.add_argument('--window', type=int, default=200, help='nt either side of the stop codon, along the transcript')
    ap.add_argument('--out', required=True, help='.bed or .bed.gz')
    a = ap.parse_args()
    fai = [(l.split('\t')[0], int(l.split('\t')[1])) for l in open(a.fai)]
    g, n_all, n_kept = read_gtf(a.gtf, a.tag)
    missing = sorted(set(g['chrom']) - {c for c, _ in fai})
    if missing:
        sys.exit(f'GTF contigs not in {a.fai}: {", ".join(missing[:10])}')
    iv, rank, info = build(g, a.window)
    P = paint(iv, fai, rank)
    display, priority = labels(a.window)
    op = gzip.open if a.out.endswith('.gz') else open
    with op(a.out, 'wt') as f:
        f.write(f'#transcript-region partition of {a.gtf}; transcripts tagged {a.tag!r} ({n_kept:,} of {n_all:,}); '
                f'window {a.window} nt; one label per position and strand\n')
        f.write('#labels\t' + ','.join(display) + '\n')
        f.write('#priority\t' + ','.join(priority) + '\n')
        P.to_csv(f, sep='\t', header=False, index=False)
    S = P.assign(bases=P['end'] - P['start']).groupby(['label', 'strand'])['bases'].sum().unstack('strand')
    S = S.reindex(display)
    S['fraction'] = S.sum(axis=1) / S.sum(axis=1).sum()
    S.to_csv(a.out.replace('.gz', '').removesuffix('.bed') + '.summary.tsv', sep='\t')
    print(f'transcripts tagged {a.tag!r}: {n_kept:,} of {n_all:,}; ' + '; '.join(f'{k} {v:,}' for k, v in info.items()),
          file=sys.stderr)
    print(f'{len(P):,} intervals on {len(fai)} contigs', file=sys.stderr)
    print(S.to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()
