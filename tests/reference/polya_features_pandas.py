#!/usr/bin/env python3
"""Poly(A) tail lengths per read, per isoform and per gene, by SALRR's rules.

    polya_features.py --reads S.reads_pt.tsv --features S.read_features.tsv \\
        --catalogue cohort.GRCh38.catalogue.gtf --classification cohort.GRCh38.classify.classification.txt.gz \\
        --mito-contig chrM --cap-nuclear 0 --cap-mito 0 --min-reads 10 --prefix S

--reads     read_id, contig, pt for every primary alignment (Dorado's pt:i; 'NA' if absent)
--features  read_id, isoform, gene from bin/read_features.py (isoform = unique assignment,
            gene = a single gene), the same reads the isoform/gene modification levels use

A read's tail is used unless, as in SALRR, it is missing (no pt), not found (pt <= 0), or
censored (pt >= its compartment's cap, the estimator's maximum). Reads on --mito-contig are
'mito', all others 'nuclear', with separate caps; 0 means no cap. Dorado's direct-RNA
estimates have none (the highest values are single reads, no pile at a maximum); SALRR's
cDNA estimates were capped at 400 nuclear / 250 mito, so set those for cDNA data.

Writes:
  <prefix>.polya_summary.tsv      per compartment: reads, excluded per reason, used, median, mean, IQR
  <prefix>.polya_hist.tsv         compartment, pt, reads: the per-sample histogram of used tails
  <prefix>.polya_reads.tsv.gz     read_id, contig, compartment, pt, status, isoform, gene
  <prefix>.polya_features.tsv.gz  every isoform and gene with >= 1 used tail: n, median, mean,
                                  q25, q75, min, max (quantiles linearly interpolated)
  <prefix>.isoform_polya.bed      isoforms with >= --min-reads tails, annotated ('#' header)
  <prefix>.gene_polya.bed         the same per gene
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

from feature_mods_pandas import annotations

TX_ID = re.compile(r'transcript_id "([^"]+)"')
GENE_ID = re.compile(r'gene_id "([^"]+)"')
STATS = ['n_tails', 'median_pt', 'mean_pt', 'q25_pt', 'q75_pt', 'min_pt', 'max_pt']
BED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature_id', 'gene_id', 'ref_gene_name',
            'ref_tx_id', 'structural_category', 'compartment'] + STATS


def spans(catalogue):
    """feature id -> (chrom, 0-based start, end, strand), for transcripts and genes."""
    out = {}
    with open(catalogue) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9 or c[2] not in ('transcript', 'gene'):
                continue
            m = (TX_ID if c[2] == 'transcript' else GENE_ID).search(c[8])
            if m:
                out[m.group(1)] = (c[0], int(c[3]) - 1, int(c[4]), c[6])
    return out


def describe(df, by):
    g = df.groupby(by)['pt']
    return pd.DataFrame({'n_tails': g.size(), 'median_pt': g.median(), 'mean_pt': g.mean(),
                         'q25_pt': g.quantile(0.25), 'q75_pt': g.quantile(0.75),
                         'min_pt': g.min(), 'max_pt': g.max()})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reads', required=True)
    ap.add_argument('--features', required=True)
    ap.add_argument('--catalogue', required=True)
    ap.add_argument('--classification', required=True)
    ap.add_argument('--mito-contig', default='chrM')
    ap.add_argument('--cap-nuclear', type=float, default=0, help='nuclear tails >= this are censored (0 = no cap)')
    ap.add_argument('--cap-mito', type=float, default=0, help='mito tails >= this are censored (0 = no cap)')
    ap.add_argument('--min-reads', type=int, default=10, help='tails a feature needs to be in the BEDs')
    ap.add_argument('--prefix', required=True)
    args = ap.parse_args()

    reads = pd.read_csv(args.reads, sep='\t', header=None, names=['read_id', 'contig', 'pt'],
                        dtype={'read_id': str, 'contig': str, 'pt': str}, keep_default_na=False, na_filter=False)
    if reads.empty:
        sys.exit(f"polya_features.py: {args.reads} has no reads")
    pt = pd.to_numeric(reads['pt'].where(reads['pt'] != 'NA'), errors='coerce')
    reads['pt'] = pt
    reads['compartment'] = np.where(reads['contig'] == args.mito_contig, 'mito', 'nuclear')
    cap = np.where(reads['compartment'] == 'mito', args.cap_mito, args.cap_nuclear)
    reads['status'] = np.select([pt.isna(), pt <= 0, (cap > 0) & (pt >= cap)],
                                ['missing', 'not_found', 'censored'], default='used')

    if os.path.getsize(args.features) > 0:
        feats = pd.read_csv(args.features, sep='\t', header=None, names=['read_id', 'isoform', 'gene'],
                            dtype=str, keep_default_na=False, na_filter=False)
    else:
        feats = pd.DataFrame(columns=['read_id', 'isoform', 'gene'])
    if feats['read_id'].duplicated().any():
        sys.exit(f"polya_features.py: {args.features} lists a read more than once")
    reads = reads.merge(feats, on='read_id', how='left').fillna({'isoform': '.', 'gene': '.'})
    if reads['read_id'].duplicated().any():
        sys.exit(f"polya_features.py: {args.reads} lists a read more than once (only primary alignments expected)")
    used = reads[reads['status'] == 'used']
    p = args.prefix

    with open(f"{p}.polya_summary.tsv", 'w') as out:
        out.write("compartment\treads\tmissing\tnot_found\tcensored\tused\tmedian_pt\tmean_pt\tq25_pt\tq75_pt\n")
        for comp in ('nuclear', 'mito'):
            r = reads[reads['compartment'] == comp]
            u = r[r['status'] == 'used']['pt']
            n = r['status'].value_counts()
            stats = '\t'.join(f"{v:.2f}" if len(u) else 'NA'
                              for v in ([u.median(), u.mean(), u.quantile(0.25), u.quantile(0.75)] if len(u) else [0] * 4))
            out.write(f"{comp}\t{len(r)}\t{n.get('missing', 0)}\t{n.get('not_found', 0)}\t{n.get('censored', 0)}\t"
                      f"{len(u)}\t{stats}\n")

    hist = used.groupby(['compartment', 'pt']).size().rename('reads').reset_index()
    hist['pt'] = hist['pt'].astype(np.int64) if (hist['pt'] % 1 == 0).all() else hist['pt']
    hist.to_csv(f"{p}.polya_hist.tsv", sep='\t', index=False)

    out = reads[['read_id', 'contig', 'compartment', 'pt', 'status', 'isoform', 'gene']].copy()
    out['pt'] = out['pt'].map(lambda v: '' if v != v else (str(int(v)) if v == int(v) else repr(v)))
    out.to_csv(f"{p}.polya_reads.tsv.gz", sep='\t', index=False, compression='gzip')

    tx_gene, iso_info, gene_names = annotations(args.catalogue, args.classification)
    where = spans(args.catalogue)
    tables = []
    for level, col in (('isoform', 'isoform'), ('gene', 'gene')):
        d = describe(used[used[col] != '.'], col).reset_index().rename(columns={col: 'feature_id'})
        d.insert(0, 'level', level)
        if level == 'isoform':
            info = d['feature_id'].map(lambda f: iso_info.get(f, ('.', '.', '.')))
            d['gene_id'] = d['feature_id'].map(lambda f: tx_gene.get(f, '.'))
            d['ref_gene_name'], d['ref_tx_id'], d['structural_category'] = info.str[0], info.str[1], info.str[2]
        else:
            d['gene_id'] = d['feature_id']
            d['ref_gene_name'] = d['feature_id'].map(lambda g: gene_names.get(g, '.') or '.')
            d['ref_tx_id'] = d['structural_category'] = '.'
        loc = d['feature_id'].map(lambda f: where.get(f, ('.', -1, -1, '.')))
        d['chrom'], d['start'], d['end'], d['strand'] = loc.str[0], loc.str[1], loc.str[2], loc.str[3]
        d['compartment'] = np.where(d['chrom'] == args.mito_contig, 'mito', 'nuclear')
        tables.append(d)
    features = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    unplaced = features[features['start'] < 0] if len(features) else features
    if len(unplaced):
        sys.exit(f"polya_features.py: {len(unplaced)} feature(s) not in {args.catalogue}, e.g. {list(unplaced['feature_id'][:3])}")
    features = features.sort_values(['chrom', 'start', 'level', 'feature_id']).reset_index(drop=True)
    features.to_csv(f"{p}.polya_features.tsv.gz", sep='\t', index=False, float_format='%.2f', compression='gzip')

    for level in ('isoform', 'gene'):
        sub = features[(features['level'] == level) & (features['n_tails'] >= args.min_reads)]
        with open(f"{p}.{level}_polya.bed", 'w') as fh:
            fh.write('#' + '\t'.join(BED_COLS) + '\n')
            for r in sub.itertuples(index=False):
                vals = [r.chrom, r.start, r.end, r.feature_id, min(int(r.n_tails), 1000), r.strand, r.feature_id,
                        r.gene_id, r.ref_gene_name, r.ref_tx_id, r.structural_category, r.compartment, int(r.n_tails)]
                vals += [f"{v:.2f}" for v in (r.median_pt, r.mean_pt, r.q25_pt, r.q75_pt, r.min_pt, r.max_pt)]
                fh.write('\t'.join(map(str, vals)) + '\n')

    counts = reads['status'].value_counts()
    print(f"polya_features.py: {len(reads)} reads; used {counts.get('used', 0)}, not found {counts.get('not_found', 0)}, "
          f"missing {counts.get('missing', 0)}, censored {counts.get('censored', 0)}; "
          f"isoforms in BED {int(((features['level'] == 'isoform') & (features['n_tails'] >= args.min_reads)).sum()) if len(features) else 0}, "
          f"genes {int(((features['level'] == 'gene') & (features['n_tails'] >= args.min_reads)).sum()) if len(features) else 0}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
