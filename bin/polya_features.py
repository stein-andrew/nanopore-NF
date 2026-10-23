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
import re
import sys

import polars as pl

from feature_mods import annotations
from tabular import read_tsv, write_tsv

TX_ID = re.compile(r'transcript_id "([^"]+)"')
GENE_ID = re.compile(r'gene_id "([^"]+)"')
STATS = ['n_tails', 'median_pt', 'mean_pt', 'q25_pt', 'q75_pt', 'min_pt', 'max_pt']
BED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'feature_id', 'gene_id', 'ref_gene_name',
            'ref_tx_id', 'structural_category', 'compartment'] + STATS
FEATURE_COLS = ['level', 'feature_id'] + STATS + ['gene_id', 'ref_gene_name', 'ref_tx_id', 'structural_category',
                                                   'chrom', 'start', 'end', 'strand', 'compartment']


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
    """Per `by`: tails used, median, mean, linearly interpolated quartiles (as pandas), min, max."""
    pt = pl.col('pt')
    return df.group_by(by).agg(pl.len().alias('n_tails'), pt.median().alias('median_pt'), pt.mean().alias('mean_pt'),
                               pt.quantile(0.25, interpolation='linear').alias('q25_pt'),
                               pt.quantile(0.75, interpolation='linear').alias('q75_pt'),
                               pt.min().alias('min_pt'), pt.max().alias('max_pt'))


def frame(mapping, key, cols):
    """dict key -> tuple of values as a DataFrame, for joining."""
    data = {key: list(mapping)}
    for i, c in enumerate(cols):
        data[c] = [v[i] for v in mapping.values()]
    return pl.DataFrame(data, schema={key: pl.Utf8, **{c: (pl.Int64 if c in ('start', 'end') else pl.Utf8) for c in cols}})


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

    reads = read_tsv(args.reads, ['read_id', 'contig', 'pt'])
    if reads.height == 0:
        sys.exit(f"polya_features.py: {args.reads} has no reads")
    pt = pl.col('pt').cast(pl.Float64, strict=False)          # 'NA' (no pt tag) -> null
    comp = pl.when(pl.col('contig') == args.mito_contig).then(pl.lit('mito')).otherwise(pl.lit('nuclear'))
    reads = reads.with_row_index('_i').with_columns(pt.alias('pt'), comp.alias('compartment'))
    cap = pl.when(pl.col('compartment') == 'mito').then(pl.lit(args.cap_mito)).otherwise(pl.lit(args.cap_nuclear))
    reads = reads.with_columns(
        pl.when(pl.col('pt').is_null()).then(pl.lit('missing'))
          .when(pl.col('pt') <= 0).then(pl.lit('not_found'))
          .when((cap > 0) & (pl.col('pt') >= cap)).then(pl.lit('censored'))
          .otherwise(pl.lit('used')).alias('status'))

    feats = read_tsv(args.features, ['read_id', 'isoform', 'gene'])
    if feats['read_id'].is_duplicated().any():
        sys.exit(f"polya_features.py: {args.features} lists a read more than once")
    if reads['read_id'].is_duplicated().any():
        sys.exit(f"polya_features.py: {args.reads} lists a read more than once (only primary alignments expected)")
    reads = reads.join(feats, on='read_id', how='left').sort('_i') \
                 .with_columns(pl.col('isoform').fill_null('.'), pl.col('gene').fill_null('.'))
    used = reads.filter(pl.col('status') == 'used')
    p = args.prefix

    with open(f"{p}.polya_summary.tsv", 'w') as out:
        out.write("compartment\treads\tmissing\tnot_found\tcensored\tused\tmedian_pt\tmean_pt\tq25_pt\tq75_pt\n")
        for c in ('nuclear', 'mito'):
            r = reads.filter(pl.col('compartment') == c)
            u = r.filter(pl.col('status') == 'used')['pt']
            n = dict(r.group_by('status').agg(pl.len()).iter_rows())
            stats = '\t'.join(f"{v:.2f}" for v in (u.median(), u.mean(), u.quantile(0.25, interpolation='linear'),
                                                    u.quantile(0.75, interpolation='linear'))) if u.len() else '\t'.join(['NA'] * 4)
            out.write(f"{c}\t{r.height}\t{n.get('missing', 0)}\t{n.get('not_found', 0)}\t{n.get('censored', 0)}\t"
                      f"{u.len()}\t{stats}\n")

    hist = used.group_by(['compartment', 'pt']).agg(pl.len().alias('reads')).sort(['compartment', 'pt'])
    if hist.height and (hist['pt'] == hist['pt'].floor()).all():
        hist = hist.with_columns(pl.col('pt').cast(pl.Int64))
    write_tsv(hist, f"{p}.polya_hist.tsv")

    pt_text = pl.when(pl.col('pt').is_null()).then(pl.lit('')) \
                .when(pl.col('pt') == pl.col('pt').floor()).then(pl.col('pt').cast(pl.Int64).cast(pl.Utf8)) \
                .otherwise(pl.col('pt').cast(pl.Utf8))
    write_tsv(reads.select('read_id', 'contig', 'compartment', pt_text.alias('pt'), 'status', 'isoform', 'gene'),
              f"{p}.polya_reads.tsv.gz")

    tx_gene, iso_info, gene_names = annotations(args.catalogue, args.classification)
    where = frame(spans(args.catalogue), 'feature_id', ['chrom', 'start', 'end', 'strand'])
    iso_df = frame(iso_info, 'feature_id', ['iso_name', 'iso_tx', 'iso_cat'])
    txg_df = frame({k: (v,) for k, v in tx_gene.items()}, 'feature_id', ['tx_gene'])
    gn_df = frame({k: (v or '.',) for k, v in gene_names.items()}, 'feature_id', ['gene_name'])
    tables = []
    for level in ('isoform', 'gene'):
        d = describe(used.filter(pl.col(level) != '.'), level).rename({level: 'feature_id'}) \
              .with_columns(pl.lit(level).alias('level'))
        if level == 'isoform':
            d = d.join(iso_df, on='feature_id', how='left').join(txg_df, on='feature_id', how='left').with_columns(
                pl.col('tx_gene').fill_null('.').alias('gene_id'), pl.col('iso_name').fill_null('.').alias('ref_gene_name'),
                pl.col('iso_tx').fill_null('.').alias('ref_tx_id'), pl.col('iso_cat').fill_null('.').alias('structural_category'))
        else:
            d = d.join(gn_df, on='feature_id', how='left').with_columns(
                pl.col('feature_id').alias('gene_id'), pl.col('gene_name').fill_null('.').alias('ref_gene_name'),
                pl.lit('.').alias('ref_tx_id'), pl.lit('.').alias('structural_category'))
        d = d.join(where, on='feature_id', how='left').with_columns(
            pl.when(pl.col('chrom') == args.mito_contig).then(pl.lit('mito')).otherwise(pl.lit('nuclear')).alias('compartment'))
        tables.append(d.select(FEATURE_COLS))
    features = pl.concat(tables, how='vertical')
    unplaced = features.filter(pl.col('start').is_null())
    if unplaced.height:
        sys.exit(f"polya_features.py: {unplaced.height} feature(s) not in {args.catalogue}, e.g. {unplaced['feature_id'].head(3).to_list()}")
    features = features.sort(['chrom', 'start', 'level', 'feature_id'])
    write_tsv(features, f"{p}.polya_features.tsv.gz", float_precision=2)

    for level in ('isoform', 'gene'):
        sub = features.filter((pl.col('level') == level) & (pl.col('n_tails') >= args.min_reads))
        with open(f"{p}.{level}_polya.bed", 'w') as fh:
            fh.write('#' + '\t'.join(BED_COLS) + '\n')
            for r in sub.iter_rows(named=True):
                vals = [r['chrom'], r['start'], r['end'], r['feature_id'], min(int(r['n_tails']), 1000), r['strand'],
                        r['feature_id'], r['gene_id'], r['ref_gene_name'], r['ref_tx_id'], r['structural_category'],
                        r['compartment'], int(r['n_tails'])]
                vals += [f"{r[k]:.2f}" for k in ('median_pt', 'mean_pt', 'q25_pt', 'q75_pt', 'min_pt', 'max_pt')]
                fh.write('\t'.join(map(str, vals)) + '\n')

    counts = dict(reads.group_by('status').agg(pl.len()).iter_rows())
    in_bed = lambda lv: features.filter((pl.col('level') == lv) & (pl.col('n_tails') >= args.min_reads)).height  # noqa: E731
    print(f"polya_features.py: {reads.height} reads; used {counts.get('used', 0)}, not found {counts.get('not_found', 0)}, "
          f"missing {counts.get('missing', 0)}, censored {counts.get('censored', 0)}; "
          f"isoforms in BED {in_bed('isoform')}, genes {in_bed('gene')}", file=sys.stderr)


if __name__ == '__main__':
    main()
