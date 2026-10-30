#!/usr/bin/env python3
"""A first look for batch effects in expression: the sample-by-sample Pearson correlation of
log2(CPM + 1) over genes whose cohort median CPM is above 1, and a batch-by-batch table of
mean correlations: within a batch on the diagonal, between two batches off it. Each pair of
batches gets its own value, so one batch that differs shows as low values in its own row and
column only (a single "within vs all the others" number per batch would drag every batch
down with it). Off-diagonal values clearly below the diagonal (by more than ~0.01) mark a
batch that differs systematically: it belongs in every between-sample comparison as a
covariate.

    batch_correlation.py <cohort.<ref>.gene_cpm.tsv.gz> <analysis_samplesheet.csv> [--column batch]

The samplesheet decides which columns are samples (its sample_id) and their batch (--column).
polars for the tables, numpy for the arithmetic.
"""
import argparse
import gzip
import sys

import numpy as np
import polars as pl


def read_table(path, sep):
    """Everything as text first, so no column's type depends on its first rows."""
    data = gzip.open(path).read() if path.endswith('.gz') else open(path, 'rb').read()
    return pl.read_csv(data, separator=sep, infer_schema_length=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cpm')
    ap.add_argument('samplesheet')
    ap.add_argument('--column', default='batch', help='samplesheet column holding the batch (default: batch)')
    ap.add_argument('--min-median-cpm', type=float, default=1.0)
    args = ap.parse_args()

    sheet = read_table(args.samplesheet, ',')
    if args.column not in sheet.columns:
        sys.exit(f'batch_correlation.py: no column {args.column!r} in {args.samplesheet}')
    batch = dict(zip(sheet['sample_id'].to_list(), sheet[args.column].to_list()))

    cpm = read_table(args.cpm, '\t')
    samples = [s for s in sheet['sample_id'].to_list() if s in cpm.columns]
    missing = [s for s in sheet['sample_id'].to_list() if s not in cpm.columns]
    if missing:
        print(f'warning: in the samplesheet but not the CPM table: {" ".join(missing)}', file=sys.stderr)
    if len(samples) < 2:
        sys.exit('batch_correlation.py: fewer than two samples to correlate')

    values = cpm.select([pl.col(s).cast(pl.Float64, strict=False) for s in samples])
    n_null = sum(values.null_count().row(0))
    if n_null:
        print(f'warning: {n_null} empty or non-numeric CPM values, counted as 0', file=sys.stderr)
    m = values.fill_null(0.0).to_numpy()

    keep = np.median(m, axis=1) > args.min_median_cpm
    x = np.log2(m[keep] + 1.0)
    r = np.corrcoef(x, rowvar=False)
    print(f'{int(keep.sum())} of {len(m)} genes with median CPM > {args.min_median_cpm:g}, {len(samples)} samples')

    # the table: one row per sample, labelled with its batch
    label = [f'{s} {batch[s]}' for s in samples]
    width = max(len(t) for t in label)
    print(' ' * width + ' ' + ' '.join(f'{s:>7}' for s in samples))
    for i, t in enumerate(label):
        print(f'{t:<{width}} ' + ' '.join(f'{v:7.3f}' for v in r[i]))

    # batch x batch: mean r within a batch (diagonal; n/a with one sample) and between batches
    groups = np.array([batch[s] for s in samples])
    names = sorted(set(groups))
    idx = {b: np.where(groups == b)[0] for b in names}
    def mean_r(a, b):
        if a == b:
            n = len(idx[a])
            return r[np.ix_(idx[a], idx[a])][np.triu_indices(n, 1)].mean() if n > 1 else float('nan')
        return r[np.ix_(idx[a], idx[b])].mean()
    w = max(len(b) for b in names) + 6
    print('\nmean r by batch (n samples):')
    print(' ' * w + ' '.join(f'{b:>8}' for b in names))
    for a in names:
        print(f'{a + f" ({len(idx[a])})":<{w}}' + ' '.join(f'{mean_r(a, b):8.3f}' if not np.isnan(mean_r(a, b)) else f'{"n/a":>8}' for b in names))

if __name__ == '__main__':
    main()
