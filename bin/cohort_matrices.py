#!/usr/bin/env python3
"""Cohort matrices for comparing samples: one row per feature, one column per sample.

    cohort_matrices.py genes --counts S1=S1.gene_counts.tsv ... --catalogue C.gtf --classification C.txt.gz --prefix P
    cohort_matrices.py mods  --sites final_sites.tsv --long site_values.tsv --samples S1 S2 ... --min-pct 20 --prefix P
    cohort_matrices.py polya --features S1=S1.polya_features.tsv.gz ... --prefix P

genes  IsoQuant requantification gene counts (catalogue genes; __ summary rows excluded) ->
       <P>.gene_counts.tsv.gz and <P>.gene_cpm.tsv.gz (CPM over each sample's assigned reads,
       as filter_models.py). A gene absent from a sample's counts is 0 there.
mods   every site that passes the final stage (corrected percent >= --min-pct, coverage >=
       the stage minimum) in at least one sample, with each sample's values there ->
       <P>.mods.{n_valid,n_mod,pct_fpcorr}.tsv.gz and <P>.mods.summary.tsv. EMPTY means not
       measured in that sample (below the coverage minimum), not unmodified; an empty
       pct_fpcorr with an n_valid means gIVT had no correction for that site.
       Inputs come from the module's awk steps: --sites (chrom start end mod strand) and
       --long (chrom start end mod strand sample n_valid n_mod pct pct_fpcorr).
polya  per-sample polya_features.tsv.gz -> <P>.polya.{isoform,gene}.{median_pt,n_tails}.tsv.gz.
       A feature with no usable tails in a sample: n_tails 0, median_pt empty.
"""
import argparse
import gzip
import sys

import polars as pl

from feature_mods import annotations
from tabular import read_tsv, read_tsv_header, write_tsv

SITE = ['chrom', 'start', 'end', 'mod', 'strand']


def pairs(items, flag):
    out = []
    for item in items:
        s, sep, p = item.partition('=')
        if not sep or not s or not p:
            sys.exit(f"cohort_matrices.py: {flag} expects SAMPLE=PATH, got {item}")
        if s in (x for x, _ in out):
            sys.exit(f"cohort_matrices.py: {flag}: sample {s} given twice")
        out.append((s, p))
    return out


def gene_counts(path):
    counts = {}
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        fh.readline()                                   # feature_id <tab> count
        for line in fh:
            if not line.strip():
                continue
            feature, count = line.rstrip('\n').split('\t')[:2]
            if not feature.startswith('__'):
                counts[feature] = float(count)
    return counts


def cmd_genes(args):
    samples = pairs(args.counts, '--counts')
    counts = {s: gene_counts(p) for s, p in samples}
    genes = sorted(set().union(*[set(c) for c in counts.values()]))
    _, _, names = annotations(args.catalogue, args.classification)
    base = {'gene_id': genes, 'gene_name': [names.get(g, '.') or '.' for g in genes]}
    m, cpm = dict(base), dict(base)
    for s, _ in samples:
        col = [counts[s].get(g, 0.0) for g in genes]
        total = sum(col)
        if total <= 0:
            sys.exit(f"cohort_matrices.py: {s}: no assigned reads, so CPM is undefined")
        m[s] = pl.Series(s, col, dtype=pl.Float64)
        if all(v == int(v) for v in col):
            m[s] = m[s].cast(pl.Int64)
        cpm[s] = pl.Series(s, [v / total * 1e6 for v in col], dtype=pl.Float64)
    write_tsv(pl.DataFrame(m), f"{args.prefix}.gene_counts.tsv.gz")
    write_tsv(pl.DataFrame(cpm), f"{args.prefix}.gene_cpm.tsv.gz", float_precision=4)
    print(f"cohort_matrices.py genes: {len(genes)} genes x {len(samples)} samples", file=sys.stderr)


def cmd_mods(args):
    sites = read_tsv(args.sites, SITE)
    long = read_tsv(args.long, SITE + ['sample', 'n_valid', 'n_mod', 'pct', 'pct_fpcorr'])
    ints = [pl.col('start').cast(pl.Int64), pl.col('end').cast(pl.Int64)]
    sites = sites.with_columns(ints).with_row_index('_i')
    long = long.with_columns(ints + [pl.col('n_valid').cast(pl.Int64), pl.col('n_mod').cast(pl.Int64),
                                     pl.col('pct_fpcorr').cast(pl.Float64, strict=False)])
    with open(f"{args.prefix}.mods.summary.tsv", 'w') as out:
        out.write("sample\tfinal_sites\tmeasured_sites\n")
        for s in args.samples:
            sub = long.filter(pl.col('sample') == s)
            final = sub.filter(pl.col('pct_fpcorr') >= args.min_pct).height
            out.write(f"{s}\t{final}\t{sub.height}\n")
    for metric, precision in (('n_valid', None), ('n_mod', None), ('pct_fpcorr', 2)):
        wide = sites
        for s in args.samples:
            part = long.filter(pl.col('sample') == s).select(SITE + [metric]).rename({metric: s})
            wide = wide.join(part, on=SITE, how='left')
        write_tsv(wide.sort('_i').drop('_i'), f"{args.prefix}.mods.{metric}.tsv.gz", float_precision=precision)
    print(f"cohort_matrices.py mods: {sites.height} sites x {len(args.samples)} samples", file=sys.stderr)


def cmd_polya(args):
    samples = pairs(args.features, '--features')
    tabs = {s: read_tsv_header(p) for s, p in samples}
    ann_cols = ['feature_id', 'gene_id', 'ref_gene_name', 'chrom', 'start', 'end', 'strand']
    for level in ('isoform', 'gene'):
        ann = pl.concat([t.filter(pl.col('level') == level).select(ann_cols) for t in tabs.values()], how='vertical') \
                .unique(subset=['feature_id'], keep='first', maintain_order=True) \
                .with_columns(pl.col('start').cast(pl.Int64)) \
                .sort(['chrom', 'start', 'feature_id']).with_row_index('_i')
        for metric, precision in (('median_pt', 2), ('n_tails', None)):
            wide = ann
            for s, t in tabs.items():
                part = t.filter(pl.col('level') == level).select(
                    'feature_id', pl.col(metric).cast(pl.Float64).alias(s))
                wide = wide.join(part, on='feature_id', how='left')
                if metric == 'n_tails':                 # a feature without tails in a sample has 0 there
                    wide = wide.with_columns(pl.col(s).fill_null(0).cast(pl.Int64))
            write_tsv(wide.sort('_i').drop('_i'), f"{args.prefix}.polya.{level}.{metric}.tsv.gz", float_precision=precision)
    print(f"cohort_matrices.py polya: {len(samples)} samples", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('genes')
    g.add_argument('--counts', action='append', required=True, metavar='SAMPLE=PATH')
    g.add_argument('--catalogue', required=True)
    g.add_argument('--classification', required=True)
    g.add_argument('--prefix', required=True)
    m = sub.add_parser('mods')
    m.add_argument('--sites', required=True)
    m.add_argument('--long', required=True)
    m.add_argument('--samples', nargs='+', required=True)
    m.add_argument('--min-pct', type=float, required=True, help='final stage: corrected percent >= this')
    m.add_argument('--prefix', required=True)
    p = sub.add_parser('polya')
    p.add_argument('--features', action='append', required=True, metavar='SAMPLE=PATH')
    p.add_argument('--prefix', required=True)
    args = ap.parse_args()
    {'genes': cmd_genes, 'mods': cmd_mods, 'polya': cmd_polya}[args.cmd](args)


if __name__ == '__main__':
    main()
