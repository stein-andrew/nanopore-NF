#!/usr/bin/env python3
"""Reference (pandas) implementation of bin/cohort_matrices.py, for tests/compare_polars_port.sh.
Same interface and outputs; see bin/cohort_matrices.py for what each subcommand does."""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

from feature_mods_pandas import annotations

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


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def write(df, path, fmt=None):
    df.to_csv(path, sep='\t', index=False, float_format=fmt, na_rep='',
              compression='gzip' if path.endswith('.gz') else None)


def gene_counts(path):
    counts = {}
    with opener(path) as fh:
        fh.readline()
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
    m = pd.DataFrame({'gene_id': genes})
    m['gene_name'] = [names.get(g, '.') or '.' for g in genes]
    cpm = m.copy()
    for s, _ in samples:
        col = np.array([counts[s].get(g, 0.0) for g in genes])
        total = col.sum()
        if total <= 0:
            sys.exit(f"cohort_matrices.py: {s}: no assigned reads, so CPM is undefined")
        m[s] = col.astype(np.int64) if np.all(col == np.floor(col)) else col
        cpm[s] = col / total * 1e6
    write(m, f"{args.prefix}.gene_counts.tsv.gz")
    write(cpm, f"{args.prefix}.gene_cpm.tsv.gz", '%.4f')
    print(f"cohort_matrices.py genes: {len(genes)} genes x {len(samples)} samples", file=sys.stderr)


def cmd_mods(args):
    sites = pd.read_csv(args.sites, sep='\t', header=None, names=SITE, dtype=str, keep_default_na=False) \
        if os.path.getsize(args.sites) else pd.DataFrame(columns=SITE)
    names = SITE + ['sample', 'n_valid', 'n_mod', 'pct', 'pct_fpcorr']
    long = pd.read_csv(args.long, sep='\t', header=None, names=names, dtype=str, keep_default_na=False) \
        if os.path.getsize(args.long) else pd.DataFrame(columns=names)
    for d in (sites, long):
        d['start'] = d['start'].astype(np.int64)
        d['end'] = d['end'].astype(np.int64)
    long['n_valid'] = long['n_valid'].astype(np.int64)
    long['n_mod'] = long['n_mod'].astype(np.int64)
    long['pct_fpcorr'] = pd.to_numeric(long['pct_fpcorr'], errors='coerce')
    with open(f"{args.prefix}.mods.summary.tsv", 'w') as out:
        out.write("sample\tfinal_sites\tmeasured_sites\n")
        for s in args.samples:
            sub = long[long['sample'] == s]
            out.write(f"{s}\t{int((sub['pct_fpcorr'] >= args.min_pct).sum())}\t{len(sub)}\n")
    for metric, fmt in (('n_valid', None), ('n_mod', None), ('pct_fpcorr', '%.2f')):
        wide = sites.copy()
        for s in args.samples:
            part = long.loc[long['sample'] == s, SITE + [metric]].rename(columns={metric: s})
            wide = wide.merge(part, on=SITE, how='left')
            if metric != 'pct_fpcorr':
                wide[s] = wide[s].astype('Int64')
        write(wide, f"{args.prefix}.mods.{metric}.tsv.gz", fmt)
    print(f"cohort_matrices.py mods: {len(sites)} sites x {len(args.samples)} samples", file=sys.stderr)


def cmd_polya(args):
    samples = pairs(args.features, '--features')
    tabs = {s: pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False) for s, p in samples}
    for level in ('isoform', 'gene'):
        ann_cols = ['feature_id', 'gene_id', 'ref_gene_name', 'chrom', 'start', 'end', 'strand']
        ann = pd.concat([t.loc[t['level'] == level, ann_cols] for t in tabs.values()]).drop_duplicates('feature_id')
        ann['start'] = ann['start'].astype(np.int64)
        ann = ann.sort_values(['chrom', 'start', 'feature_id']).reset_index(drop=True)
        for metric, fmt in (('median_pt', '%.2f'), ('n_tails', None)):
            wide = ann.copy()
            for s, t in tabs.items():
                part = t.loc[t['level'] == level, ['feature_id', metric]].rename(columns={metric: s})
                part[s] = pd.to_numeric(part[s])
                wide = wide.merge(part, on='feature_id', how='left')
                if metric == 'n_tails':                     # a feature without tails in a sample has 0 there
                    wide[s] = wide[s].fillna(0).astype(np.int64)
            write(wide, f"{args.prefix}.polya.{level}.{metric}.tsv.gz", fmt)
    print(f"cohort_matrices.py polya: {len(samples)} samples", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
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
    m.add_argument('--min-pct', type=float, required=True)
    m.add_argument('--prefix', required=True)
    p = sub.add_parser('polya')
    p.add_argument('--features', action='append', required=True, metavar='SAMPLE=PATH')
    p.add_argument('--prefix', required=True)
    args = ap.parse_args()
    {'genes': cmd_genes, 'mods': cmd_mods, 'polya': cmd_polya}[args.cmd](args)


if __name__ == '__main__':
    main()
