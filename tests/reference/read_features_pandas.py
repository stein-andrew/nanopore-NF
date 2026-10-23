#!/usr/bin/env python3
"""Which catalogue isoform and gene each read belongs to, from IsoQuant's read assignments.

    read_features.py S.read_assignments.tsv.gz > S.read_features.tsv

Output, one row per read that has either, no header: read_id <tab> isoform <tab> gene,
with '.' for none.
  isoform  the read's single assignment when its type is unique or unique_minor_difference,
           i.e. the reads IsoQuant's unique_only transcript counts are made of;
  gene     the one gene_id across all the read's rows ('.' if none or several), so reads
           ambiguous between isoforms of the same gene still count for that gene.

IsoQuant writes one row per candidate isoform for ambiguous reads; the header line starts
with '#'. Columns are found by name (read_id, isoform_id, gene_id, assignment_type).
"""
import gzip
import sys

import pandas as pd

UNIQUE = {'unique', 'unique_minor_difference'}
NEEDED = ['read_id', 'isoform_id', 'gene_id', 'assignment_type']


def header(path):
    """(names, number of lines before the data) from the '#read_id ...' header."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for n, line in enumerate(fh):
            if not line.startswith('#'):
                break
            names = line.lstrip('#').rstrip('\n').split('\t')
            if all(c in names for c in NEEDED):
                return names, n + 1
    sys.exit(f"read_features.py: {path}: no '#' header line with {', '.join(NEEDED)}")


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <read_assignments.tsv[.gz]> > read_features.tsv")
    path = sys.argv[1]
    names, skip = header(path)
    cols = [names.index(c) for c in NEEDED]
    ra = pd.read_csv(path, sep='\t', header=None, skiprows=skip, usecols=cols, dtype=str,
                     keep_default_na=False, na_filter=False, comment=None)
    ra.columns = [names[i] for i in sorted(cols)]

    # vectorized: a read keeps its isoform only if it has exactly one unique-type row, and its
    # gene only if all its rows name the same gene
    uniq = ra[ra['assignment_type'].isin(UNIQUE) & (ra['isoform_id'] != '.')]
    once = uniq[uniq['read_id'].map(uniq['read_id'].value_counts()) == 1].set_index('read_id')['isoform_id']

    g = ra[ra['gene_id'] != '.'].drop_duplicates(['read_id', 'gene_id'])
    genes = g[g['read_id'].map(g['read_id'].value_counts()) == 1].set_index('read_id')['gene_id']

    out = pd.concat([once.rename('isoform'), genes.rename('gene')], axis=1).fillna('.')
    out = out[(out['isoform'] != '.') | (out['gene'] != '.')]
    out.to_csv(sys.stdout, sep='\t', header=False)
    print(f"read_features.py: {ra['read_id'].nunique()} reads; isoform for {(out['isoform'] != '.').sum()}, "
          f"gene for {(out['gene'] != '.').sum()}", file=sys.stderr)


if __name__ == '__main__':
    main()
