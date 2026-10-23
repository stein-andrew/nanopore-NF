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
with '#'. Columns are found by name (read_id, isoform_id, gene_id, assignment_type) and taken
by position. Output is sorted by read_id.
"""

import gzip
import io
import sys

import polars as pl

UNIQUE = ['unique', 'unique_minor_difference']
NEEDED = ['read_id', 'isoform_id', 'gene_id', 'assignment_type']


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <read_assignments.tsv[.gz]> > read_features.tsv")
    path = sys.argv[1]
    with (gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')) as fh:
        data = fh.read()

    # the '#' header that names the columns, and where the data starts
    names, start, offset = None, 0, 0
    for line in io.BytesIO(data):
        if not line.startswith(b'#'):
            break
        offset += len(line)
        fields = line.decode().lstrip('#').rstrip('\n').split('\t')
        if all(c in fields for c in NEEDED):
            names, start = fields, offset
    if names is None:
        sys.exit(f"read_features.py: {path}: no '#' header line with {', '.join(NEEDED)}")
    body = data[start:]
    if not body.strip():
        ra = pl.DataFrame({c: pl.Series(c, [], dtype=pl.Utf8) for c in NEEDED})
    else:
        df = pl.read_csv(io.BytesIO(body), separator='\t', has_header=False, infer_schema_length=0, quote_char=None,
                         truncate_ragged_lines=True)       # IsoQuant rows carry a field more than its header names
        ra = pl.DataFrame({c: df.to_series(names.index(c)) for c in NEEDED})

    # a read keeps its isoform only if it has exactly one unique-type row, and its gene only if
    # all its rows name the same gene
    once = ra.filter(pl.col('assignment_type').is_in(UNIQUE) & (pl.col('isoform_id') != '.')) \
             .group_by('read_id').agg(pl.len().alias('n'), pl.col('isoform_id').first().alias('isoform')) \
             .filter(pl.col('n') == 1).select('read_id', 'isoform')
    genes = ra.filter(pl.col('gene_id') != '.').unique(subset=['read_id', 'gene_id']) \
              .group_by('read_id').agg(pl.len().alias('n'), pl.col('gene_id').first().alias('gene')) \
              .filter(pl.col('n') == 1).select('read_id', 'gene')
    out = once.join(genes, on='read_id', how='full', coalesce=True) \
              .with_columns(pl.col('isoform').fill_null('.'), pl.col('gene').fill_null('.')) \
              .sort('read_id')
    sys.stdout.write(out.write_csv(separator='\t', include_header=False, quote_style='never'))
    print(f"read_features.py: {ra['read_id'].n_unique()} reads; isoform for {(out['isoform'] != '.').sum()}, "
          f"gene for {(out['gene'] != '.').sum()}", file=sys.stderr)


if __name__ == '__main__':
    main()
