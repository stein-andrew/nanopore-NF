"""Shared polars table I/O for the bin/ scripts, so they all read and write TSVs one way.

    read_tsv(path, names)          headerless TSV (plain or .gz) -> DataFrame of strings named
                                   `names`; empty file -> empty frame; wrong width -> exit
    read_tsv_header(path)          TSV with a header line (plain or .gz) -> DataFrame of strings
    write_tsv(df, path, header=True, float_precision=None)
                                   .gz when the path ends so; never quoted; nulls written empty
    codes_of(series)               an integer code per value and the distinct values, in
                                   first-seen order (what pandas.factorize gave)

Everything is read as text (infer_schema_length=0) and cast by the caller, so no column changes
type because of what polars guesses from the first rows. Empty fields read as null.
"""
import gzip
import io
import os
import sys

import numpy as np
import polars as pl


def _bytes(path):
    with (gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')) as fh:
        return fh.read()


def _empty(names):
    return pl.DataFrame({n: pl.Series(n, [], dtype=pl.Utf8) for n in names})


def read_tsv(path, names):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return _empty(names)
    data = _bytes(path)
    if not data.strip():
        return _empty(names)
    df = pl.read_csv(io.BytesIO(data), separator='\t', has_header=False, infer_schema_length=0, quote_char=None)
    if df.width != len(names):
        sys.exit(f"{os.path.basename(sys.argv[0])}: {path} has {df.width} columns, expected {len(names)} "
                 f"({' '.join(names)})")
    return df.rename(dict(zip(df.columns, names)))


def read_tsv_header(path):
    return pl.read_csv(io.BytesIO(_bytes(path)), separator='\t', has_header=True, infer_schema_length=0, quote_char=None)


def write_tsv(df, path, header=True, float_precision=None):
    text = df.write_csv(separator='\t', include_header=header, quote_style='never', float_precision=float_precision)
    if path.endswith('.gz'):
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        with open(path, 'w') as fh:
            fh.write(text)


def codes_of(series):
    values = series.to_list()
    names = list(dict.fromkeys(values))
    index = {v: i for i, v in enumerate(names)}
    return np.fromiter((index[v] for v in values), dtype=np.int64, count=len(values)), names
