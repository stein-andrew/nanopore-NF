#!/usr/bin/env python3
"""Vectorized drop-in for gIVT's IVT_fp_correction.py: same arguments, same corrected values.

    fp_correct.py -m S.cov20.bed -ref genome.fa -e errortable.tsv.gz \\
        -mt a,0.81 m,0.80 ... -o S.cov20_fpcorr.bed \\
        [--valid_only | --filter_kmer --filter_mismatch] [-n MIN_COV] [--engine polars|pandas]

Same as IVT_fp_correction.py (checked row for row by tests/compare_fpcorr.py):
  * column 19 = max(0, percent_modified - 100 * fp), fp from the error table for the site's
    9mer (reverse-complemented on '-') and mod code, in the column gIVT's mod_to_col()
    picks for that code's threshold. That function is copied as is, including that a
    threshold >= 0.99 falls back to the 0.70 column; the pipeline refuses those first;
  * which sites get no correction (empty column 19), or are dropped by the --filter_*
    flags, and the order the checks apply in: no 9mer (contig not in the reference, or
    within 4 bp of an end) -> centre base is not the mod's base (after reverse complement;
    case-sensitive, so a soft-masked lowercase base is a mismatch, as in gIVT) -> 9mer
    and mod code not in the table;
  * row order, and --min_coverage on column 10.

Different on purpose:
  * columns 1-18 are each input line's own text (gIVT round-trips them through pandas,
    which rewrites e.g. 50.00 as 50.0);
  * no stdout line per uncorrectable site, one count per reason instead;
  * numeric contig names (1, 22) work; gIVT's pandas read makes them ints that never match;
  * the error table may only hold uppercase ACGT 9mers (gIVT's has 65,536 per mod code),
    because sites are looked up by a 2-bit encoding; anything else is refused;
  * the reference must be an uncompressed FASTA with its .fai next to it (the pipeline's is).

Why it is fast:
  * the reference is never loaded: the FASTA is memory-mapped and each site's 9 bases are
    read at offsets computed from the .fai (gIVT holds the whole genome as Python strings);
  * no per-site Python loop: 9mers are gathered for all sites at once, encoded as 18-bit
    integers and looked up in a dense (mod code x 4^9) array;
  * only the six columns the correction needs are parsed, by polars (default; threads
    from POLARS_MAX_THREADS) or pandas; each output line is its input line plus column 19.
"""
import argparse
import csv
import gzip
import os
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd

MOD_TO_BASE = {'17802': 'T', '19229': 'G', '19228': 'C', '19227': 'T',
               '69426': 'A', 'a': 'A', 'm': 'C', '17596': 'A'}          # gIVT's mapping
NCOL = 18                                                                # modkit bedMethyl
OFFSETS = np.arange(-4, 5)
SHIFTS = (2 * np.arange(8, -1, -1)).astype(np.int64)                    # first base most significant

COMPLEMENT = np.arange(256, dtype=np.uint8)                              # only ACGT change, as str.translate
for _a, _b in zip(b'ACGT', b'TGCA'):
    COMPLEMENT[_a] = _b
BASE_CODE = np.full(256, 255, dtype=np.uint8)                            # 255 = not A/C/G/T
for _i, _b in enumerate(b'ACGT'):
    BASE_CODE[_b] = _i


def log(msg):
    print(f"fp_correct.py: {msg}", file=sys.stderr, flush=True)


def mod_to_col(mod_code_to_threshold, line):
    """gIVT's column choice, unchanged: the column before the first threshold > the mod's."""
    m_to_c = defaultdict(lambda: 2)
    for mod in mod_code_to_threshold:
        for i, l in enumerate(line.strip().split('\t')):
            if "false_positive" not in l:
                continue
            if float(l.split('_')[2]) > mod_code_to_threshold[mod]:
                if i == 2:
                    m_to_c[mod] = 2
                else:
                    m_to_c[mod] = i - 1
                break
    return m_to_c


def encode_kmers(kmers):
    """(n, 9) uint8 bases -> 18-bit codes, and which rows are all A/C/G/T."""
    v = BASE_CODE[kmers]
    return (v.astype(np.int64) << SHIFTS).sum(axis=1), (v != 255).all(axis=1)


def load_table(path, thresholds):
    """mod code -> row, plus (n_mods, 4^9) false-positive rates and a presence mask.

    Presence is kept apart from the value so a NaN rate behaves as in gIVT
    (max(0, pct - NaN) is 0 there), not as "not in table". round_trip parsing gives the
    same doubles as gIVT's float(); a repeated (mod, 9mer) keeps the last row, as in gIVT."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        header = fh.readline()
    mod_col = mod_to_col(thresholds, header)
    # Only the columns some mod code will read: those mod_to_col chose, plus column 2,
    # which every code without a threshold reads. round_trip parsing is exact but slow,
    # so all 30 threshold columns would cost seconds for nothing.
    used = sorted({2, *mod_col.values()})
    tab = pd.read_csv(path, sep='\t', header=None, skiprows=1, usecols=[0, 1, *used], dtype={0: str, 1: str},
                      keep_default_na=False, na_values=['nan', 'NaN'], float_precision='round_trip')
    mers, mods = tab[0].to_numpy(dtype=object), tab[1].to_numpy(dtype=object)
    lengths = tab[0].str.len().to_numpy()
    raw = np.asarray(mers, dtype='S9').view(np.uint8).reshape(-1, 9)   # longer 9mers are caught by `lengths`
    codes, acgt = encode_kmers(raw)
    bad = np.flatnonzero((lengths != 9) | ~acgt)
    if len(bad):
        sys.exit(f"fp_correct.py: {path} row {bad[0] + 2}: 9mer {mers[bad[0]]!r} is not 9 uppercase A/C/G/T")

    mod_codes, names = pd.factorize(mods)
    index = {str(m): i for i, m in enumerate(names)}
    rate = np.empty(len(tab))
    for i, m in enumerate(names):
        rows = mod_codes == i
        rate[rows] = tab[mod_col[str(m)]].to_numpy(dtype=np.float64)[rows]
    last = ~pd.DataFrame({'m': mod_codes, 'k': codes}).duplicated(keep='last').to_numpy()

    fps = np.zeros((len(names), 4 ** 9))
    present = np.zeros((len(names), 4 ** 9), dtype=bool)
    fps[mod_codes[last], codes[last]] = rate[last]
    present[mod_codes[last], codes[last]] = True
    return index, fps, present


def read_fai(fasta):
    if fasta.endswith(('.gz', '.bgz')):
        sys.exit(f"fp_correct.py: {fasta}: needs an uncompressed FASTA (it is memory-mapped)")
    fai = fasta + '.fai'
    if not os.path.exists(fai):
        sys.exit(f"fp_correct.py: no {fai} (samtools faidx {fasta})")
    out = {}
    with open(fai) as fh:
        for line in fh:
            name, length, offset, linebases, linewidth = line.split('\t')[:5]
            out[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return out


def read_sites(path, engine):
    """Raw lines (written back unchanged) and the six columns the correction needs."""
    with open(path, 'rb') as fh:
        data = fh.read()
    nlines = data.count(b'\n') + (0 if data.endswith(b'\n') else 1)
    if data.count(b'\t') != (NCOL - 1) * nlines:
        sys.exit(f"fp_correct.py: {path}: not every line has {NCOL} tab-separated columns (modkit bedMethyl)")
    lines = data.decode('ascii').split('\n')
    del data
    if lines[-1] == '':
        lines.pop()
    if engine == 'polars':
        import polars as pl
        # All 18 columns as text, taken by position: the names polars generates for a
        # headerless file with a column selection differ between versions.
        df = pl.read_csv(path, separator='\t', has_header=False, infer_schema_length=0, quote_char=None)
        if df.width != NCOL:
            sys.exit(f"fp_correct.py: {path} has {df.width} columns, expected {NCOL} (modkit bedMethyl)")
        col = df.to_series
        chrom_codes, chroms = pd.factorize(col(0).to_numpy())
        mod_codes, mods = pd.factorize(col(3).to_numpy())
        sites = dict(minus=(col(5) == '-').to_numpy(),
                     start=col(1).cast(pl.Int64).to_numpy(),
                     cov=col(9).cast(pl.Int64).to_numpy(),
                     pct=col(10).cast(pl.Float64).to_numpy())
    else:
        df = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 3, 5, 9, 10], quoting=csv.QUOTE_NONE,
                         keep_default_na=False, na_filter=False, float_precision='round_trip',
                         dtype={0: 'category', 1: np.int64, 3: 'category', 5: 'category', 9: np.int64, 10: np.float64})
        chrom_codes, chroms = df[0].cat.codes.to_numpy(), list(df[0].cat.categories)
        mod_codes, mods = df[3].cat.codes.to_numpy(), list(df[3].cat.categories)
        sites = dict(minus=(df[5] == '-').to_numpy(), start=df[1].to_numpy(),
                     cov=df[9].to_numpy(), pct=df[10].to_numpy())
    sites.update(chrom_codes=np.asarray(chrom_codes, dtype=np.int64), chroms=[str(c) for c in chroms],
                 mod_codes=np.asarray(mod_codes, dtype=np.int64), mods=[str(m) for m in mods])
    if len(sites['start']) != len(lines):
        sys.exit(f"fp_correct.py: {path}: parsed {len(sites['start'])} rows but the file has {len(lines)} lines")
    return lines, sites


def write_sites(lines, keep, adjusted, path, chunk=2_000_000):
    """Each kept input line plus column 19 (empty where uncorrected), in input order."""
    rows = np.flatnonzero(keep)
    with open(path, 'w') as out:
        for i in range(0, len(rows), chunk):
            part = rows[i:i + chunk].tolist()
            col19 = ['' if v != v else repr(v) for v in adjusted[part].tolist()]
            kept = [lines[r] for r in part]
            out.write('\n'.join(map('\t'.join, zip(kept, col19))))
            out.write('\n')


def gather_kmers(codes, names, start, fasta):
    """(n, 9) uint8 forward 9mers read straight from the memory-mapped FASTA, and whether
    each site has one at all (gIVT: seq[pos-4:pos+5] must exist and be 9 long)."""
    fai = read_fai(fasta)
    mm = np.memmap(fasta, dtype=np.uint8, mode='r')
    n = len(start)
    kmers = np.zeros((n, 9), dtype=np.uint8)
    have = np.zeros(n, dtype=bool)
    order = np.argsort(codes, kind='stable')
    bounds = np.searchsorted(codes[order], np.arange(len(names) + 1))
    for ci, name in enumerate(names):
        if name not in fai:
            continue
        length, offset, linebases, linewidth = fai[name]
        rows = order[bounds[ci]:bounds[ci + 1]]
        p = start[rows]
        ok = (p >= 4) & (p + 5 <= length)
        q = p[ok, None] + OFFSETS
        kmers[rows[ok]] = mm[offset + (q // linebases) * linewidth + q % linebases]
        have[rows[ok]] = True
    return kmers, have


def correct(sites, fasta, index, fps, present):
    """Corrected percent per site (NaN where none) and a reason code per site:
    0 corrected, 1 no 9mer, 2 centre base mismatch, 3 9mer/mod not in table."""
    n = len(sites['start'])
    kmers, have = gather_kmers(sites['chrom_codes'], sites['chroms'], sites['start'], fasta)

    minus = sites['minus']
    kmers[minus] = COMPLEMENT[kmers[minus]][:, ::-1]

    mcodes, mods = sites['mod_codes'], sites['mods']
    expected = np.array([ord(MOD_TO_BASE.get(m, '.')) for m in mods], dtype=np.uint8)[mcodes]
    centre_ok = have & (kmers[:, 4] == expected)

    code, acgt = encode_kmers(kmers)
    trow = np.array([index.get(m, -1) for m in mods], dtype=np.int64)[mcodes]

    look = centre_ok & acgt & (trow >= 0)
    found = np.zeros(n, dtype=bool)
    found[look] = present[trow[look], code[look]]

    adjusted = np.full(n, np.nan)
    x = sites['pct'][found] - 100 * fps[trow[found], code[found]]
    adjusted[found] = np.where(x > 0, x, 0.0)          # gIVT's max(0, x), NaN rate included

    reason = np.full(n, 3, dtype=np.int8)
    reason[~have] = 1
    reason[have & ~centre_ok] = 2
    reason[found] = 0
    return adjusted, reason


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--modkit', '-m', required=True, help='modkit bedMethyl to correct')
    ap.add_argument('--reference', '-ref', required=True, help='uncompressed reference FASTA, with .fai')
    ap.add_argument('--errortable', '-e', required=True, help='gIVT false-positive rates per 9mer (can be gzipped)')
    ap.add_argument('--mod_threshold', '-mt', nargs='+', help='mod code and threshold, e.g. m,0.7 a,0.8')
    ap.add_argument('--outpath', '-o', required=True, help='output bedMethyl with corrected percent as column 19')
    ap.add_argument('--min_coverage', '-n', type=int, default=1, help='only write rows with Nvalid_cov >= this')
    ap.add_argument('--valid_only', action='store_true', help='drop rows that cannot be corrected (= both --filter_*)')
    ap.add_argument('--filter_kmer', action='store_true', help='drop rows without a 9mer')
    ap.add_argument('--filter_mismatch', action='store_true', help="drop rows whose centre base isn't the mod's base")
    ap.add_argument('--engine', choices=['polars', 'pandas'], default='polars', help='input parsing (default polars)')
    args = ap.parse_args()

    filter_kmer = args.filter_kmer or args.valid_only
    filter_mismatch = args.filter_mismatch or args.valid_only
    thresholds = defaultdict(lambda: 0.7)
    for pair in args.mod_threshold or []:
        mod, value = pair.split(',')
        thresholds[mod] = float(value)

    if args.engine == 'polars':
        try:
            import polars  # noqa: F401
        except ImportError:
            sys.exit("fp_correct.py: polars is not installed here; use --engine pandas")

    t = time.time()
    index, fps, present = load_table(args.errortable, thresholds)
    log(f"error table: {len(index)} mod code(s) in {time.time() - t:.1f}s")

    if os.path.getsize(args.modkit) == 0:
        open(args.outpath, 'w').close()
        log("empty input: wrote an empty output")
        return

    t = time.time()
    lines, sites = read_sites(args.modkit, args.engine)
    log(f"read {len(lines)} sites ({args.engine}) in {time.time() - t:.1f}s")

    t = time.time()
    adjusted, reason = correct(sites, args.reference, index, fps, present)
    log(f"corrected in {time.time() - t:.1f}s")

    drop = ((reason == 1) & filter_kmer) | ((reason == 2) & filter_mismatch) | ((reason == 3) & args.valid_only)
    keep = ~drop & (sites['cov'] >= args.min_coverage)

    t = time.time()
    write_sites(lines, keep, adjusted, args.outpath)
    log(f"wrote {int(keep.sum())} rows in {time.time() - t:.1f}s")

    counts = np.bincount(reason, minlength=4)
    print(f"sites {len(reason)}: corrected {counts[0]}, no 9mer {counts[1]}, base mismatch {counts[2]}, "
          f"9mer/mod not in table {counts[3]}; written {int(keep.sum())}")


if __name__ == '__main__':
    main()
