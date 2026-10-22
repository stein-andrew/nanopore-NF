#!/usr/bin/env python3
"""One set of per-base call thresholds, written out for both tools that need it:
modkit pileup (calls a modification at a base when its probability >= the threshold)
and gIVT's IVT_fp_correction.py (subtracts the false-positive rate measured at that same
threshold). gIVT's correction is only right if the two agree, hence one source.

    mod_thresholds.py --sample-probs S.sample_probs.txt --percentile 0.1 --prefix S
    mod_thresholds.py --sample-probs S.sample_probs.txt --percentile 0.1 --prefix S \
        --fixed 'A:0.9 C:0.85 G:0.95 T:0.9'

--sample-probs: modkit sample-probs' threshold table (stdout or its --out-dir file). Columns
    are found by header name (base, threshold, optional percentile), so tab, space or
    box-drawn layouts all parse. These are the thresholds modkit pileup would estimate
    itself; they become the ones used unless --fixed is given.
--fixed: one threshold for each of A, C, G, T, used for every sample instead. The
    sample-probs estimate is still recorded next to it (or NA if the table won't parse).

Writes:
  <prefix>.thresholds.tsv   base, estimated, used, source
  <prefix>.modkit_flags.txt --filter-threshold A:x --filter-threshold C:x ... (one line)
  <prefix>.givt_flags.txt   every RNA mod code with its base's threshold: a,x 17596,x ...

Every used threshold must be in [0.70, 0.99): gIVT's error table starts at 0.70, and in
its mod_to_col() a threshold >= 0.99 finds no larger column and silently falls back to
the 0.70 false-positive rates (below 0.70 it also gets the 0.70 rates). Out of range is
an error here, not a quiet miscorrection there.
"""
import argparse
import re
import sys

BASES = 'ACGT'
# RNA mod codes Dorado/modkit report, by the canonical base they modify (as in gIVT).
CODES = {'A': ['a', '17596', '69426'], 'C': ['m', '19228'], 'G': ['19229'], 'T': ['17802', '19227']}
LOW, HIGH = 0.70, 0.99


def tokens(line):
    """Split a table line on tabs, spaces and box-drawing/pipe characters."""
    return [t for t in re.split(r'[\s|│┃]+', line.strip()) if t and not re.fullmatch(r'[-─━+=┼]+', t)]


def parse_sample_probs(path, percentile):
    """base -> threshold from a sample-probs table, or raise ValueError saying why not."""
    lines = open(path).read().splitlines()
    head = None
    for i, line in enumerate(lines):
        low = [t.lower() for t in tokens(line)]
        if ('base' in low or 'primary_base' in low) and 'threshold' in low:
            head = i, low
            break
    if head is None:
        raise ValueError(f"{path}: no header with base and threshold columns")
    i0, cols = head
    b = cols.index('base') if 'base' in cols else cols.index('primary_base')
    t = cols.index('threshold')
    p = cols.index('percentile') if 'percentile' in cols else None
    want = {round(percentile, 6), round(percentile * 100, 6)}      # 0.1 or 10, either style

    found = {}
    for line in lines[i0 + 1:]:
        row = tokens(line)
        if len(row) <= max(b, t) or row[b] not in BASES:
            continue
        if p is not None:
            try:
                if round(float(row[p]), 6) not in want:
                    continue
            except (ValueError, IndexError):
                continue
        value = float(row[t])
        if row[b] in found and found[row[b]] != value:
            raise ValueError(f"{path}: several thresholds for base {row[b]} at this percentile "
                             f"(per-mod-code rows?); set --fixed")
        found[row[b]] = value
    missing = [x for x in BASES if x not in found]
    if missing:
        raise ValueError(f"{path}: no threshold for base(s) {','.join(missing)}")
    return found


def parse_fixed(text):
    found = {}
    for item in text.split():
        base, sep, value = item.partition(':')
        if not sep or base not in BASES:
            sys.exit(f"--fixed expects 'A:0.9 C:0.85 G:0.95 T:0.9', got item: {item}")
        found[base] = float(value)
    missing = [x for x in BASES if x not in found]
    if missing:
        sys.exit(f"--fixed: no threshold for base(s) {','.join(missing)}; all of A C G T are needed")
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sample-probs', required=True, help='modkit sample-probs threshold table')
    ap.add_argument('--percentile', type=float, required=True, help='the --percentiles value given to sample-probs')
    ap.add_argument('--fixed', default='', help="'A:x C:x G:x T:x' to use instead of the estimate")
    ap.add_argument('--prefix', required=True)
    args = ap.parse_args()

    try:
        estimated = parse_sample_probs(args.sample_probs, args.percentile)
    except (OSError, ValueError) as e:
        if not args.fixed:
            sys.exit(f"mod_thresholds.py: {e}")
        print(f"mod_thresholds.py: warning: {e}; recording the estimate as NA", file=sys.stderr)
        estimated = {}

    used = parse_fixed(args.fixed) if args.fixed else estimated
    source = 'fixed' if args.fixed else 'sample-probs'
    bad = {b: v for b, v in used.items() if not LOW <= v < HIGH}
    if bad:
        sys.exit(f"mod_thresholds.py: threshold(s) {bad} outside [{LOW}, {HIGH}): gIVT's error table cannot "
                 f"correct calls made there (see this script's docstring). Set mod_filter_thresholds explicitly.")

    with open(f"{args.prefix}.thresholds.tsv", 'w') as out:
        out.write("base\testimated\tused\tsource\n")
        for b in BASES:
            est = f"{estimated[b]:.4f}" if b in estimated else 'NA'
            out.write(f"{b}\t{est}\t{used[b]:.4f}\t{source}\n")
    with open(f"{args.prefix}.modkit_flags.txt", 'w') as out:
        out.write(' '.join(f"--filter-threshold {b}:{used[b]}" for b in BASES) + '\n')
    with open(f"{args.prefix}.givt_flags.txt", 'w') as out:
        out.write(' '.join(f"{code},{used[b]}" for b in BASES for code in CODES[b]) + '\n')
    print(f"thresholds ({source}): " + ' '.join(f"{b}={used[b]}" for b in BASES))


if __name__ == '__main__':
    main()
