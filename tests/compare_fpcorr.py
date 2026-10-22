#!/usr/bin/env python3
"""Check fp_correct.py's output against gIVT's IVT_fp_correction.py, run on the same input.

    compare_fpcorr.py original.bed new.bed

Row for row: columns 1-10 and 12-18 must be identical text; column 11 equal as a number
(gIVT round-trips it through pandas, which rewrites e.g. 50.00 as 50.0); column 19 empty in
both or equal as a number in both. Prints a count per column and the first few mismatches;
exits 1 if there are any.
"""
import sys
from itertools import zip_longest


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} original.bed new.bed")
    orig, new = sys.argv[1], sys.argv[2]
    rows, bad, shown = 0, {}, 0
    with open(orig) as fa, open(new) as fb:
        for n, (a, b) in enumerate(zip_longest(fa, fb), 1):
            if a is None or b is None:
                sys.exit(f"row counts differ: {orig if a is None else new} ends at row {n - 1}")
            rows += 1
            x, y = a.rstrip('\n').split('\t'), b.rstrip('\n').split('\t')
            if len(x) != 19 or len(y) != 19:
                sys.exit(f"row {n}: {len(x)} vs {len(y)} columns, expected 19")
            for c in range(19):
                if c == 10:
                    same = float(x[c]) == float(y[c])
                elif c == 18:
                    same = (x[c] == '' and y[c] == '') or (x[c] != '' and y[c] != '' and float(x[c]) == float(y[c]))
                else:
                    same = x[c] == y[c]
                if not same:
                    bad[c + 1] = bad.get(c + 1, 0) + 1
                    if shown < 5:
                        print(f"row {n} column {c + 1}: original {x[c]!r}  new {y[c]!r}")
                        shown += 1
    if bad:
        print(f"MISMATCH in {rows} rows; per column: {dict(sorted(bad.items()))}")
        sys.exit(1)
    print(f"identical: {rows} rows (column 11 compared as numbers, column 19 as numbers or both empty)")


if __name__ == '__main__':
    main()
