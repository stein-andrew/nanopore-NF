#!/usr/bin/env python3
"""Check fp_correct.py's output against gIVT's IVT_fp_correction.py, run on the same input.

    compare_fpcorr.py original.bed new.bed
    compare_fpcorr.py --run input.bed reference.fa errortable.tsv.gz 'a,0.7 m,0.85 ...' [--image mods.sif]

Row for row: columns 1-10 and 12-18 must be identical text; column 11 equal as a number
(gIVT round-trips it through pandas, which rewrites e.g. 50.00 as 50.0); column 19 empty in
both or equal as a number in both. Prints a count per column and the first few mismatches;
exits 1 if there are any.
--run runs both scripts from bin/ on the input, once with the default flags and once with
--filter_mismatch (the flag MOD_STAGES uses), and compares each pair. fp_correct.py compares
bases upper-cased and gIVT case-sensitively, so gIVT gets an upper-cased copy of the reference
(same line lengths, so the same .fai): on a soft-masked reference this checks the upper-casing
too. --image runs both inside an apptainer image (the pipeline's mods.sif has pysam, pandas and
polars).
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from itertools import zip_longest

BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'bin')
CASES = ['', '--filter_mismatch']


def compare(orig, new):
    """True if the two outputs match as described above; prints the first few differences."""
    rows, bad, shown = 0, {}, 0
    with open(orig) as fa, open(new) as fb:
        for n, (a, b) in enumerate(zip_longest(fa, fb), 1):
            if a is None or b is None:
                print(f"row counts differ: {orig if a is None else new} ends at row {n - 1}")
                return False
            rows += 1
            x, y = a.rstrip('\n').split('\t'), b.rstrip('\n').split('\t')
            if len(x) != 19 or len(y) != 19:
                print(f"row {n}: {len(x)} vs {len(y)} columns, expected 19")
                return False
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
        return False
    print(f"identical: {rows} rows (column 11 compared as numbers, column 19 as numbers or both empty)")
    return True


def run(a):
    tmp = tempfile.mkdtemp(prefix='compare_fpcorr.')
    upper = os.path.join(tmp, 'upper.fa')
    with open(a.reference) as fi, open(upper, 'w') as fo:
        for line in fi:
            fo.write(line if line.startswith('>') else line.upper())
    shutil.copy(a.reference + '.fai', upper + '.fai')
    pre = ['apptainer', 'exec', '--bind', '/projects/Genometechlab,/scratch', a.image] if a.image else []
    env = {**os.environ, 'PYTHONPATH': BIN}
    ok, rows = True, {}
    for flags in CASES:
        tag = flags.lstrip('-') or 'default'
        out = {}
        for script, ref in (('IVT_fp_correction.py', upper), ('fp_correct.py', a.reference)):
            out[script] = os.path.join(tmp, f'{script[:-3]}.{tag}.bed')
            cmd = pre + ['python3', os.path.join(BIN, script), '-m', a.input, '-ref', ref, '-e', a.errortable,
                         '-mt', *a.mod_threshold.split(), '-o', out[script], *flags.split()]
            r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, env=env)
            if r.returncode:
                print(f"[{tag}] {script} failed ({r.returncode}):\n{r.stderr[-2000:]}")
                return False
        print(f"[{tag}] ", end='')
        ok &= compare(out['IVT_fp_correction.py'], out['fp_correct.py'])
        rows[tag] = sum(1 for _ in open(out['fp_correct.py']))
    print(f"rows written: {rows}; --filter_mismatch dropped {rows['default'] - rows['filter_mismatch']}; outputs in {tmp}")
    return ok


def main():
    if len(sys.argv) == 3 and not sys.argv[1].startswith('-'):
        sys.exit(0 if compare(sys.argv[1], sys.argv[2]) else 1)
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--run', nargs=4, required=True, metavar=('INPUT', 'REFERENCE', 'ERRORTABLE', 'MOD_THRESHOLDS'))
    ap.add_argument('--image', help='run both scripts in this apptainer image')
    a = ap.parse_args()
    a.input, a.reference, a.errortable, a.mod_threshold = a.run
    sys.exit(0 if run(a) else 1)


if __name__ == '__main__':
    main()
