#!/usr/bin/env python3
"""pick_subset.py: choose a subset of each sample's chunks for a quick model-comparison run.

    python pick_subset.py SPLIT_DIR [--reads 1500000] [--min-chunks 2] [--min-per-gb 1.4]

Reads SPLIT_DIR/samplesheet.csv (one row per chunk, written by make_split.sh), counts the reads in each
chunk's pod5 files (file metadata only, so it is quick), and for each sample picks chunks spread evenly
across the run (first to last) until the sample has about --reads reads. Spreading matters: pod5 files in
name order follow the run in time, and pores, read lengths and quality change over a run.

Writes, in SPLIT_DIR:
  samplesheet.subset.csv   the chosen rows, same format (GM12878 first if present), for basecall.nf
  subset_files.txt         the chosen pod5 files (real paths), for copying elsewhere:
                           rsync -av --files-from=subset_files.txt / host:/dest/
and prints per sample: chunks chosen, reads, GB, and estimated H200 hours per run (--min-per-gb is
minutes of H200 time per GB of pod5 with all four modification models: 1.4 measured for dorado 2.1.2
with the v6 set; dorado 1.3.0 may be slower, and an A100 roughly half as fast).
"""
import argparse
import csv
import math
import os
import sys

import pod5


def chunk_files(path):
    files = []
    for d, _, fs in os.walk(path, followlinks=True):
        files += [os.path.realpath(os.path.join(d, f)) for f in fs if f.endswith('.pod5')]
    return sorted(files)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('split_dir')
    ap.add_argument('--reads', type=float, default=1.5e6, help='target reads per sample (default 1.5 million)')
    ap.add_argument('--min-chunks', type=int, default=2, help='at least this many chunks per sample (default 2)')
    ap.add_argument('--min-per-gb', type=float, default=1.4, help='H200 minutes per GB of pod5 (default 1.4)')
    a = ap.parse_args()

    sheet = os.path.join(a.split_dir, 'samplesheet.csv')
    with open(sheet, newline='') as fh:
        rd = csv.DictReader(fh)
        header, rows = rd.fieldnames, list(rd)
    for k in ('sample_id', 'run_id', 'pod5'):
        if k not in header:
            sys.exit(f'pick_subset.py: {sheet} has no {k} column')

    samples = {}
    for r in rows:
        samples.setdefault(r['sample_id'], []).append(r)

    chosen, files_out, report = [], [], []
    for s, rs in samples.items():
        rs = sorted(rs, key=lambda r: r['run_id'])
        info = []
        for r in rs:
            fs = chunk_files(r['pod5'])
            n = 0
            for f in fs:
                with pod5.Reader(f) as rdr:
                    n += rdr.num_reads
            info.append((r, fs, n, sum(os.path.getsize(f) for f in fs)))
            print(f'  counted {r["run_id"]}: {n:,} reads', file=sys.stderr)
        total = sum(x[2] for x in info)
        mean = total / len(info) if info else 0
        def spaced(m):                                              # m chunks evenly spaced, first to last
            if m == 1:
                return [len(info) // 2]
            return sorted({round(i * (len(info) - 1) / (m - 1)) for i in range(m)})
        m = len(info) if mean == 0 else min(len(info), max(a.min_chunks, math.ceil(a.reads / mean)))
        idx = spaced(m)
        while sum(info[i][2] for i in idx) < a.reads and m < len(info):   # chunk sizes vary: top up
            m += 1
            idx = spaced(m)
        pick = [info[i] for i in idx]
        chosen += [x[0] for x in pick]
        files_out += [f for x in pick for f in x[1]]
        reads, gb = sum(x[2] for x in pick), sum(x[3] for x in pick) / 1e9
        report.append((s, len(pick), len(info), reads, total, gb, gb * a.min_per_gb / 60))

    chosen.sort(key=lambda r: (r['sample_id'] != 'GM12878', r['sample_id'], r['run_id']))
    out = os.path.join(a.split_dir, 'samplesheet.subset.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=header, lineterminator='\n')
        w.writeheader()
        w.writerows(chosen)
    with open(os.path.join(a.split_dir, 'subset_files.txt'), 'w') as fh:
        fh.write('\n'.join(files_out) + '\n')

    print(f'{"sample":28s} {"chunks":>9s} {"reads":>12s} {"of all":>7s} {"GB":>7s} {"H200 h/run":>11s}')
    for s, k, n, reads, total, gb, h in sorted(report, key=lambda x: (x[0] != 'GM12878', x[0])):
        print(f'{s:28s} {k:>4d} of {n:<3d} {reads:>12,} {100 * reads / max(total, 1):>6.0f}% {gb:>7.1f} {h:>11.1f}')
    tg, th = sum(x[5] for x in report), sum(x[6] for x in report)
    print(f'{"total":28s} {len(chosen):>4d} chunks {"":>13s} {tg:>7.1f} {th:>11.1f}')
    print(f'wrote {out} and subset_files.txt ({len(files_out)} pod5 files)')


if __name__ == '__main__':
    main()
