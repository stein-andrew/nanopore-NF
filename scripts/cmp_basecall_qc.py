#!/usr/bin/env python3
"""cmp_basecall_qc.py: basecall accuracy and poly(A) estimates per run, paired on the same reads.

    python cmp_basecall_qc.py --sample GM12878 --out basecall_qc/GM12878 [--frac 0.1] [--seed 42] \\
        --run A_d212_v6=out_A/bam/GRCh38/GM12878.GRCh38.bam --run B_d130_v52=... --run C_d212_v52=...

Primary alignments only (samtools view -F 2308), and the same read-name-hash sample in every run (samtools view -s SEED.FRAC).
Per read, from the CIGAR (M/=/X, I, D; N and clips excluded) and NM (mismatches + inserted + deleted bases):
  identity = matches / (matches + mismatches + insertions + deletions), as in Ghohabi Esfahani et al. 2026, Fig. 1;
  mismatch, insertion and deletion rates; aligned length; MAPQ; dorado's poly(A) estimate (pt:i; < 0 or absent = not called).
Outputs in --out:
  tables/per_read_<run>.tsv.gz   per-read values
  tables/per_run.tsv             medians per run
  tables/paired.tsv              same reads, run x vs run y: identity difference, poly(A) agreement
  figures/identity.{pdf,png}     identity per run (density, medians) and paired identity differences
  figures/polya.{pdf,png}        poly(A) estimates, run x vs run y, on the same reads
"""
import argparse
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

SHORT = {'A_d212_v6': 'v6 / 2.1.2', 'B_d130_v52': 'v5.2 / 1.3.0', 'C_d212_v52': 'v5.2 / 2.1.2'}
COLORS = {'A_d212_v6': '#1b9e77', 'B_d130_v52': '#7570b3', 'C_d212_v52': '#d95f02'}
PAIRS = [('A_d212_v6', 'C_d212_v52'), ('C_d212_v52', 'B_d130_v52'), ('A_d212_v6', 'B_d130_v52')]
CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def per_read(bam, frac, seed):
    sub = f'{seed}' + f'{frac:.8f}'.rstrip('0')[1:]
    p = subprocess.Popen(['samtools', 'view', '-F', '2308', '-s', sub, bam], stdout=subprocess.PIPE, text=True)
    rows = []
    for line in p.stdout:
        f = line.rstrip('\n').split('\t')
        m = i = d = 0
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in 'M=X':
                m += n
            elif op == 'I':
                i += n
            elif op == 'D':
                d += n
        nm = pt = None
        for t in f[11:]:
            if t.startswith('NM:i:'):
                nm = int(t[5:])
            elif t.startswith('pt:i:'):
                pt = int(t[5:])
        if nm is None or m == 0:
            continue
        mis = max(nm - i - d, 0)
        tot = m + i + d
        rows.append((f[0], (m - mis) / tot, mis / tot, i / tot, d / tot, m + i, int(f[4]),
                     pt if (pt is not None and pt > 0) else np.nan))
    if p.wait() != 0:
        sys.exit(f'cmp_basecall_qc.py: samtools view failed on {bam}')
    return pd.DataFrame(rows, columns=['read_id', 'identity', 'mismatch_rate', 'insertion_rate', 'deletion_rate',
                                       'aligned_query_bases', 'mapq', 'polya'])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--run', action='append', required=True, metavar='NAME=BAM')
    ap.add_argument('--out', required=True)
    ap.add_argument('--frac', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=42)
    a = ap.parse_args()
    runs, bams = [], {}
    for item in a.run:
        name, _, path = item.partition('=')
        if not os.path.exists(path):
            sys.exit(f'cmp_basecall_qc.py: --run {item}: no such file')
        runs.append(name); bams[name] = path
    os.makedirs(f'{a.out}/tables', exist_ok=True); os.makedirs(f'{a.out}/figures', exist_ok=True)

    t = {}
    for r in runs:
        t[r] = per_read(bams[r], a.frac, a.seed)
        t[r].to_csv(f'{a.out}/tables/per_read_{r}.tsv.gz', sep='\t', index=False, float_format='%.5g')
        print(f'{r}: {len(t[r]):,} primary alignments sampled', file=sys.stderr)
    per_run = pd.DataFrame([{'sample': a.sample, 'run': r, 'reads': len(d), 'median_identity': d['identity'].median(),
                             'mean_identity': d['identity'].mean(), 'median_mismatch_rate': d['mismatch_rate'].median(),
                             'median_insertion_rate': d['insertion_rate'].median(), 'median_deletion_rate': d['deletion_rate'].median(),
                             'median_aligned_bases': d['aligned_query_bases'].median(), 'mapq60_share': (d['mapq'] >= 60).mean(),
                             'polya_called_share': d['polya'].notna().mean(), 'median_polya': d['polya'].median()} for r, d in t.items()])
    per_run.to_csv(f'{a.out}/tables/per_run.tsv', sep='\t', index=False, float_format='%.5g')

    pairs = [(x, y) for x, y in PAIRS if x in runs and y in runs] or \
            [(runs[i], runs[j]) for i in range(len(runs)) for j in range(i + 1, len(runs))]
    paired, merged = [], {}
    for x, y in pairs:
        m = t[x].merge(t[y], on='read_id', suffixes=('_x', '_y')); merged[(x, y)] = m
        dd = m['identity_x'] - m['identity_y']
        pa = m.dropna(subset=['polya_x', 'polya_y'])
        paired.append({'sample': a.sample, 'pair': f'{x} vs {y}', 'shared_reads': len(m),
                       'median_identity_x': m['identity_x'].median(), 'median_identity_y': m['identity_y'].median(),
                       'median_identity_diff_x_minus_y': dd.median(), 'share_reads_x_higher': (dd > 0).mean(),
                       'share_reads_equal': (dd == 0).mean(), 'polya_both_called': len(pa),
                       'polya_called_x_only': int((m['polya_x'].notna() & m['polya_y'].isna()).sum()),
                       'polya_called_y_only': int((m['polya_y'].notna() & m['polya_x'].isna()).sum()),
                       'polya_pearson_r': pa['polya_x'].corr(pa['polya_y']) if len(pa) > 2 else np.nan,
                       'polya_median_diff_x_minus_y': (pa['polya_x'] - pa['polya_y']).median() if len(pa) else np.nan,
                       'polya_median_abs_diff': (pa['polya_x'] - pa['polya_y']).abs().median() if len(pa) else np.nan})
    pd.DataFrame(paired).to_csv(f'{a.out}/tables/paired.tsv', sep='\t', index=False, float_format='%.5g')

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.8), layout='constrained')
    bins = np.linspace(0.80, 1.0, 81)
    for r in runs:
        v = t[r]['identity']
        a1.hist(v.clip(0.80, 1.0), bins=bins, density=True, histtype='step', lw=1.5, color=COLORS.get(r), label=f'{SHORT.get(r, r)} (median {v.median():.4f})')
        a1.axvline(v.median(), color=COLORS.get(r), ls='--', lw=0.8)
    a1.set_xlabel('read identity'); a1.set_ylabel('density'); a1.legend(fontsize=7, frameon=False, loc='upper left')
    a1.set_title(f'{a.sample}: alignment identity, primary alignments', fontsize=10)
    db = np.linspace(-0.03, 0.03, 61)
    for x, y in pairs:
        m = merged[(x, y)]; dd = (m['identity_x'] - m['identity_y']).clip(-0.03, 0.03)
        a2.hist(dd, bins=db, histtype='step', lw=1.5, label=f'{SHORT.get(x, x)} − {SHORT.get(y, y)} (median {dd.median():+.4f})')
    a2.axvline(0, color='0.4', lw=0.8); a2.set_xlabel('identity difference, same read'); a2.set_ylabel('reads')
    a2.legend(fontsize=7, frameon=False); a2.set_title('paired on the same reads', fontsize=10)
    fig.savefig(f'{a.out}/figures/identity.pdf'); fig.savefig(f'{a.out}/figures/identity.png', dpi=200); plt.close(fig)

    fig, axes = plt.subplots(1, len(pairs), figsize=(4.2 * len(pairs), 3.8), squeeze=False, layout='constrained')
    for ax, (x, y) in zip(axes[0], pairs):
        pa = merged[(x, y)].dropna(subset=['polya_x', 'polya_y'])
        if len(pa):
            hi = np.nanpercentile(np.r_[pa['polya_x'], pa['polya_y']], 99.5)
            ax.hist2d(pa['polya_y'].clip(0, hi), pa['polya_x'].clip(0, hi), bins=80, norm=LogNorm(), cmap='viridis')
            ax.plot([0, hi], [0, hi], 'w--', lw=0.8)
            ax.set_title(f'r = {pa["polya_x"].corr(pa["polya_y"]):.3f}, n = {len(pa):,}', fontsize=9)
        ax.set_xlabel(f'poly(A), {SHORT.get(y, y)}'); ax.set_ylabel(f'poly(A), {SHORT.get(x, x)}')
    fig.suptitle(f'{a.sample}: dorado poly(A) estimates on the same reads', fontsize=10)
    fig.savefig(f'{a.out}/figures/polya.pdf'); fig.savefig(f'{a.out}/figures/polya.png', dpi=200); plt.close(fig)
    print(f'wrote {a.out}/tables and figures', file=sys.stderr)


if __name__ == '__main__':
    main()
