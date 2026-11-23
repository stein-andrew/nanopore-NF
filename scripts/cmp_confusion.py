#!/usr/bin/env python3
"""cmp_confusion.py: read-level modification-call confusion between basecalling setups on the SAME reads.

    python cmp_confusion.py --sample GM12878 \\
        --run A_d212_v6=out_A/bam/GRCh38/GM12878.GRCh38.bam \\
        --run B_d130_v52=out_B/bam/GRCh38/GM12878.GRCh38.bam \\
        --run C_d212_v52=out_C/bam/GRCh38/GM12878.GRCh38.bam \\
        --out confusion/GM12878 [--frac 0.01] [--seed 42] [--threads 8]

For each run, the same random fraction of reads is taken (samtools view -s SEED.FRAC selects by a hash of the
read name, so every run keeps the same reads), and `modkit extract calls --no-filtering --mapped-only` writes
one row per aligned base: the model's most likely class (canonical '-' or a modification code) and its
probability. No threshold is applied, so the comparison does not depend on any filtering rule. Bases dorado
skipped (inferred canonical, '.' MM mode) are kept as canonical calls and counted separately.

Each pair of runs is joined on (read, reference position), so a cell counts the same base of the same molecule
called by both models. Outputs in --out:
  tables/confusion_<x>_vs_<y>.tsv      base, call_x, call_y, n, frac_of_x_call (row-normalised)
  tables/agreement_by_confidence.tsv   per pair and base: share of identical calls by run x's call probability
  tables/summary.tsv                   per pair and base: positions compared, share identical, top changes
  figures/confusion.pdf                one page per pair: a row-normalised matrix per base (A, C, G, U)
Extraction tables are cached in <out>/extract/ and reused when rerun.
"""
import argparse
import os
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

CODE_NAMES = {'-': 'canonical', 'a': 'm6A', '17596': 'Inosine', '69426': '2OmeA', 'm': 'm5C', '19228': '2OmeC',
              '17802': 'Psi', '19227': '2OmeU', '19229': '2OmeG'}
BASE_CODES = {'A': ['-', 'a', '17596', '69426'], 'C': ['-', 'm', '19228'], 'G': ['-', '19229'], 'T': ['-', '17802', '19227']}
SHORT = {'A_d212_v6': 'v6 / 2.1.2', 'B_d130_v52': 'v5.2 / 1.3.0', 'C_d212_v52': 'v5.2 / 2.1.2'}
PREFERRED_PAIRS = [('A_d212_v6', 'C_d212_v52', 'models (same dorado)'),
                   ('C_d212_v52', 'B_d130_v52', 'dorado version (same models)'),
                   ('A_d212_v6', 'B_d130_v52', 'new vs production')]
KEEP = ['read_id', 'ref_position', 'chrom', 'call_prob', 'call_code', 'canonical_base', 'inferred', 'within_alignment']
CONF_BINS = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001]


def extract(run, bam, a, image, binds):
    """modkit extract calls on a read-name-hash sample of the BAM; cached as a gzipped TSV of the kept columns."""
    out = os.path.join(a.out, 'extract', f'{run}.calls.tsv.gz')
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print(f'{run}: using cached {out}', file=sys.stderr)
        return out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    frac = f'{a.seed}' + f'{a.frac:.8f}'.rstrip('0')[1:]          # samtools -s SEED.FRACTION, e.g. 42.01
    tmp = out + '.partial'
    cmd = (f"samtools view -u -@ 2 -s {frac} '{bam}' "
           f"| modkit extract calls - - --mapped-only --no-filtering --threads {a.threads} "
           f"--log-filepath '{out}.log' --suppress-progress")
    wrapped = ['apptainer', 'exec'] + sum([['-B', b] for b in binds], []) + [image, 'bash', '-c', cmd] if image else ['bash', '-c', cmd]
    print(f'{run}: extracting ({a.frac:g} of reads, seed {a.seed})', file=sys.stderr)
    with subprocess.Popen(wrapped, stdout=subprocess.PIPE) as p:
        header = p.stdout.readline().decode().rstrip('\n').split('\t')
        missing = [c for c in KEEP if c not in header]
        if missing:
            p.kill()
            sys.exit(f'cmp_confusion.py: modkit extract calls output lacks columns {missing}; header: {header}')
        idx = [header.index(c) for c in KEEP]
        import gzip
        with gzip.open(tmp, 'wt', compresslevel=3) as fh:
            fh.write('\t'.join(KEEP) + '\n')
            for line in p.stdout:
                f = line.decode().rstrip('\n').split('\t')
                fh.write('\t'.join(f[i] for i in idx) + '\n')
        if p.wait() != 0:
            sys.exit(f'cmp_confusion.py: extraction failed for {run}; see {out}.log')
    os.replace(tmp, out)
    return out


def load(path, read_codes):
    df = pd.read_csv(path, sep='\t', dtype={'read_id': str, 'ref_position': np.int64, 'chrom': str, 'call_prob': np.float64,
                                            'call_code': str, 'canonical_base': str, 'inferred': str, 'within_alignment': str})
    df = df[(df['ref_position'] >= 0) & (df['within_alignment'] != 'false')]
    # one integer key per (read, chromosome, position); read codes shared across runs
    new = pd.Index(df['read_id'].unique()).difference(read_codes.index)
    if len(new):
        read_codes = pd.concat([read_codes, pd.Series(np.arange(len(read_codes), len(read_codes) + len(new)), index=new)])
    df['key'] = read_codes.reindex(df['read_id']).to_numpy().astype(np.int64) * (1 << 34) + df['ref_position'].to_numpy()
    df = df.drop_duplicates('key')
    return df[['key', 'chrom', 'canonical_base', 'call_code', 'call_prob', 'inferred']], read_codes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--run', action='append', required=True, metavar='NAME=BAM')
    ap.add_argument('--out', required=True)
    ap.add_argument('--frac', type=float, default=0.01, help='fraction of reads to sample (same reads in every run)')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--image', default='/scratch/stein.an/ont-pipeline/containers/dorado-2.1.2-dna.sif',
                    help="image with modkit 0.6.x and samtools; '' to use tools on PATH")
    ap.add_argument('--bind', action='append', default=['/scratch', '/projects'])
    a = ap.parse_args()
    if not 0 < a.frac < 1:
        sys.exit('cmp_confusion.py: --frac must be between 0 and 1')
    if a.image and not shutil.which('apptainer'):
        sys.exit('cmp_confusion.py: apptainer not found; pass --image "" to use modkit/samtools on PATH')

    runs, bams = [], {}
    for item in a.run:
        name, _, path = item.partition('=')
        if not os.path.exists(path):
            sys.exit(f'cmp_confusion.py: --run {item}: no such file')
        runs.append(name); bams[name] = path
    short = {r: SHORT.get(r, r) for r in runs}
    pairs = [(x, y, w) for x, y, w in PREFERRED_PAIRS if x in runs and y in runs] or \
            [(runs[i], runs[j], '') for i in range(len(runs)) for j in range(i + 1, len(runs))]
    os.makedirs(f'{a.out}/tables', exist_ok=True); os.makedirs(f'{a.out}/figures', exist_ok=True)

    tables, read_codes = {}, pd.Series(dtype=np.int64)
    for r in runs:
        tables[r], read_codes = load(extract(r, bams[r], a, a.image, a.bind), read_codes)
        print(f'{r}: {len(tables[r]):,} aligned base calls', file=sys.stderr)

    summ, conf_rows = [], []
    with PdfPages(f'{a.out}/figures/confusion.pdf') as pdf:
        for x, y, why in pairs:
            m = tables[x].merge(tables[y][['key', 'call_code', 'call_prob', 'inferred', 'canonical_base']],
                                on='key', suffixes=('_x', '_y'))
            m = m[m['canonical_base_x'] == m['canonical_base_y']]          # same read base in both basecalls
            rows = []
            fig, axes = plt.subplots(1, 4, figsize=(19, 5.2), layout='constrained')
            for ax, base in zip(axes, ['A', 'C', 'G', 'T']):
                d = m[m['canonical_base_x'] == base]
                codes = [c for c in BASE_CODES[base] if c in set(d['call_code_x']) | set(d['call_code_y'])] or ['-']
                tab = pd.crosstab(d['call_code_x'], d['call_code_y']).reindex(index=codes, columns=codes, fill_value=0)
                frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0)
                for cx in codes:
                    for cy in codes:
                        rows.append({'base': base, 'call_x': CODE_NAMES.get(cx, cx), 'call_y': CODE_NAMES.get(cy, cy),
                                     'n': int(tab.loc[cx, cy]), 'frac_of_x_call': frac.loc[cx, cy]})
                same = (d['call_code_x'] == d['call_code_y']).mean() if len(d) else np.nan
                changes = (d[d['call_code_x'] != d['call_code_y']].groupby(['call_code_x', 'call_code_y']).size()
                           .sort_values(ascending=False).head(3))
                summ.append({'pair': f'{x} vs {y}', 'comparison': why, 'base': base, 'positions': len(d),
                             'identical_share': same,
                             'inferred_x': (d['inferred_x'] == 'true').mean() if len(d) else np.nan,
                             'inferred_y': (d['inferred_y'] == 'true').mean() if len(d) else np.nan,
                             'top_changes': '; '.join(f'{CODE_NAMES.get(i, i)}->{CODE_NAMES.get(j, j)}: {n:,}'
                                                      for (i, j), n in changes.items())})
                if len(d):
                    b = pd.cut(d['call_prob_x'], CONF_BINS, right=False)
                    g = (d['call_code_x'] == d['call_code_y']).groupby(b, observed=False).agg(['mean', 'size'])
                    for interval, row in g.iterrows():
                        conf_rows.append({'pair': f'{x} vs {y}', 'base': base, 'call_prob_x_bin': str(interval),
                                          'positions': int(row['size']), 'identical_share': row['mean']})
                im = ax.imshow(frac.to_numpy(dtype=float), cmap='viridis', vmin=0, vmax=1)
                labels = [CODE_NAMES.get(c, c) for c in codes]
                ax.set_xticks(range(len(codes))); ax.set_xticklabels(labels, rotation=45, ha='right')
                ax.set_yticks(range(len(codes))); ax.set_yticklabels(labels)
                for i, cx in enumerate(codes):
                    for j, cy in enumerate(codes):
                        v = frac.loc[cx, cy]
                        ax.text(j, i, f"{v:.1%}\n{int(tab.loc[cx, cy]):,}" if np.isfinite(v) else '', ha='center', va='center',
                                fontsize=7, color='white' if (np.isfinite(v) and v < 0.6) else 'black')
                ax.set_xlabel(f'{short[y]} call'); ax.set_ylabel(f'{short[x]} call')
                ax.set_title(f"{'U' if base == 'T' else base}: {len(d):,} positions, {same:.2%} identical", fontsize=10)
            fig.colorbar(im, ax=axes.tolist(), shrink=0.7, label=f'share of {short[x]} calls')
            fig.suptitle(f'{a.sample}: {why or short[x] + " vs " + short[y]} (same reads, same positions, no threshold)', fontsize=12)
            pdf.savefig(fig); plt.close(fig)
            pd.DataFrame(rows).to_csv(f'{a.out}/tables/confusion_{x}_vs_{y}.tsv', sep='\t', index=False, float_format='%.5g')
    pd.DataFrame(summ).to_csv(f'{a.out}/tables/summary.tsv', sep='\t', index=False, float_format='%.5g')
    pd.DataFrame(conf_rows).to_csv(f'{a.out}/tables/agreement_by_confidence.tsv', sep='\t', index=False, float_format='%.5g')
    print(f'wrote {a.out}/tables and {a.out}/figures/confusion.pdf', file=sys.stderr)


if __name__ == '__main__':
    main()
