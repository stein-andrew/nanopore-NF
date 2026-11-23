#!/usr/bin/env python3
"""cmp_motif.py: are the m6A sites one model calls and the other doesn't real? Sequence context and GLORI support.

    python cmp_motif.py --sample 08_01_25_RNA004_HEK293 --v6 A.bed --v52 C.bed --fasta GRCh38.fa --mask-bed clusters.bed \\
        --out motif/08_01_nofilter [--ortho-dir /projects/Genometechlab/andrew/orthogonal]

Sites with >= --min-cov valid reads for m6A in both pileups, split at >= --min-pct occupancy into: called by both, v6 only,
v5.2 only; plus a random sample of all covered A positions as background. For each class: share in a DRACH motif
(5-mer [AGT][AG]AC[ACT] on the RNA strand, centred on the A), share in GLORI 1.0 ∩ 2.0 (HEK293 only, with --ortho-dir;
loaded and checked against the reference exactly as cmp_ortho.py does), median occupancy in each model, and the
commonest 5-mers. Real m6A is strongly DRACH-enriched; a v6-only class that looks like the background is likely noise.
Needs cmp_ortho.py in the same directory. Outputs: tables/motif.tsv, figures/motif.{pdf,png}
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmp_ortho as co  # noqa: E402

DRACH = re.compile(r'^[AGT][AG]AC[ACT]$')
COMP = str.maketrans('ACGTN', 'TGCAN')


def kmers(fasta, chrom, pos, strand):
    cols = [fasta.bases(chrom, pos + o) for o in (-2, -1, 0, 1, 2)]
    out = np.array([''.join(b) for b in zip(*cols)], dtype=object)
    minus = strand == '-'
    out[minus] = [k.translate(COMP)[::-1] for k in out[minus]]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--v6', required=True, help='v6 pileup (bedMethyl)')
    ap.add_argument('--v52', required=True, help='v5.2 pileup (bedMethyl), same dorado')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ortho-dir')
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-cov', type=int, default=20)
    ap.add_argument('--min-pct', type=float, default=20)
    ap.add_argument('--background', type=int, default=200_000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--label', default='', help='arm label for the figure title')
    a = ap.parse_args()
    for f in (a.v6, a.v52, a.fasta, a.fasta + '.fai'):
        if not os.path.exists(f):
            sys.exit(f'cmp_motif.py: {f} not found')
    os.makedirs(f'{a.out}/tables', exist_ok=True); os.makedirs(f'{a.out}/figures', exist_ok=True)
    mask = pd.read_csv(a.mask_bed, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                       comment='#') if a.mask_bed else None
    t = {k: co.read_bed(p, ['a'], a.min_cov, a.contig_regex, mask) for k, p in (('v6', a.v6), ('v52', a.v52))}
    chroms = sorted(set(t['v6']['chrom']) | set(t['v52']['chrom']))
    ids = {c: i for i, c in enumerate(chroms)}
    for k in t:
        t[k]['key'] = co.keys(t[k]['chrom'], t[k]['pos'].to_numpy(), t[k]['strand'], ids)
    m = t['v6'][['key', 'chrom', 'pos', 'strand', 'pct']].merge(t['v52'][['key', 'pct']], on='key', suffixes=('_v6', '_v52'))
    hit6, hit52 = m['pct_v6'] >= a.min_pct, m['pct_v52'] >= a.min_pct
    rng = np.random.default_rng(a.seed)
    classes = {'both': m[hit6 & hit52], 'v6 only': m[hit6 & ~hit52], 'v5.2 only': m[~hit6 & hit52],
               'background (all covered A)': m.iloc[rng.choice(len(m), min(a.background, len(m)), replace=False)]}

    glori, support = None, {}
    if a.ortho_dir:
        listed, _ = co.load_sets(a.ortho_dir, '10ng_mRNA_input')
        listed = listed[listed['mod'] == 'm6A']
        truth, _ = co.check_against_reference(listed, co.Fasta(a.fasta))
        truth = co.add_consensus(truth)
        ks = lambda sets_: set(co.keys(*(lambda g: (g['chrom'], g['pos'].to_numpy(), g['strand']))(truth[truth['set'].isin(sets_)]), ids))
        for label, members in [('GLORI 1.0 ∩ 2.0', ['GLORI1.0∩2.0']), ('m6A-SAC-seq', ['m6A-SAC-seq']),
                               ('miCLIP (either antibody)', ['miCLIP-CIMS', 'miCLIP-CITS']), ('DART-seq', ['DART-seq']),
                               ('2+ methods', ['m6A_2+methods'])]:
            if truth['set'].isin(members).any():
                support[label] = ks(members)
        glori = support.get('GLORI 1.0 ∩ 2.0')
    fasta = co.Fasta(a.fasta)
    rows = []
    for name, d in classes.items():
        if d.empty:
            rows.append({'sample': a.sample, 'class': name, 'sites': 0}); continue
        k = kmers(fasta, d['chrom'].to_numpy(), d['pos'].to_numpy(), d['strand'].to_numpy())
        top = pd.Series(k).value_counts().head(6)
        rows.append({'sample': a.sample, 'class': name, 'sites': len(d),
                     'drach_share': float(np.mean([bool(DRACH.match(x)) for x in k])),
                     'central_base_A_share': float(np.mean([len(x) == 5 and x[2] == 'A' for x in k])),
                     'glori_share': float(np.isin(d['key'].to_numpy(), list(glori)).mean()) if glori is not None else np.nan,
                     **{f'share_{lab}': float(np.isin(d['key'].to_numpy(), list(v)).mean()) for lab, v in support.items()},
                     'median_pct_v6': d['pct_v6'].median(), 'median_pct_v52': d['pct_v52'].median(),
                     'top_5mers': ', '.join(f'{x} {n / len(d):.1%}' for x, n in top.items())})
    res = pd.DataFrame(rows)
    res.to_csv(f'{a.out}/tables/motif.tsv', sep='\t', index=False, float_format='%.4g')
    print(res.drop(columns=['top_5mers']).to_string(index=False), file=sys.stderr)

    cols = ['drach_share'] + [f'share_{lab}' for lab in support]
    labels = {'drach_share': 'in DRACH', **{f'share_{lab}': f'in {lab}' for lab in support}}
    palette = ['#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#fdbf6f', '#e31a1c']
    fig, ax = plt.subplots(figsize=(8.5, 3.8), layout='constrained')
    x = np.arange(len(res)); w = 0.8 / len(cols)
    for j, c in enumerate(cols):
        ax.bar(x + (j - (len(cols) - 1) / 2) * w, res[c], w, label=labels[c], color=palette[j % len(palette)])
    ax.set_xticks(x); ax.set_xticklabels([f'{c}\n(n={n:,})' for c, n in zip(res['class'], res['sites'])], fontsize=8)
    ax.set_ylim(0, 1); ax.set_ylabel('share of sites'); ax.legend(fontsize=8, frameon=False)
    ax.set_title(f'{a.sample}: m6A sites >= {a.min_pct:g}% by model {("(" + a.label + ")") if a.label else ""}', fontsize=10)
    fig.savefig(f'{a.out}/figures/motif.pdf'); fig.savefig(f'{a.out}/figures/motif.png', dpi=200); plt.close(fig)


if __name__ == '__main__':
    main()
