#!/usr/bin/env python3
"""cmp_truthsets.py: compare every truth set with every other, and with each nanopore run, per modification.

    python cmp_truthsets.py --sample 08_01_25_RNA004_HEK293 --ortho-dir /projects/Genometechlab/andrew/orthogonal \\
        --fasta GRCh38.fa --mask-bed clusters.bed --out truthsets/08_01_nofilter --label threshold-free \\
        --run A_d212_v6=A.bed --run C_d212_v52=C.bed --run B_d130_v52=B.bed [--tolerance 1]

Truth sets are loaded and checked against the reference exactly as cmp_ortho.py does (same loaders, same consensus
sets). Each nanopore run becomes one more "set": its sites with >= --min-cov valid reads and >= --min-pct occupancy.

Two universes:
  listed      every site as published (nanopore: covered in that run). This is the paper-style comparison and is
              dominated by what each method could detect at all.
  detectable  only positions with >= --min-cov valid reads in every nanopore run, for truth sets and runs alike, so
              every method is judged on the same positions. This is the fair comparison.
Matching: same chromosome and strand, positions within +-tolerance (0 = exact; 1 helps U-runs and CLIP offsets).
A multi-position site (PRAISE ranges) matches if any of its candidate positions does.
Per pair: sites in each, matched sites each way, Jaccard, overlap coefficient (matched / smaller set),
share of A found in B and of B found in A, and Pearson r of levels on exact matches where both report a level.
Outputs: tables/pairs.tsv, tables/sets.tsv, figures/matrix_<mod>.{pdf,png} (share of row set found in column set).
Needs cmp_ortho.py in the same directory.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmp_ortho as co  # noqa: E402

RUNLAB = {'A_d212_v6': 'nanopore v6', 'C_d212_v52': 'nanopore v5.2 (2.1.2)', 'B_d130_v52': 'nanopore v5.2 (1.3.0)'}


def set_table(df):
    """site_id, key, level for one set (one row per candidate position)."""
    return df[['site_id', 'key', 'level']].reset_index(drop=True)


def matched_sites(a, b_keys, tol):
    """site_ids of set a with any candidate within +-tol positions of a key in b (keys step 2 per position)."""
    hit = np.zeros(len(a), dtype=bool)
    for d in range(-tol, tol + 1):
        hit |= np.isin(a['key'].to_numpy() + 2 * d, b_keys)
    return pd.unique(a['site_id'].to_numpy()[hit])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--ortho-dir', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--run', action='append', default=[], metavar='NAME=BED')
    ap.add_argument('--label', default='')
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-cov', type=int, default=20)
    ap.add_argument('--min-pct', type=float, default=20)
    ap.add_argument('--tolerance', type=int, default=0)
    ap.add_argument('--glori-sheet', default='10ng_mRNA_input')
    a = ap.parse_args()
    os.makedirs(f'{a.out}/tables', exist_ok=True); os.makedirs(f'{a.out}/figures', exist_ok=True)

    listed, notes = co.load_sets(a.ortho_dir, a.glori_sheet)
    truth, qc = co.check_against_reference(listed, co.Fasta(a.fasta))
    truth = co.add_consensus(truth)
    mask = pd.read_csv(a.mask_bed, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                       comment='#') if a.mask_bed else None
    runs, tables = [], {}
    for item in a.run:
        name, _, path = item.partition('=')
        runs.append(name)
        tables[name] = co.read_bed(path, sorted({co.CODE[m] for m in truth['mod'].unique()}), a.min_cov, a.contig_regex, mask)
    chroms = sorted(set(truth['chrom']) | set().union(*[set(t['chrom']) for t in tables.values()]) if tables else set(truth['chrom']))
    ids = {c: i for i, c in enumerate(chroms)}
    truth['key'] = co.keys(truth['chrom'], truth['pos'].to_numpy(), truth['strand'], ids)

    pair_rows, set_rows = [], []
    for mod in [m for m in co.MODS if m in set(truth['mod'])]:
        code = co.CODE[mod]
        sets = {}
        for name, d in truth[truth['mod'] == mod].groupby('set', sort=False):
            sets[name] = set_table(d)
        run_keys = {}
        for r in runs:
            t = tables[r][tables[r]['code'] == code]
            k = co.keys(t['chrom'], t['pos'].to_numpy(), t['strand'], ids)
            run_keys[r] = (k, t['pct'].to_numpy())
        universe = None
        if runs:
            universe = run_keys[runs[0]][0]
            for r in runs[1:]:
                universe = np.intersect1d(universe, run_keys[r][0])
        for uni in (['listed', 'detectable'] if runs else ['listed']):
            S = {}
            for name, d in sets.items():
                if uni == 'detectable':
                    keep = d.loc[np.isin(d['key'], universe), 'site_id'].unique()
                    d = d[d['site_id'].isin(keep)]
                S[name] = d
            for r in runs:
                k, pct = run_keys[r]
                sel = pct >= a.min_pct
                if uni == 'detectable':
                    sel &= np.isin(k, universe)
                S[RUNLAB.get(r, r)] = pd.DataFrame({'site_id': [f'{r}:{x}' for x in k[sel]], 'key': k[sel], 'level': pct[sel]})
            names = list(S)
            for n in names:
                set_rows.append({'mod': mod, 'universe': uni, 'set': n, 'sites': S[n]['site_id'].nunique(),
                                 'kind': 'nanopore' if n.startswith('nanopore') else 'orthogonal'})
            for i, x in enumerate(names):
                for y in names[i + 1:]:
                    A, B = S[x], S[y]
                    na, nb = A['site_id'].nunique(), B['site_id'].nunique()
                    if na == 0 or nb == 0:
                        continue
                    ma = len(matched_sites(A, B['key'].to_numpy(), a.tolerance))
                    mb = len(matched_sites(B, A['key'].to_numpy(), a.tolerance))
                    both = min(ma, mb)
                    j = A[['key', 'level']].dropna().drop_duplicates('key').merge(
                        B[['key', 'level']].dropna().drop_duplicates('key'), on='key')
                    r_lv = j['level_x'].corr(j['level_y']) if len(j) >= 10 and j['level_x'].std() > 0 and j['level_y'].std() > 0 else np.nan
                    pair_rows.append({'mod': mod, 'universe': uni, 'set_a': x, 'set_b': y, 'sites_a': na, 'sites_b': nb,
                                      'a_found_in_b': ma, 'b_found_in_a': mb, 'share_a_in_b': ma / na, 'share_b_in_a': mb / nb,
                                      'jaccard': both / (na + nb - both), 'overlap_coef': both / min(na, nb),
                                      'levels_n': len(j), 'levels_r': r_lv})
    pairs, setsdf = pd.DataFrame(pair_rows), pd.DataFrame(set_rows)
    pairs.to_csv(f'{a.out}/tables/pairs.tsv', sep='\t', index=False, float_format='%.4g')
    setsdf.to_csv(f'{a.out}/tables/sets.tsv', sep='\t', index=False)
    qc.to_csv(f'{a.out}/tables/orthogonal_qc.tsv', sep='\t', index=False, float_format='%.4f')

    uni = 'detectable' if runs else 'listed'
    for mod in pairs['mod'].unique() if not pairs.empty else []:
        p = pairs[(pairs['mod'] == mod) & (pairs['universe'] == uni)]
        names = [n for n in dict.fromkeys(list(p['set_a']) + list(p['set_b']))
                 if '∩' not in n and not n.endswith('_2of3')]          # derived consensus sets: tables only
        p = p[p['set_a'].isin(names) & p['set_b'].isin(names)]
        if len(names) < 2:
            continue
        M = np.full((len(names), len(names)), np.nan); N = np.zeros_like(M)
        for _, r in p.iterrows():
            i, j = names.index(r['set_a']), names.index(r['set_b'])
            M[i, j], M[j, i] = r['share_a_in_b'], r['share_b_in_a']
        sz = setsdf[(setsdf['mod'] == mod) & (setsdf['universe'] == uni)].set_index('set')['sites']
        fig, ax = plt.subplots(figsize=(1.1 * len(names) + 3, 0.9 * len(names) + 2), layout='constrained')
        im = ax.imshow(M * 100, cmap='viridis', vmin=0, vmax=100)
        for i in range(len(names)):
            for j in range(len(names)):
                if i != j and np.isfinite(M[i, j]):
                    ax.text(j, i, f'{100 * M[i, j]:.0f}', ha='center', va='center', fontsize=8,
                            color='white' if M[i, j] < 0.6 else 'black')
        ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=40, ha='right', fontsize=8)
        ax.set_yticks(range(len(names))); ax.set_yticklabels([f'{n} (n={sz.get(n, 0):,})' for n in names], fontsize=8)
        ax.set_xlabel('found in this set'); ax.set_ylabel('sites of this set')
        fig.colorbar(im, ax=ax, shrink=0.8, label='% of row set found in column set')
        ax.set_title(f'{mod}: {a.sample}, {uni} positions{", " + a.label if a.label else ""}'
                     f'{", +-" + str(a.tolerance) + " nt" if a.tolerance else ""}', fontsize=10)
        fig.savefig(f'{a.out}/figures/matrix_{mod}.pdf'); fig.savefig(f'{a.out}/figures/matrix_{mod}.png', dpi=200); plt.close(fig)
    print(f'wrote {a.out}/tables and figures', file=sys.stderr)


if __name__ == '__main__':
    main()
