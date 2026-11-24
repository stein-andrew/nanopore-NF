#!/usr/bin/env python3
"""cmp_rules.py: compare site-calling rules on no-filter modkit pileups (no new basecalling or pileups needed).

    python cmp_rules.py --mod m6A --label v5.2 \\
        --bed 08_01=<nofilter>/C_d212_v52/08_01_25_RNA004_HEK293.bed --bed 09_09=... --bed 10_06=... \\
        [--negative-bed gIVT=<no-filter pileup of IVT RNA, same model>] \\
        --ortho-dir /projects/Genometechlab/andrew/orthogonal --fasta GRCh38.fa --mask-bed clusters.bed \\
        --out-dir rules/v5.2_m6A

Every rule sees the same positions: >= --min-valid valid reads, masked loci and alt/decoy contigs removed.
Rules (each also requires >= --min-valid valid reads):
  pct20        percent modified >= 20 (the current rule)
  or_k20       pct20, or >= 20 modified reads at any percent
  or_k20_f5    pct20, or >= 20 modified reads and >= 5%
  or_k20_f10   pct20, or >= 20 modified reads and >= 10%
  wilson10     95% Wilson lower bound of the modified fraction >= 10% (depth-aware, ignores background)
  bg5          modified reads significantly above the site's 5-mer background (one-sided binomial test,
               Benjamini-Hochberg q <= --q across the run's sites) and percent modified >= background + --min-effect
Background per 5-mer (RNA strand, centred on the site), estimated from each run itself: the pooled modified
fraction of that context's sites after dropping the top (1 - --trim) of sites by percent modified, so modified
sites barely move it. A 5-mer with fewer than --min-context-sites sites uses its 3-mer (the site and one base either
side), and a sparse 3-mer the modification's pooled value; background_<run>.tsv records which level each used.
The binomial test assumes no extra site-to-site variation within a context; the --min-effect floor is what
keeps that from flooding the high-depth bands.
Reports, per run and rule: sites called; precision and recall against the modification's main truth set
(precision is a lower bound: the maps are incomplete); the same by depth band; agreement between runs at
positions covered in both (share of one run's calls also called in the other, and Jaccard), which needs no
truth set; and, for --negative-bed, calls per million covered positions (every call there is false).
Outputs in --out-dir: summary.tsv, by_depth.tsv, between_runs.tsv, background_<run>.tsv, rules.pdf/png.
Needs cmp_ortho.py in the same directory (truth sets, position keys).
"""
import argparse
import itertools
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmp_ortho as co  # noqa: E402

TRUTH = {'m6A': 'GLORI1.0∩2.0', 'Inosine': 'SLIC-seq_2of3', 'Psi': 'BID∩PRAISE', 'm5C': 'UBS-seq'}
RULES = ['pct20', 'or_k20', 'or_k20_f5', 'or_k20_f10', 'wilson10', 'bg5']
BANDS = [20, 50, 100, 300, 1000, 10000, np.inf]
COMP = np.zeros(256, dtype=np.uint8) + ord('N')
for a_, b_ in zip(b'ACGTN', b'TGCAN'):
    COMP[a_] = b_


class FastaMM:
    """vectorised reference lookups through a memory map of the FASTA and its .fai."""
    def __init__(self, path):
        self.mm = np.memmap(path, dtype=np.uint8, mode='r')
        self.idx = {}
        for line in open(path + '.fai'):
            n, length, off, lb, lw = line.split('\t')[:5]
            self.idx[n] = (int(length), int(off), int(lb), int(lw))

    def kmers(self, chrom, pos1, strand, k=5):
        """k-mer centred on each 1-based position, read on the RNA strand, as an (n, k) uint8 array ('N' off-contig)."""
        h = k // 2
        out = np.full((len(pos1), k), ord('N'), dtype=np.uint8)
        for c in pd.unique(chrom):
            if c not in self.idx:
                continue
            length, off, lb, lw = self.idx[c]
            sel = np.flatnonzero(chrom == c)
            p0 = pos1[sel] - 1
            for j, d in enumerate(range(-h, h + 1)):
                q = p0 + d
                ok = (q >= 0) & (q < length)
                b = np.full(len(sel), ord('N'), dtype=np.uint8)
                b[ok] = self.mm[off + (q[ok] // lb) * lw + q[ok] % lb]
                out[sel, j] = b
        low = (out >= 97) & (out <= 122)
        out[low] -= 32
        minus = np.asarray(strand) == '-'
        out[minus] = COMP[out[minus][:, ::-1]]
        return out


def read_pileup(path, code, min_valid, contig_re, mask):
    parts = []
    for ch in pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 3, 5, 9, 10, 11], chunksize=4_000_000,
                          names=['chrom', 'start', 'code', 'strand', 'nvalid', 'pct', 'nmod'],
                          dtype={'chrom': str, 'start': np.int64, 'code': str, 'strand': str, 'nvalid': np.int64,
                                 'pct': np.float64, 'nmod': np.int64}):
        ch = ch[(ch['code'] == code) & (ch['nvalid'] >= min_valid)]
        if contig_re:
            ch = ch[~ch['chrom'].str.contains(contig_re, regex=True)]
        parts.append(ch)
    d = pd.concat(parts, ignore_index=True)
    d['pos'] = d['start'] + 1
    if mask is not None and len(d):                      # same interval logic as cmp_ortho.read_bed
        drop = np.zeros(len(d), dtype=bool)
        for c, iv in mask.groupby('chrom'):
            sel = np.flatnonzero(d['chrom'].to_numpy() == c)
            if not len(sel):
                continue
            s, e = iv['start'].to_numpy(), iv['end'].to_numpy()
            o = np.argsort(s); s, e = s[o], np.maximum.accumulate(e[o])
            z = d['start'].to_numpy()[sel]
            j = np.searchsorted(s, z, side='right') - 1
            drop[sel] = (j >= 0) & (z < e[np.clip(j, 0, None)])
        d = d[~drop]
    return d.reset_index(drop=True)


def background(kmer, nvalid, nmod, pct, trim, min_sites):
    """per-site background rate: the 5-mer's pooled rate if the context has >= min_sites sites, else the 3-mer's
    (the site and one base either side), else the modification's pooled rate. Each pool drops the top (1 - trim) of
    its sites by percent modified first. Returns the per-site rates and the per-5-mer table (with the level used)."""
    df = pd.DataFrame({'k5': kmer, 'nvalid': nvalid, 'nmod': nmod, 'pct': pct})
    df['k3'] = df['k5'].str[1:4]

    def pooled(key):
        cut = df.groupby(key)['pct'].quantile(trim)
        keep = df['pct'] <= df[key].map(cut)
        g = df[keep].groupby(key).agg(nmod=('nmod', 'sum'), nvalid=('nvalid', 'sum'))
        g['sites'] = df.groupby(key).size()
        g['p0'] = (g['nmod'] / g['nvalid'].clip(lower=1)).clip(lower=1e-6)
        return g

    g5, g3 = pooled('k5'), pooled('k3')
    allk = df['pct'] <= df['pct'].quantile(trim)
    overall = max(df.loc[allk, 'nmod'].sum() / max(df.loc[allk, 'nvalid'].sum(), 1), 1e-6)
    s5, s3 = df['k5'].map(g5['sites']).fillna(0).to_numpy(), df['k3'].map(g3['sites']).fillna(0).to_numpy()
    p5, p3 = df['k5'].map(g5['p0']).to_numpy(), df['k3'].map(g3['p0']).to_numpy()
    p0 = np.where(s5 >= min_sites, p5, np.where(s3 >= min_sites, p3, overall))
    t = g5.reset_index().rename(columns={'k5': 'kmer'})
    t['k3'] = t['kmer'].str[1:4]
    t['p0_3mer'] = t['k3'].map(g3['p0']); t['sites_3mer'] = t['k3'].map(g3['sites'])
    t['level_used'] = np.where(t['sites'] >= min_sites, '5mer', np.where(t['sites_3mer'] >= min_sites, '3mer', 'overall'))
    t['p0_used'] = np.where(t['level_used'] == '5mer', t['p0'], np.where(t['level_used'] == '3mer', t['p0_3mer'], overall))
    return p0, t.assign(overall=overall)


def bh(p):
    n = len(p)
    if not n:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def wilson_lb(k, n, z=1.96):
    p = k / n
    return (p + z * z / (2 * n) - z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n)


def call(d, p0, a):
    pct, nmod, nvalid = d['pct'].to_numpy(), d['nmod'].to_numpy(), d['nvalid'].to_numpy()
    r0 = pct >= 20
    q = bh(binom.sf(nmod - 1, nvalid, p0))
    calls = {'pct20': r0,
             'or_k20': r0 | (nmod >= 20),
             'or_k20_f5': r0 | ((nmod >= 20) & (pct >= 5)),
             'or_k20_f10': r0 | ((nmod >= 20) & (pct >= 10)),
             'wilson10': wilson_lb(nmod, nvalid) >= 0.10,
             'bg5': (q <= a.q) & (pct >= 100 * p0 + a.min_effect)}
    return {r: calls[r] for r in a.rules}


def pairs(arg):
    out = []
    for s in arg or []:
        name, _, path = s.partition('=')
        if not path:
            sys.exit(f'expected NAME=PATH, got {s}')
        out.append((name, path))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--mod', required=True, choices=co.MODS)
    ap.add_argument('--label', default='', help='model label for titles, e.g. v5.2')
    ap.add_argument('--bed', action='append', required=True, metavar='NAME=PATH', help='no-filter pileup of one run (repeat)')
    ap.add_argument('--negative-bed', action='append', metavar='NAME=PATH', help='no-filter pileup of unmodified RNA, same model')
    ap.add_argument('--ortho-dir')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--truth-tsv', help='testing: chrom, pos, strand[, level] instead of the cmp_ortho truth set')
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-valid', type=int, default=20)
    ap.add_argument('--q', type=float, default=0.01)
    ap.add_argument('--min-effect', type=float, default=5, help='percentage points above background, for bg5')
    ap.add_argument('--trim', type=float, default=0.9)
    ap.add_argument('--min-context-sites', type=int, default=30)
    ap.add_argument('--rules', nargs='+', default=RULES, choices=RULES)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    runs, negs = pairs(a.bed), pairs(a.negative_bed)

    truth_keys = None
    st = TRUTH.get(a.mod)
    if a.truth_tsv:
        truth = pd.read_csv(a.truth_tsv, sep='\t')
    elif st and a.ortho_dir:
        listed, _ = co.load_sets(a.ortho_dir, '10ng_mRNA_input')
        listed = listed[listed['mod'] == a.mod]
        truth, _ = co.check_against_reference(listed, co.Fasta(a.fasta))
        truth = co.add_consensus(truth)
        truth = truth[truth['set'] == st]
    else:
        truth = None
    mask = pd.read_csv(a.mask_bed, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                       comment='#') if a.mask_bed else None
    fa = FastaMM(a.fasta)
    ids = {c: i for i, c in enumerate(fa.idx)}
    if truth is not None and len(truth):
        truth_keys = np.unique(co.keys(truth['chrom'], truth['pos'].to_numpy(), truth['strand'], ids))

    summary, by_depth, store = [], [], {}
    labels = [f'{BANDS[i]}-{BANDS[i + 1] - 1}' if np.isfinite(BANDS[i + 1]) else f'{BANDS[i]}+' for i in range(len(BANDS) - 1)]
    for kind, (name, path) in [('run', r) for r in runs] + [('negative', n) for n in negs]:
        d = read_pileup(path, co.CODE[a.mod], a.min_valid, a.contig_regex, mask)
        km = fa.kmers(d['chrom'].to_numpy(), d['pos'].to_numpy(), d['strand'].to_numpy())
        kmer = km.view('S5').ravel().astype(str)
        p0, bgt = background(kmer, d['nvalid'].to_numpy(), d['nmod'].to_numpy(), d['pct'].to_numpy(), a.trim, a.min_context_sites)
        bgt.to_csv(f'{a.out_dir}/background_{name}.tsv', sep='\t', index=False, float_format='%.5g')
        calls = call(d, p0, a)
        keys = co.keys(d['chrom'], d['pos'].to_numpy(), d['strand'], ids)
        in_truth = np.isin(keys, truth_keys) if (truth_keys is not None and kind == 'run') else None
        band = np.digitize(d['nvalid'].to_numpy(), BANDS[1:-1])
        cov_truth = int(in_truth.sum()) if in_truth is not None else 0
        for r, c in calls.items():
            row = {'mod': a.mod, 'label': a.label, 'run': name, 'kind': kind, 'rule': r, 'covered': len(d),
                   'called': int(c.sum()), 'called_per_million_covered': 1e6 * c.sum() / max(len(d), 1),
                   'precision': np.nan, 'recall': np.nan, 'truth_covered': cov_truth}
            if in_truth is not None and cov_truth:
                row['precision'] = in_truth[c].mean() if c.any() else np.nan
                row['recall'] = in_truth[c].sum() / cov_truth
            summary.append(row)
            for b in range(len(labels)):
                sel = band == b
                cb = c & sel
                by_depth.append({'mod': a.mod, 'run': name, 'kind': kind, 'rule': r, 'depth_band': labels[b],
                                 'covered': int(sel.sum()), 'called': int(cb.sum()),
                                 'precision': (in_truth[cb].mean() if (in_truth is not None and cb.any()) else np.nan)})
        if kind == 'run':
            store[name] = (keys, calls)
        print(f'{name}: {len(d):,} positions with >= {a.min_valid} valid reads; median background '
              f'{100 * np.median(p0):.2f}%', file=sys.stderr)

    between = []
    for x, y in itertools.combinations(store, 2):
        kx, cx = store[x]; ky, cy = store[y]
        both, ix, iy = np.intersect1d(kx, ky, return_indices=True)
        for r in a.rules:
            A, Bc = cx[r][ix], cy[r][iy]
            inter = int((A & Bc).sum())
            between.append({'mod': a.mod, 'rule': r, 'run_a': x, 'run_b': y, 'co_covered': len(both),
                            'called_a': int(A.sum()), 'called_b': int(Bc.sum()), 'both': inter,
                            'share_a_in_b': inter / A.sum() if A.any() else np.nan,
                            'share_b_in_a': inter / Bc.sum() if Bc.any() else np.nan,
                            'jaccard': inter / (A | Bc).sum() if (A | Bc).any() else np.nan})
    S, D, Bt = pd.DataFrame(summary), pd.DataFrame(by_depth), pd.DataFrame(between)
    S.to_csv(f'{a.out_dir}/summary.tsv', sep='\t', index=False, float_format='%.4g')
    D.to_csv(f'{a.out_dir}/by_depth.tsv', sep='\t', index=False, float_format='%.4g')
    Bt.to_csv(f'{a.out_dir}/between_runs.tsv', sep='\t', index=False, float_format='%.4g')

    panels = [p for p, ok in (('pr', S['precision'].notna().any()), ('rep', len(Bt) > 0), ('neg', bool(negs))) if ok]
    if panels:
        fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 3.2), layout='constrained', squeeze=False)
        cols = dict(zip(RULES, ['#4d4d4d', '#e41a1c', '#ff7f00', '#984ea3', '#377eb8', '#4daf4a']))
        mk = dict(zip([n for n, _ in runs], 'os^Dv<>p'))
        for ax, p in zip(axes[0], panels):
            if p == 'pr':
                for (run, r), g in S[S['kind'] == 'run'].groupby(['run', 'rule']):
                    ax.scatter(g['recall'], g['precision'], marker=mk[run], s=26, color=cols[r], label=r if run == runs[0][0] else None)
                for run, m_ in mk.items():
                    ax.scatter([], [], marker=m_, color='0.6', s=22, label=f'run {run}')
                ax.set_xlabel('recall'); ax.set_ylabel('precision (lower bound)')
                ax.set_title(f"{a.mod} against {st.replace('∩', ' and ') if st else 'truth'}", fontsize=8)
                ax.legend(frameon=False, fontsize=6)
            elif p == 'rep':
                m = Bt.groupby('rule')[['share_a_in_b', 'share_b_in_a']].mean().mean(axis=1).reindex(a.rules)
                ax.bar(range(len(m)), m.values, color='0.45'); ax.set_xticks(range(len(m))); ax.set_xticklabels(m.index, rotation=40, ha='right', fontsize=7)
                ax.set_ylabel('calls repeated in another run'); ax.set_ylim(0, 1); ax.set_title('Agreement between runs', fontsize=8)
            else:
                m = S[S['kind'] == 'negative'].groupby('rule')['called_per_million_covered'].mean().reindex(a.rules)
                ax.bar(range(len(m)), m.values, color='#b04a4a'); ax.set_xticks(range(len(m))); ax.set_xticklabels(m.index, rotation=40, ha='right', fontsize=7)
                ax.set_ylabel('calls per million covered positions'); ax.set_title('Negative control (all false)', fontsize=8)
        fig.suptitle(f'{a.mod}{", " + a.label if a.label else ""}', fontsize=9)
        fig.savefig(f'{a.out_dir}/rules.pdf'); fig.savefig(f'{a.out_dir}/rules.png', dpi=200); plt.close(fig)
    with pd.option_context('display.width', 220, 'display.max_rows', 200):
        print(S.drop(columns=['mod', 'label']).to_string(index=False, float_format=lambda v: f'{v:.4g}'), file=sys.stderr)
        if len(Bt):
            print(Bt.groupby('rule')[['share_a_in_b', 'share_b_in_a', 'jaccard']].mean().round(3).to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()