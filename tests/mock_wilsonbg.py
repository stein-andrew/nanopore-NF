#!/usr/bin/env python3
"""mock_wilsonbg.py: synthetic m6A no-filter pileups with planted 5-mer backgrounds, to test the wilsonbg rules of
scripts/wilson_sweep.py (--background empirical) before they see real data.

    python tests/mock_wilsonbg.py make --out-dir mock/context [--scenario context|flat] [--seed 1]
    python scripts/wilson_sweep.py --mod m6A --bed 08_01=mock/context/08_01.bed --bed 09_09=... --bed 10_06=... \\
        --fasta mock/context/ref.fa --truth-tsv mock/context/truth.tsv --background empirical --out-dir mock/context/sweep
    python tests/mock_wilsonbg.py summary --mock-dir mock/context --sweep-dir mock/context/sweep

make: random contigs (ref.fa + .fai); every A on either strand is a site. Each A-centred 5-mer (RNA strand) gets a
planted background rate: 'context' draws it log-normally around 2% and lifts 10 contexts to 10-15.5% (real m6A
backgrounds reach 15.6%); 'flat' gives every context 3%. Sites vary around their context's rate (beta, concentration
50) and keep their rate in every run, so false positives repeat between runs as in real data. 2% of sites are truly
modified (DRACH contexts 10 times as likely), at 10-90% stoichiometry. The three runs share sites and site expression;
depth per run is Poisson around expression x a run factor; modified reads ~ Binomial(depth, s + (1 - s) p0_site).
Writes <run>.bed (bedMethyl, code 'a' only) for 08_01, 09_09, 10_06; truth.tsv (chrom, pos, strand, level);
planted_background.tsv (kmer, p0_planted, sites, true_sites, drach, hot); mock_givt07.tsv.gz, a gIVT-format 9-mer table
whose false_positive_0.7_threshold is the planted rate of the 9-mer's centre 5-mer, with support < 200 for about 5% of
9-mers (these must fall back to the empirical background) and decoy rows for code 'm' (rate 0.9, must be ignored).
Run the sweep with --background empirical --background givt07:mock/context/mock_givt07.tsv.gz.
summary: (1) with background 0, wilsonbg<c>_d<d> must equal wilson<c>_t<d>, and every wilsonbg<c>_t10_d<d> must equal
wilson<c>_t10; (1b) cmp_rules.read_pileup(fail_as_canonical=True) on a handmade bedMethyl; (2) on the mock itself, wilsonbg<c>_t10_d<d> with d <= 5 against wilson<c>_t10 (equal wherever every
background is <= 5%, as in the flat mock); (3) estimated against planted background per 5-mer; (4) the gIVT lookup
against an independent one (9-mers cut from the sequence directly); (5) each family against the plain-cutoff curve and
the Wilson 95% curve at matched recall (truth is complete here, so precision is exact); (6) where false calls fall.
"""
import argparse
import itertools
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts'))

RUNS = {'08_01': 1.0, '09_09': 0.9, '10_06': 1.4}
COMP = str.maketrans('ACGT', 'TGCA')


def drach(k):
    return k[0] in 'AGT' and k[1] in 'AG' and k[3] == 'C' and k[4] in 'ACT'


def make(a):
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.out_dir, exist_ok=True)
    seqs = {f'chr{i + 1}': ''.join(rng.choice(list('ACGT'), a.length)) for i in range(a.contigs)}
    with open(f'{a.out_dir}/ref.fa', 'w') as f, open(f'{a.out_dir}/ref.fa.fai', 'w') as fai:
        off = 0
        for name, s in seqs.items():
            head = f'>{name}\n'
            body = ''.join(s[i:i + 60] + '\n' for i in range(0, len(s), 60))
            fai.write(f'{name}\t{len(s)}\t{off + len(head)}\t60\t61\n')
            f.write(head + body)
            off += len(head) + len(body)

    rows = []
    for name, s in seqs.items():
        for p in range(2, len(s) - 2):
            if s[p] == 'A':
                rows.append((name, p, '+', s[p - 2:p + 3]))
            elif s[p] == 'T':
                rows.append((name, p, '-', s[p - 2:p + 3].translate(COMP)[::-1]))
    sites = pd.DataFrame(rows, columns=['chrom', 'start', 'strand', 'kmer'])

    kmers = [x[0] + x[1] + 'A' + x[2] + x[3] for x in itertools.product('ACGT', repeat=4)]
    if a.scenario == 'context':
        p_ctx = np.clip(np.exp(rng.normal(np.log(0.02), 0.6, len(kmers))), 0.002, 0.08)
        hot = rng.choice(len(kmers), 10, replace=False)
        p_ctx[hot] = rng.uniform(0.10, 0.155, len(hot))
    else:
        p_ctx, hot = np.full(len(kmers), 0.03), np.array([], dtype=int)
    planted = pd.DataFrame({'kmer': kmers, 'p0_planted': p_ctx, 'drach': [drach(k) for k in kmers],
                            'hot': np.isin(np.arange(len(kmers)), hot)})
    mu = sites['kmer'].map(planted.set_index('kmer')['p0_planted']).to_numpy()
    p_site = rng.beta(mu * 50, (1 - mu) * 50)

    n = len(sites)
    w = np.where(sites['kmer'].map(planted.set_index('kmer')['drach']).to_numpy(), 10.0, 1.0)
    true = rng.choice(n, int(round(0.02 * n)), replace=False, p=w / w.sum())
    band = rng.choice(3, len(true), p=[0.30, 0.35, 0.35])
    lo, hi = np.array([0.10, 0.20, 0.40])[band], np.array([0.20, 0.40, 0.90])[band]
    stoich = np.zeros(n)
    stoich[true] = rng.uniform(lo, hi)
    expr = np.exp(rng.normal(np.log(60), 1.0, n))

    g = sites.groupby('kmer').size()
    planted['sites'] = planted['kmer'].map(g).fillna(0).astype(int)
    planted['true_sites'] = planted['kmer'].map(sites.iloc[true].groupby('kmer').size()).fillna(0).astype(int)
    planted.to_csv(f'{a.out_dir}/planted_background.tsv', sep='\t', index=False, float_format='%.5g')
    pd.DataFrame({'chrom': sites['chrom'].iloc[true].to_numpy(), 'pos': sites['start'].iloc[true].to_numpy() + 1,
                  'strand': sites['strand'].iloc[true].to_numpy(), 'level': 100 * stoich[true]}
                 ).to_csv(f'{a.out_dir}/truth.tsv', sep='\t', index=False, float_format='%.2f')

    for run, factor in RUNS.items():
        depth = np.minimum(rng.poisson(expr * factor), 20000)
        nmod = rng.binomial(depth, stoich + (1 - stoich) * p_site)
        keep = depth > 0
        d = sites[keep]
        nv, nm = depth[keep], nmod[keep]
        bed = pd.DataFrame({'chrom': d['chrom'], 'start': d['start'], 'end': d['start'] + 1, 'code': 'a', 'score': nv,
                            'strand': d['strand'], 's2': d['start'], 'e2': d['start'] + 1, 'rgb': '255,0,0',
                            'nvalid': nv, 'pct': np.round(100 * nm / nv, 2), 'nmod': nm, 'ncanon': nv - nm,
                            'nother': 0, 'ndel': 0, 'nfail': 0, 'ndiff': 0, 'nnocall': 0})
        bed.sort_values(['chrom', 'start', 'strand']).to_csv(f'{a.out_dir}/{run}.bed', sep='\t', header=False,
                                                             index=False, float_format='%.2f')

    g_rng = np.random.default_rng(a.seed + 1)            # own stream: the sites and runs above stay as they were
    mers = [x[:4] + 'A' + x[4:] for x in map(''.join, itertools.product('ACGT', repeat=8))]
    fp = pd.Series(mers).str[2:7].map(planted.set_index('kmer')['p0_planted']).to_numpy()
    low = g_rng.random(len(mers)) < 0.05
    sup = np.where(low, g_rng.integers(20, 200, len(mers)), g_rng.integers(200, 100_000, len(mers)))
    pd.concat([pd.DataFrame({'mer': mers, 'mod': 'a', 'false_positive_0.7_threshold': fp,
                             'false_positive_0.71_threshold': 0.9 * fp, 'support': sup}),
               pd.DataFrame({'mer': mers, 'mod': 'm', 'false_positive_0.7_threshold': 0.9,
                             'false_positive_0.71_threshold': 0.9, 'support': 10_000})]
              ).to_csv(f'{a.out_dir}/mock_givt07.tsv.gz', sep='\t', index=False, float_format='%.6g')
    print(f'{a.scenario}: {n:,} sites, {len(true):,} true ({np.isin(true, np.flatnonzero(w > 1)).mean():.0%} in DRACH), '
          f'{len(hot)} hot contexts -> {a.out_dir}', file=sys.stderr)


def at_matched(rows, curve):
    """precision minus the curve's precision at the same recall (NaN outside the curve's recall range)."""
    c = curve.sort_values('recall')
    return rows['precision'] - np.interp(rows['recall'], c['recall'], c['precision'], left=np.nan, right=np.nan)


def summary(a):
    import cmp_ortho as co
    import cmp_rules as cr
    import wilson_sweep as ws

    # (1) wilsonbg with a zero background is the plain Wilson rule at t = d
    rng = np.random.default_rng(0)
    nv = rng.integers(20, 2000, 50_000)
    d = pd.DataFrame({'nvalid': nv, 'nmod': rng.binomial(nv, rng.uniform(0, 0.4, len(nv)))})
    d['pct'] = 100 * d['nmod'] / d['nvalid']
    r = ws.rules(d, {'wilsonbg': np.zeros(len(d))})
    pairs_ = [(f'wilsonbg{c}_d{x:g}', f'wilson{c}_t{x:g}') for c in ws.CONF for x in ws.DELTAS if x in ws.TARGETS]
    pairs_ += [(f'wilsonbg{c}_t10_d{x:g}', f'wilson{c}_t10') for c in ws.T10_CONF for x in ws.T10_DELTAS]
    bad = [p for p in pairs_ if not np.array_equal(r[p[0]], r[p[1]])]
    print(f'(1) zero background: {len(pairs_) - len(bad)}/{len(pairs_)} wilsonbg rules identical to the Wilson rule '
          f'they reduce to{"; MISMATCH " + str(bad) if bad else ""}; rules(d) without p0 has no wilsonbg: '
          f'{not any(k.startswith("wilsonbg") for k in ws.rules(d))}')

    # (1b) fail_as_canonical: depth N_valid + N_fail, occupancy N_mod / (N_valid + N_fail)
    rows = [('chr1', 100, 'a', 20, 25.0, 5, 20),      # default: 20 reads, 25%; fail as canonical: 40 reads, 12.5%
            ('chr1', 200, 'a', 10, 30.0, 3, 15),      # default: below 20 reads; fail as canonical: 25 reads, 12%
            ('chr1', 300, 'm', 50, 50.0, 25, 0)]      # another code: never read for 'a'
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
        for c, s, code, nv, pct, nm, nf in rows:
            f.write('\t'.join(map(str, [c, s, s + 1, code, nv, '+', s, s + 1, '255,0,0', nv, f'{pct:.2f}', nm,
                                        nv - nm, 0, 0, nf, 0, 0])) + '\n')
    d0 = cr.read_pileup(f.name, 'a', 20, None, None)
    d1 = cr.read_pileup(f.name, 'a', 20, None, None, fail_as_canonical=True)
    os.unlink(f.name)
    ok = (d0[['start', 'nvalid', 'pct']].values.tolist() == [[100, 20, 25.0]]
          and d1[['start', 'nvalid', 'pct']].values.tolist() == [[100, 40, 12.5], [200, 25, 12.0]])
    print(f'(1b) fail_as_canonical reading: {"as expected" if ok else "WRONG"} (default {d0[["start", "nvalid", "pct"]].values.tolist()}, '
          f'fail as canonical {d1[["start", "nvalid", "pct"]].values.tolist()})')

    # (2) the t10 variants on the mock: equal to wilson<c>_t10 wherever background + d <= 10%
    fa = cr.FastaMM(f'{a.mock_dir}/ref.fa')
    S = pd.read_csv(f'{a.sweep_dir}/sweep.tsv', sep='\t').set_index(['run', 'rule'])
    print('\n(2) wilsonbg<c>_t10_d<d> (d <= 5) against wilson<c>_t10, per run: sites differing (recomputed) and calls in sweep.tsv')
    out, bg = [], {}
    for run in RUNS:
        d = cr.read_pileup(f'{a.mock_dir}/{run}.bed', 'a', 20, None, None)
        km = fa.kmers(d['chrom'].to_numpy(), d['pos'].to_numpy(), d['strand'].to_numpy())
        p0, _ = cr.background(km.view('S5').ravel().astype(str), d['nvalid'].to_numpy(), d['nmod'].to_numpy(),
                              d['pct'].to_numpy(), 0.9, 30)
        bg[run] = (d, km, p0)
        r = ws.rules(d, {'wilsonbg': p0})
        for c in ws.T10_CONF:
            for x in [x for x in ws.T10_DELTAS if x <= 5]:
                rn, ref = f'wilsonbg{c}_t10_d{x:g}', f'wilson{c}_t10'
                out.append({'run': run, 'rule': rn, 'against': ref, 'max_background_pct': 100 * p0.max(),
                            'sites_differing': int((r[rn] != r[ref]).sum()), 'identical': bool(np.array_equal(r[rn], r[ref])),
                            'called_sweep': S.loc[(run, rn), 'called'], 'called_ref_sweep': S.loc[(run, ref), 'called']})
    print(pd.DataFrame(out).round(2).to_string(index=False))

    # (3) background recovery
    pl = pd.read_csv(f'{a.mock_dir}/planted_background.tsv', sep='\t')
    print('\n(3) estimated against planted background per 5-mer, percentage points')
    out = []
    for run in RUNS:
        b = pd.read_csv(f'{a.sweep_dir}/background_{run}.tsv', sep='\t').merge(pl, on='kmer', how='inner')
        e = 100 * (b['p0_used'] - b['p0_planted'])
        out.append({'run': run, 'contexts': len(b), 'level_5mer': (b['level_used'] == '5mer').mean(),
                    'pearson_r': np.corrcoef(b['p0_used'], b['p0_planted'])[0, 1] if b['p0_planted'].nunique() > 1 else np.nan,
                    'median_abs_err': e.abs().median(), 'bias_all': e.mean(),
                    'bias_drach': e[b['drach']].mean(), 'bias_non_drach': e[~b['drach']].mean(),
                    'hot_planted': 100 * b.loc[b['hot'], 'p0_planted'].mean() if b['hot'].any() else np.nan,
                    'hot_estimated': 100 * b.loc[b['hot'], 'p0_used'].mean() if b['hot'].any() else np.nan})
    print(pd.DataFrame(out).round(3).to_string(index=False))

    # (4) gIVT lookup against an independent one: 9-mers cut from the sequence, planted rate where support >= 200
    gpath = f'{a.mock_dir}/mock_givt07.tsv.gz'
    gt = pd.read_csv(gpath, sep='\t')
    ga = gt[gt['mod'] == 'a'].set_index('mer')
    seqs, cur = {}, None
    for line in open(f'{a.mock_dir}/ref.fa'):
        if line.startswith('>'):
            cur = line[1:].strip(); seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    seqs = {c: ''.join(v) for c, v in seqs.items()}
    table = ws.load_givt(gpath, 'a')
    print(f'\n(4) gIVT lookup: load_givt kept {len(table[0]):,} rows (code a; table has {len(gt):,}), '
          f'largest rate {table[0].max():.3f} (decoy rows are 0.9)')
    out = []
    for run, (d, km, p0) in bg.items():
        k9 = []
        for c, p, s in zip(d['chrom'], d['pos'] - 1, d['strand']):
            m = seqs[c][p - 4:p + 5] if p >= 4 else ''
            k9.append(m if s == '+' else m.translate(COMP)[::-1])
        k9 = pd.Series(k9)
        exp_used = (k9.map(ga['support']).fillna(0) >= 200).to_numpy()
        expected = np.where(exp_used, k9.map(ga['false_positive_0.7_threshold']).to_numpy(dtype=float), p0)
        planted5 = k9.str[2:7].map(pl.set_index('kmer')['p0_planted']).to_numpy(dtype=float)
        got, used = ws.givt_p0(table, fa.kmers(d['chrom'].to_numpy(), d['pos'].to_numpy(), d['strand'].to_numpy(), k=9)
                               .view('S9').ravel().astype(str), p0)
        out.append({'run': run, 'sites': len(d), 'table_used': used.mean(), 'fallback': (~used).mean(),
                    'used_mask_identical': bool(np.array_equal(used, exp_used)),
                    'rates_identical': bool(np.array_equal(got, expected)),
                    'max_diff_from_planted_5mer': np.abs(got[used] - planted5[used]).max()})   # table has 6 significant digits
    print(pd.DataFrame(out).round(4).to_string(index=False))
    bs = f'{a.sweep_dir}/background_summary.tsv'
    if os.path.exists(bs):
        print('\nbackground_summary.tsv (percent):')
        b = pd.read_csv(bs, sep='\t')
        for col in ['median_drach', 'median_other', 'median_all', 'max']:
            b[col] = 100 * b[col]
        print(b.round(3).to_string(index=False))

    # (5) matched recall: each family against the plain-cutoff curve and the Wilson 95% curve
    S = pd.read_csv(f'{a.sweep_dir}/sweep.tsv', sep='\t')
    rows = []
    for run, g in S.groupby('run'):
        g = g.copy()
        g['gain_vs_pct'] = at_matched(g, g[g['family'] == 'pct'])
        g['gain_vs_wilson95'] = at_matched(g, g[g['family'] == 'wilson95'])
        rows.append(g)
    M = pd.concat(rows)
    keep = ['pct_15', 'pct_20', 'pct_25', 'wilson95_t7.5', 'wilson95_t10', 'wilson95_t12.5', 'wilsonbg95_d2.5',
            'wilsonbg95_d5', 'wilsonbg95_d7.5', 'wilsonbg95_d10', 'wilsonbg90_d5', 'wilsonbg99_d5', 'wilsonbg95_t10_d2.5',
            'wilsonbg95_t10_d5', 'wilsonbg95_t10_d10', 'wilsonbg_givt07_95_d5', 'wilsonbg_givt07_95_d7.5',
            'wilsonbg_givt07_95_t10_d5']
    t = M.groupby('rule')[['called', 'recall', 'precision', 'gain_vs_pct', 'gain_vs_wilson95']].mean()
    print('\n(5) mean of the three runs; gains are precision points at matched recall')
    print(t.reindex([k for k in keep if k in t.index]).assign(
        gain_vs_pct=lambda x: 100 * x['gain_vs_pct'], gain_vs_wilson95=lambda x: 100 * x['gain_vs_wilson95']).round(3).to_string())
    fam = M[M['family'].str.startswith('wilsonbg')].groupby(['family', 'run'])[['gain_vs_pct', 'gain_vs_wilson95']].mean()
    print('\nwilsonbg families, mean gain over their rules (points), per run:')
    print((100 * fam).round(2).unstack('run').to_string())

    # (6) where false calls fall, run 08_01 (recomputed with the same functions the sweep uses)
    hot = set(pl.loc[pl['hot'], 'kmer'])
    d, km, p0 = bg['08_01']
    km = km.view('S5').ravel().astype(str)
    calls = ws.rules(d, {'wilsonbg': p0})
    tr = pd.read_csv(f'{a.mock_dir}/truth.tsv', sep='\t')
    ids = {c: i for i, c in enumerate(fa.idx)}
    in_t = np.isin(co.keys(d['chrom'], d['pos'].to_numpy(), d['strand'], ids), co.keys(tr['chrom'], tr['pos'].to_numpy(), tr['strand'], ids))
    is_hot = np.isin(km, list(hot))
    print(f'\n(6) run 08_01: {len(d):,} sites with >= 20 reads, {is_hot.mean():.1%} in hot contexts; false calls by context')
    out = []
    for rname in ['pct_20', 'wilson95_t10', 'wilsonbg95_d5', 'wilsonbg95_d7.5', 'wilsonbg95_t10_d5']:
        c = calls[rname]
        fp = c & ~in_t
        out.append({'rule': rname, 'called': int(c.sum()), 'false': int(fp.sum()),
                    'false_in_hot': int((fp & is_hot).sum()), 'share_false_in_hot': (fp & is_hot).sum() / max(fp.sum(), 1),
                    'precision_hot': in_t[c & is_hot].mean() if (c & is_hot).any() else np.nan,
                    'precision_other': in_t[c & ~is_hot].mean() if (c & ~is_hot).any() else np.nan})
    print(pd.DataFrame(out).round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    m = sub.add_parser('make')
    m.add_argument('--out-dir', required=True)
    m.add_argument('--scenario', choices=['context', 'flat'], default='context')
    m.add_argument('--seed', type=int, default=1)
    m.add_argument('--contigs', type=int, default=4)
    m.add_argument('--length', type=int, default=60_000)
    s = sub.add_parser('summary')
    s.add_argument('--mock-dir', required=True)
    s.add_argument('--sweep-dir', required=True)
    a = ap.parse_args()
    if a.cmd == 'make':
        make(a)
    else:
        summary(a)


if __name__ == '__main__':
    main()
