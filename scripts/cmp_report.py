#!/usr/bin/env python3
"""cmp_report.py: gather the v6 vs v5.2 comparison outputs into report tables, summary figures and a markdown digest.

    python cmp_report.py --base /scratch/stein.an/basecall_cmp --out /projects/Genometechlab/andrew/basecall_cmp_results/report

Reads (missing pieces are reported in the digest, not fatal):
  ortho/<sample>_<arm>/tables/{summary,ci,stoichiometry}.tsv       accuracy against orthogonal sets (HEK293)
  plots/<sample>_<arm>/tables/pairwise_metrics.tsv                  site-level agreement between runs
  plots/HEK293_replicates_<run>_<arm>/tables/pairwise_metrics.tsv   agreement between HEK293 samples, per model
  confusion/<sample>/tables/confusion_<x>_vs_<y>.tsv                read-level call changes
  pileup064/<run>/<sample>.thresholds.tsv, sampleprobs_050/..., sp064_diag/*/thresholds.tsv
  ortho/<sample>_<arm>/tables/pr_curves.tsv                       precision-recall curves, operating points
  basecall_qc/<sample>/tables/{per_run,paired}.tsv                read identity and poly(A) (cmp_basecall_qc.py)
  motif/<sample>_<arm>/tables/motif.tsv                           m6A sequence context by class (cmp_motif.py)
  tower_ubam/tower_basecall.log                                   basecalling runtimes
Arms: nofilter (every call counted as its most likely class) and lab050 (thresholds estimated by modkit 0.5.0).
Writes to --out: tables/*.tsv, figures/*.{pdf,png}, digest.md
"""
import argparse
import glob
import os
import re
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HEK = ['08_01_25_RNA004_HEK293', '09_09_25_RNA004_HEK293', '10_06_25_RNA004_HEK293_CP']
SAMPLES = ['GM12878'] + HEK
SLAB = {'GM12878': 'GM12878', '08_01_25_RNA004_HEK293': 'HEK293 08/01', '09_09_25_RNA004_HEK293': 'HEK293 09/09',
        '10_06_25_RNA004_HEK293_CP': 'HEK293 10/06 CP'}
RUNS = ['A_d212_v6', 'C_d212_v52', 'B_d130_v52']
RLAB = {'A_d212_v6': 'v6 / dorado 2.1.2', 'C_d212_v52': 'v5.2 / dorado 2.1.2', 'B_d130_v52': 'v5.2 / dorado 1.3.0'}
ARMS = {'nofilter': 'threshold-free', 'lab050': 'lab rule (modkit 0.5.0 thresholds)'}
PAIRS = {'A_d212_v6 - C_d212_v52': 'v6 vs v5.2 (models)', 'C_d212_v52 - B_d130_v52': 'dorado 2.1.2 vs 1.3.0 (same models)',
         'A_d212_v6 - B_d130_v52': 'v6 vs production'}
KEY_SETS = [('m6A', 'GLORI1.0∩2.0'), ('Inosine', 'SLIC-seq_2of3'), ('Psi', 'BID∩PRAISE'), ('Psi', 'PRAISE'), ('m5C', 'UBS-seq')]
MODS = ['m6A', 'Inosine', '2OmeA', 'm5C', '2OmeC', 'Psi', '2OmeU', '2OmeG']
CODE = {'m6A': 'a', 'Inosine': '17596', '2OmeA': '69426', 'm5C': 'm', '2OmeC': '19228', 'Psi': '17802', '2OmeU': '19227', '2OmeG': '19229'}
BASE_OF = {'m6A': 'A', 'Inosine': 'A', '2OmeA': 'A', 'm5C': 'C', '2OmeC': 'C', 'Psi': 'T', '2OmeU': 'T', '2OmeG': 'G'}
COL = {'A_d212_v6': '#1b9e77', 'C_d212_v52': '#d95f02', 'B_d130_v52': '#7570b3'}
SCOL = {'08_01_25_RNA004_HEK293': '#1f78b4', '09_09_25_RNA004_HEK293': '#33a02c', '10_06_25_RNA004_HEK293_CP': '#e31a1c',
        'GM12878': '#6a3d9a'}
missing = []


def read(path, **kw):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        missing.append(path)
        return None
    return pd.read_csv(path, sep='\t', **kw)


def thresholds(path):
    """{base: value} from any modkit/cmp_pileup thresholds file (with or without a header)."""
    if not os.path.exists(path):
        missing.append(path)
        return {}
    out = {}
    for line in open(path):
        t = line.split()
        if t and t[0] in 'ACGT' and len(t[0]) == 1:
            try:
                out[t[0]] = float(t[-1])
            except ValueError:
                pass
    return out


def md_table(df, floatfmt='{:.3g}'):
    cols = list(df.columns)
    fmt = lambda v: floatfmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v) else ('' if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
    lines = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '---|' * len(cols)]
    lines += ['| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in df.itertuples(index=False)]
    return '\n'.join(lines)


def save(fig, out, name):
    fig.savefig(f'{out}/figures/{name}.pdf'); fig.savefig(f'{out}/figures/{name}.png', dpi=200); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--base', default='/scratch/stein.an/basecall_cmp')
    ap.add_argument('--out', required=True)
    ap.add_argument('--style')
    a = ap.parse_args()
    if a.style:
        plt.style.use(a.style)
    B = a.base
    os.makedirs(f'{a.out}/tables', exist_ok=True); os.makedirs(f'{a.out}/figures', exist_ok=True)
    md = ['# v6 vs v5.2 RNA004 modification calling: report tables', '',
          'Runs: ' + '; '.join(f'**{k}** = {v}' for k, v in RLAB.items()) + '. Arms: ' +
          '; '.join(f'**{k}** = {v}' for k, v in ARMS.items()) + '.', '']

    # ---------------------------------------------------------------- 1. accuracy against orthogonal sets
    acc, diffs, stoich, qc = [], [], [], None
    for arm in ARMS:
        for s in HEK:
            d = f'{B}/ortho/{s}_{arm}/tables'
            ci, sm, st = read(f'{d}/ci.tsv'), read(f'{d}/summary.tsv'), read(f'{d}/stoichiometry.tsv')
            if qc is None:
                qc = read(f'{d}/orthogonal_qc.tsv')
            if sm is not None:
                cov = sm.drop_duplicates(['set', 'mod']).set_index(['set', 'mod'])['sites_covered_all_runs']
            if ci is not None:
                ci = ci.assign(arm=arm, sample=s)
                one = ci[ci['comparison'].isin(RUNS)].copy()
                one['sites_covered'] = [cov.get((x, y), np.nan) for x, y in zip(one['set'], one['mod'])] if sm is not None else np.nan
                acc.append(one); diffs.append(ci[~ci['comparison'].isin(RUNS)])
            if st is not None:
                stoich.append(st.assign(arm=arm, sample=s))
    acc = pd.concat(acc, ignore_index=True) if acc else pd.DataFrame()
    diffs = pd.concat(diffs, ignore_index=True) if diffs else pd.DataFrame()
    stoich = pd.concat(stoich, ignore_index=True) if stoich else pd.DataFrame()
    for name, df in [('accuracy', acc), ('accuracy_differences', diffs), ('stoichiometry', stoich)]:
        df.to_csv(f'{a.out}/tables/{name}.tsv', sep='\t', index=False, float_format='%.4g')
    if qc is not None:
        qc.to_csv(f'{a.out}/tables/orthogonal_qc.tsv', sep='\t', index=False, float_format='%.4f')
        md += ['## Orthogonal sets after the reference check', '',
               md_table(qc[['set', 'mod', 'sites_listed', 'match_offset_0', 'status', 'sites_kept']], '{:.4f}'), '']

    if not diffs.empty:
        # replication verdict: per arm x set x metric, how many samples show v6 - v5.2 above / below zero with 95% CI
        v = diffs[diffs['comparison'] == 'A_d212_v6 - C_d212_v52'].copy()
        v['call'] = np.select([v['ci_low'] > 0, v['ci_high'] < 0], ['v6 higher', 'v6 lower'], 'no clear difference')
        rep = (v.groupby(['arm', 'mod', 'set', 'metric'])
                .agg(samples=('call', 'size'), v6_higher=('call', lambda c: int((c == 'v6 higher').sum())),
                     v6_lower=('call', lambda c: int((c == 'v6 lower').sum())), median_difference=('estimate', 'median'))
                .reset_index())                                    # works on pandas 1.x and 2.x
        rep.to_csv(f'{a.out}/tables/replication_v6_vs_v52.tsv', sep='\t', index=False, float_format='%.4g')
        allp = diffs.copy()
        allp['call'] = np.select([allp['ci_low'] > 0, allp['ci_high'] < 0], ['x higher', 'x lower'], 'no clear difference')
        rep_all = (allp.groupby(['comparison', 'arm', 'mod', 'set', 'metric'])
                   .agg(samples=('call', 'size'), x_higher=('call', lambda c: int((c == 'x higher').sum())),
                        x_lower=('call', lambda c: int((c == 'x lower').sum())), median_difference=('estimate', 'median'))
                   .reset_index())
        rep_all.to_csv(f'{a.out}/tables/replication_all_pairs.tsv', sep='\t', index=False, float_format='%.4g')
        keyrep = rep_all[rep_all.apply(lambda r: (r['mod'], r['set']) in KEY_SETS, axis=1)]
        for comp, label in PAIRS.items():
            blk = keyrep[keyrep['comparison'] == comp].drop(columns='comparison')
            if not blk.empty:
                md += [f'## Replication across HEK293 samples: {label} (x - y = {comp})',
                       'Counts of samples whose 95% CI for x - y lies above zero (x_higher) or below zero (x_lower).', '',
                       md_table(blk.assign(arm=blk['arm'].map(ARMS))), '']
        md += ['## Accuracy: v6 minus v5.2 (same dorado), per HEK293 sample, 95% paired block-bootstrap CI', '']
        for arm in ARMS:
            vv = v[(v['arm'] == arm) & v.apply(lambda r: (r['mod'], r['set']) in KEY_SETS, axis=1)].copy()
            if vv.empty:
                continue
            vv['sample'] = vv['sample'].map(SLAB)
            vv['difference [95% CI]'] = [f'{e:+.3f} [{lo:+.3f}, {hi:+.3f}]' for e, lo, hi in zip(vv['estimate'], vv['ci_low'], vv['ci_high'])]
            piv = vv.pivot_table(index=['mod', 'set', 'metric'], columns='sample', values='difference [95% CI]', aggfunc='first')
            md += [f'### {ARMS[arm]}', '', md_table(piv.reset_index()), '']
        if not acc.empty:
            aa = acc[acc.apply(lambda r: (r['mod'], r['set']) in KEY_SETS, axis=1)].copy()
            aa['value [95% CI]'] = [f'{e:.3f} [{lo:.3f}, {hi:.3f}]' for e, lo, hi in zip(aa['estimate'], aa['ci_low'], aa['ci_high'])]
            aa['run'] = aa['comparison'].map(RLAB); aa['sample'] = aa['sample'].map(SLAB)
            for arm in ARMS:
                x = aa[aa['arm'] == arm]
                if x.empty:
                    continue
                piv = x.pivot_table(index=['mod', 'set', 'metric', 'sample'], columns='run', values='value [95% CI]', aggfunc='first')
                md += [f'### Per-run values, {ARMS[arm]}', '', md_table(piv.reset_index()), '']

        # figure 1: v6 - v5.2 differences with CIs, per key set and metric, both arms, three samples
        metrics = ['recall_at_20', 'precision_at_20', 'auprc']
        fig, axes = plt.subplots(len(KEY_SETS), len(metrics), figsize=(12, 2.0 * len(KEY_SETS) + 0.8), squeeze=False, layout='constrained')
        for i, (mod, st_) in enumerate(KEY_SETS):
            for j, m in enumerate(metrics):
                ax = axes[i, j]; y = 0
                for arm, mk in zip(ARMS, ['o', 's']):
                    for s in HEK:
                        r = v[(v['arm'] == arm) & (v['set'] == st_) & (v['mod'] == mod) & (v['metric'] == m) & (v['sample'] == s)]
                        if not r.empty:
                            e, lo, hi = r.iloc[0][['estimate', 'ci_low', 'ci_high']]
                            ax.errorbar(e, y, xerr=[[e - lo], [hi - e]], fmt=mk, color=SCOL[s], ms=5, capsize=2,
                                        label=f'{SLAB[s]}, {ARMS[arm]}' if (i == 0 and j == 0) else None)
                        y += 1
                    y += 0.6
                ax.axvline(0, color='0.4', lw=0.8, ls='--'); ax.set_yticks([])
                if i == 0:
                    ax.set_title(m.replace('_', ' ').replace('at 20', 'at 20%'), fontsize=10)
                if j == 0:
                    ax.set_ylabel(f'{mod}\n{st_}', fontsize=9)
        fig.legend(loc='outside lower center', ncol=3, fontsize=7, frameon=False)
        fig.suptitle('v6 minus v5.2 (same dorado 2.1.2): difference with 95% paired block-bootstrap CI\n'
                     'circles: threshold-free; squares: lab rule (0.5.0 thresholds)', fontsize=11)
        save(fig, a.out, 'fig1_accuracy_differences')

    # ---------------------------------------------------------------- 2. read-level call changes
    conf = []
    for s in SAMPLES:
        for x, y in [('A_d212_v6', 'C_d212_v52'), ('C_d212_v52', 'B_d130_v52')]:
            t = read(f'{B}/confusion/{s}/tables/confusion_{x}_vs_{y}.tsv')
            if t is None:
                continue
            nx = t.groupby(['base', 'call_x'])['n'].sum(); ny = t.groupby(['base', 'call_y'])['n'].sum()
            both = t[t['call_x'] == t['call_y']].set_index(['base', 'call_x'])['n']
            tot = t.groupby('base')['n'].sum(); same = t[t['call_x'] == t['call_y']].groupby('base')['n'].sum()
            for (base, mod), n in both.items():
                if mod == 'canonical':
                    continue
                conf.append({'sample': s, 'pair': f'{x} vs {y}', 'base': base, 'mod': mod, 'calls_x': nx[(base, mod)],
                             'calls_y': ny.get((base, mod), np.nan), 'shared': n, 'ratio_x_over_y': nx[(base, mod)] / ny.get((base, mod), np.nan),
                             'shared_pct_of_x': 100 * n / nx[(base, mod)], 'shared_pct_of_y': 100 * n / ny.get((base, mod), np.nan),
                             'identical_calls_pct_this_base': 100 * same[base] / tot[base]})
    conf = pd.DataFrame(conf)
    conf.to_csv(f'{a.out}/tables/readlevel_changes.tsv', sep='\t', index=False, float_format='%.4g')
    if not conf.empty:
        md += ['## Read-level call changes (same reads, same positions, no threshold)', '']
        for pair in conf['pair'].unique():
            c = conf[conf['pair'] == pair]
            piv = c.pivot_table(index='mod', columns='sample', values='ratio_x_over_y').reindex([m for m in MODS if m in set(c['mod'])])
            piv.columns = [SLAB[x] for x in piv.columns]
            kept = c.pivot_table(index='mod', columns='sample', values='shared_pct_of_y').reindex(piv.index)
            kept.columns = [SLAB[x] for x in kept.columns]
            md += [f'### {pair}: calls_x / calls_y', '', md_table(piv.reset_index(), '{:.2f}'), '',
                   f'### {pair}: % of y\'s modification calls also made by x', '', md_table(kept.reset_index(), '{:.0f}'), '']
        fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
        c = conf[conf['pair'] == 'A_d212_v6 vs C_d212_v52']; f = conf[conf['pair'] == 'C_d212_v52 vs B_d130_v52']
        mods = [m for m in MODS if m in set(c['mod'])]; w = 0.8 / max(len(SAMPLES), 1)
        for k, s in enumerate(SAMPLES):
            r = c[c['sample'] == s].set_index('mod').reindex(mods)
            ax.bar(np.arange(len(mods)) + (k - (len(SAMPLES) - 1) / 2) * w, np.log2(r['ratio_x_over_y']), w, color=SCOL[s], label=SLAB[s])
            rf = f[f['sample'] == s].set_index('mod').reindex(mods)
            ax.plot(np.arange(len(mods)) + (k - (len(SAMPLES) - 1) / 2) * w, np.log2(rf['ratio_x_over_y']), 'k_', ms=8,
                    label='same models, dorado 2.1.2 vs 1.3.0' if k == 0 else None)
        ax.axhline(0, color='0.3', lw=0.8); ax.set_xticks(range(len(mods))); ax.set_xticklabels(mods)
        ax.set_ylabel('log2(v6 calls / v5.2 calls)'); ax.legend(fontsize=7, frameon=False, ncol=3)
        ax.set_title('Per-read modification calls on identical reads: v6 vs v5.2 (bars); dorado version only (black ticks)', fontsize=10)
        save(fig, a.out, 'fig2_readlevel_redistribution')

    # ---------------------------------------------------------------- 3. site-level agreement and reproducibility
    agr, reps = [], []
    for arm in ARMS:
        for s in SAMPLES:
            t = read(f'{B}/plots/{s}_{arm}/tables/pairwise_metrics.tsv')
            if t is not None:
                agr.append(t.assign(arm=arm, sample=s))
        for r in RUNS:
            t = read(f'{B}/plots/HEK293_replicates_{r}_{arm}/tables/pairwise_metrics.tsv')
            if t is not None:
                reps.append(t.assign(arm=arm, model=r))
    agr = pd.concat(agr, ignore_index=True) if agr else pd.DataFrame()
    reps = pd.concat(reps, ignore_index=True) if reps else pd.DataFrame()
    agr.to_csv(f'{a.out}/tables/site_agreement.tsv', sep='\t', index=False, float_format='%.4g')
    reps.to_csv(f'{a.out}/tables/replicate_agreement.tsv', sep='\t', index=False, float_format='%.4g')
    if not agr.empty:
        md += ['## Site-level agreement (sites with >= 20 reads in both runs)', '']
        for arm in ARMS:
            x = agr[agr['arm'] == arm]
            if x.empty:
                continue
            x = x.assign(pair=x['run_x'] + ' vs ' + x['run_y'], sample=x['sample'].map(SLAB))
            for metric in ['icc', 'jaccard', 'r_all']:
                if metric in x.columns:
                    piv = x.pivot_table(index=['pair', 'modification'], columns='sample', values=metric)
                    md += [f'### {ARMS[arm]}: {metric}', '', md_table(piv.reset_index(), '{:.2f}'), '']
    if not reps.empty:
        g = reps.groupby(['arm', 'model', 'modification'])[['icc', 'jaccard', 'r_all']].mean().reset_index()
        md += ['## Reproducibility across the three HEK293 samples (mean of the three pairwise comparisons)', '']
        for arm in ARMS:
            x = g[g['arm'] == arm]
            if not x.empty:
                piv = x.pivot_table(index='modification', columns='model', values='icc')
                piv.columns = [RLAB.get(c, c) for c in piv.columns]
                md += [f'### {ARMS[arm]}: ICC(2,1)', '', md_table(piv.reset_index(), '{:.2f}'), '']
        # figure 3: ICC by modification: v6 vs v5.2 within sample, dorado floor, replicates per model (threshold-free)
        arm = 'nofilter'
        x = agr[(agr['arm'] == arm) & agr['sample'].isin(HEK)] if not agr.empty else pd.DataFrame()
        gg = g[g['arm'] == arm]
        if not x.empty and not gg.empty:
            mods = [m for m in MODS if m in set(x['modification'])]
            series = {'v6 vs v5.2, same reads': x[(x['run_x'] == 'A_d212_v6') & (x['run_y'] == 'C_d212_v52')].groupby('modification')['icc'].mean(),
                      'dorado 2.1.2 vs 1.3.0, same reads': x[(x['run_x'] == 'C_d212_v52') & (x['run_y'] == 'B_d130_v52')].groupby('modification')['icc'].mean(),
                      'v6 across HEK293 samples': gg[gg['model'] == 'A_d212_v6'].set_index('modification')['icc'],
                      'v5.2 across HEK293 samples': gg[gg['model'] == 'C_d212_v52'].set_index('modification')['icc']}
            fig, ax = plt.subplots(figsize=(9, 4), layout='constrained'); w = 0.2
            for k, (lab, ser) in enumerate(series.items()):
                ax.bar(np.arange(len(mods)) + (k - 1.5) * w, ser.reindex(mods), w, label=lab,
                       color=['#1b9e77', '0.6', '#66c2a5', '#fc8d62'][k])
            ax.set_xticks(range(len(mods))); ax.set_xticklabels(mods); ax.set_ylim(0, 1.05); ax.set_ylabel('ICC(2,1), mean over HEK293 samples')
            ax.legend(fontsize=7, frameon=False, ncol=2)
            ax.set_title('Site-level agreement (threshold-free): model change vs dorado version vs resequencing', fontsize=10)
            save(fig, a.out, 'fig3_agreement_scale')

    # ---------------------------------------------------------------- 3b. precision-recall curves and operating points
    curves = []
    for arm in ARMS:
        for s in HEK:
            c = read(f'{B}/ortho/{s}_{arm}/tables/pr_curves.tsv')
            if c is not None:
                curves.append(c.assign(arm=arm, sample=s))
    curves = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    if not curves.empty:
        ops = []
        for (arm, s, st_, mod), g in curves.groupby(['arm', 'sample', 'set', 'mod']):
            ref = g[(g['run'] == 'C_d212_v52') & (g['cutoff_pct'] == 20)]
            v6 = g[g['run'] == 'A_d212_v6'].sort_values('cutoff_pct')
            if ref.empty or v6.empty:
                continue
            p20, r20, n20 = ref.iloc[0][['precision', 'recall', 'n_called']]
            v20 = v6[v6['cutoff_pct'] == 20].iloc[0]
            okp = v6[v6['precision'] >= p20]
            cp = okp.iloc[0] if not okp.empty else None
            cn = v6.iloc[int(np.abs(v6['n_called'].to_numpy() - n20).argmin())]
            ops.append({'arm': arm, 'sample': s, 'mod': mod, 'set': st_, 'v52_precision_20': p20, 'v52_recall_20': r20,
                        'v52_sites_20': n20, 'v6_precision_20': v20['precision'], 'v6_recall_20': v20['recall'],
                        'v6_sites_20': v20['n_called'],
                        'v6_cutoff_matching_v52_precision': cp['cutoff_pct'] if cp is not None else np.nan,
                        'v6_recall_there': cp['recall'] if cp is not None else np.nan,
                        'v6_sites_there': cp['n_called'] if cp is not None else np.nan,
                        'v6_cutoff_matching_v52_sites': cn['cutoff_pct'], 'v6_precision_at_matched_sites': cn['precision'],
                        'v6_recall_at_matched_sites': cn['recall']})
        ops = pd.DataFrame(ops)
        ops.to_csv(f'{a.out}/tables/operating_points.tsv', sep='\t', index=False, float_format='%.4g')
        k = ops[ops.apply(lambda r: (r['mod'], r['set']) in KEY_SETS, axis=1)] if not ops.empty else ops
        if not k.empty:
            md += ['## Operating points: what v6 cutoff matches v5.2 at its 20% cutoff?',
                   'v6_cutoff_matching_v52_precision: lowest v6 occupancy cutoff whose precision reaches v5.2\'s precision at 20%; '
                   'v6_cutoff_matching_v52_sites: the v6 cutoff calling about as many sites as v5.2 at 20%.', '',
                   md_table(k.assign(arm=k['arm'].map(ARMS), sample=k['sample'].map(SLAB))[
                       ['arm', 'sample', 'mod', 'set', 'v52_precision_20', 'v52_recall_20', 'v6_precision_20', 'v6_recall_20',
                        'v6_cutoff_matching_v52_precision', 'v6_recall_there', 'v6_cutoff_matching_v52_sites',
                        'v6_precision_at_matched_sites']]), '']
        fig, axes = plt.subplots(len(KEY_SETS), len(HEK), figsize=(4.0 * len(HEK), 3.0 * len(KEY_SETS)), squeeze=False, layout='constrained')
        for i, (mod, st_) in enumerate(KEY_SETS):
            for j, s in enumerate(HEK):
                ax = axes[i, j]
                for arm, ls in zip(ARMS, ['-', '--']):
                    for r in RUNS:
                        c = curves[(curves['arm'] == arm) & (curves['sample'] == s) & (curves['set'] == st_) &
                                   (curves['mod'] == mod) & (curves['run'] == r)].sort_values('cutoff_pct')
                        if c.empty:
                            continue
                        ax.plot(c['recall'], c['precision'], ls, color=COL[r], lw=1.3,
                                label=f'{RLAB[r]}, {ARMS[arm]}' if (i == 0 and j == 0) else None)
                        p = c[c['cutoff_pct'] == 20]
                        ax.plot(p['recall'], p['precision'], 'o', color=COL[r], ms=4, mec='k', mew=0.4)
                ax.set_xlim(0, 1); ax.set_ylim(bottom=0)
                if i == 0:
                    ax.set_title(SLAB[s], fontsize=10)
                if j == 0:
                    ax.set_ylabel(f'{mod} vs {st_}\nprecision (lower bound)', fontsize=8)
                if i == len(KEY_SETS) - 1:
                    ax.set_xlabel('recall')
        fig.legend(loc='outside lower center', ncol=3, fontsize=7, frameon=False)
        fig.suptitle('Precision-recall against orthogonal sets, occupancy cutoff 1-100% (dots: 20%)', fontsize=11)
        save(fig, a.out, 'fig5_pr_curves')

    if not stoich.empty:
        st = stoich[stoich['n_sites'] >= 10]
        sets = [x for x in ['GLORI1.0∩2.0', 'GLORI1.0', 'GLORI2.0', 'GLORI1.0_2022', 'BID-seq', 'PRAISE', 'BID∩PRAISE', 'UBS-seq'] if x in set(st['set'])]
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 0.55 * len(sets) * len(HEK) + 1.5), layout='constrained')
        y, ticks = 0, []
        for st_ in sets:
            for s in HEK:
                for k2, r in enumerate(RUNS):
                    for arm, mk in zip(ARMS, ['o', 's']):
                        row = st[(st['set'] == st_) & (st['sample'] == s) & (st['run'] == r) & (st['arm'] == arm)]
                        if row.empty:
                            continue
                        row = row.iloc[0]; yy = y + (k2 - 1) * 0.22
                        a1.errorbar(row['mean_diff_nanopore_minus_set'], yy, fmt=mk, color=COL[r], ms=4,
                                    xerr=[[row['mean_diff_nanopore_minus_set'] - row['loa_low']], [row['loa_high'] - row['mean_diff_nanopore_minus_set']]],
                                    elinewidth=0.6, alpha=0.9)
                        a2.plot(row['pearson_r'], yy, mk, color=COL[r], ms=4)
                ticks.append((y, f'{st_}, {SLAB[s]}')); y += 1
        for ax in (a1, a2):
            ax.set_yticks([t_[0] for t_ in ticks]); ax.set_yticklabels([t_[1] for t_ in ticks], fontsize=7); ax.invert_yaxis()
        hs = [plt.Line2D([], [], marker='o', ls='', color=COL[r], label=RLAB[r]) for r in RUNS] + \
             [plt.Line2D([], [], marker=mk, ls='', color='0.3', mfc='none', label=ARMS[arm]) for arm, mk in zip(ARMS, ['o', 's'])]
        a1.legend(handles=hs, fontsize=7, frameon=False, loc='lower left')
        a1.axvline(0, color='0.4', lw=0.8, ls='--'); a1.set_xlabel('mean occupancy difference, nanopore - set (points), with 95% limits of agreement')
        a2.set_xlabel('Pearson r, nanopore occupancy vs set level'); a2.set_yticklabels([])
        fig.suptitle('Stoichiometry at covered truth sites (colours: runs; circles threshold-free, squares lab rule)', fontsize=10)
        save(fig, a.out, 'fig6_stoichiometry')

    # ---------------------------------------------------------------- basecall accuracy, poly(A), motif context, runtimes
    bq_run, bq_pair = [], []
    for s in SAMPLES:
        x, y = read(f'{B}/basecall_qc/{s}/tables/per_run.tsv'), read(f'{B}/basecall_qc/{s}/tables/paired.tsv')
        if x is not None:
            bq_run.append(x)
        if y is not None:
            bq_pair.append(y)
        for f in ('identity', 'polya'):
            src = f'{B}/basecall_qc/{s}/figures/{f}.png'
            if os.path.exists(src):
                shutil.copy(src, f'{a.out}/figures/fig7_{f}_{s}.png')
    if bq_run:
        x = pd.concat(bq_run, ignore_index=True); x.to_csv(f'{a.out}/tables/basecall_accuracy.tsv', sep='\t', index=False, float_format='%.5g')
        md += ['## Basecall accuracy and poly(A), primary alignments (same read sample in every run)', '',
               md_table(x.assign(sample=x['sample'].map(SLAB), run=x['run'].map(RLAB))[
                   ['sample', 'run', 'reads', 'median_identity', 'median_mismatch_rate', 'median_insertion_rate',
                    'median_deletion_rate', 'polya_called_share', 'median_polya']], '{:.4g}'), '']
    if bq_pair:
        y = pd.concat(bq_pair, ignore_index=True); y.to_csv(f'{a.out}/tables/basecall_paired.tsv', sep='\t', index=False, float_format='%.5g')
        md += ['### Same reads, run x vs run y', '', md_table(y.assign(sample=y['sample'].map(SLAB))[
            ['sample', 'pair', 'shared_reads', 'median_identity_diff_x_minus_y', 'share_reads_x_higher', 'polya_pearson_r',
             'polya_median_diff_x_minus_y', 'polya_median_abs_diff']], '{:.4g}'), '']
    mot = []
    for arm in ARMS:
        for s in SAMPLES:
            x = read(f'{B}/motif/{s}_{arm}/tables/motif.tsv')
            if x is not None:
                mot.append(x.assign(arm=arm))
            src = f'{B}/motif/{s}_{arm}/figures/motif.png'
            if os.path.exists(src):
                shutil.copy(src, f'{a.out}/figures/fig8_motif_{s}_{arm}.png')
    if mot:
        x = pd.concat(mot, ignore_index=True); x.to_csv(f'{a.out}/tables/m6a_motif.tsv', sep='\t', index=False, float_format='%.4g')
        md += ['## m6A sites by model: DRACH context and GLORI support', '',
               md_table(x.assign(arm=x['arm'].map(ARMS), sample=x['sample'].map(SLAB))[
                   ['arm', 'sample', 'class', 'sites', 'drach_share', 'glori_share', 'median_pct_v6', 'median_pct_v52']], '{:.3g}'), '']
    logf = f'{B}/tower_ubam/tower_basecall.log'
    if os.path.exists(logf):
        rt = [m_.groups() for m_ in re.finditer(r'done\s+(\S+)\s+(\S+)\s+in\s+(\d+)\s+min\s+\(([^)]+)\)', open(logf).read())]
        if rt:
            rt = pd.DataFrame(rt, columns=['run', 'sample', 'minutes', 'ubam_size']).astype({'minutes': int})
            rt.to_csv(f'{a.out}/tables/runtimes.tsv', sep='\t', index=False)
            piv = rt.pivot_table(index='sample', columns='run', values='minutes').rename(columns=RLAB)
            piv.index = [SLAB.get(i, i) for i in piv.index]
            batch = [l.strip() for l in open(logf) if 'batch' in l.lower()][:3]
            md += ['## Basecalling time on the tower (4 x A100 80 GB, minutes per sample subset)', '',
                   md_table(piv.reset_index().rename(columns={'index': 'sample'}), '{:.0f}'), '',
                   'Batch-size notes from the log: ' + ' | '.join(batch), '']
    else:
        missing.append(logf)

    # ---------------------------------------------------------------- 4. thresholds
    rows = []
    for s in SAMPLES:
        for r in RUNS:
            t64 = thresholds(f'{B}/pileup064/{r}/{s}.thresholds.tsv')
            t50 = thresholds(f'{B}/sampleprobs_050/{r}/thresholds.tsv' if s == 'GM12878' else f'{B}/sampleprobs_050/{r}/{s}/thresholds.tsv')
            for base in 'ACGT':
                rows.append({'sample': s, 'run': r, 'base': 'U' if base == 'T' else base,
                             'modkit_0.6.4_as_run': t64.get(base, np.nan), 'modkit_0.5.0_lab_method': t50.get(base, np.nan)})
    thr = pd.DataFrame(rows)
    for name in ['like050', 'nosampling']:
        t = thresholds(f'{B}/sp064_diag/{name}/thresholds.tsv')
        for base, val in t.items():
            thr = pd.concat([thr, pd.DataFrame([{'sample': 'GM12878', 'run': f'B_d130_v52 (0.6.4 {name})', 'base': 'U' if base == 'T' else base,
                                                 'modkit_0.6.4_as_run': val}])], ignore_index=True)
    thr.to_csv(f'{a.out}/tables/thresholds.tsv', sep='\t', index=False, float_format='%.4f')
    md += ['## Per-base 10th-percentile thresholds', '', md_table(thr.assign(sample=thr['sample'].map(lambda x: SLAB.get(x, x))), '{:.3f}'), '']
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharey=True, layout='constrained')
    for ax, base in zip(axes, 'ACGU'):
        t = thr[(thr['base'] == base) & thr['run'].isin(RUNS)]
        for k, r in enumerate(RUNS):
            u = t[t['run'] == r]
            ax.scatter(np.full(len(u), k - 0.12), u['modkit_0.6.4_as_run'], color=COL[r], marker='o', s=22, label='0.6.4 as run' if (k == 0 and base == 'A') else None)
            ax.scatter(np.full(len(u), k + 0.12), u['modkit_0.5.0_lab_method'], color=COL[r], marker='^', s=26, facecolors='none',
                       label='0.5.0 (lab method)' if (k == 0 and base == 'A') else None)
        ax.set_xticks(range(len(RUNS))); ax.set_xticklabels([RLAB[r].replace(' / ', '\n') for r in RUNS], fontsize=7)
        ax.set_title(base); ax.axhline(0.998, color='0.5', lw=0.7, ls=':')
    axes[0].set_ylabel('threshold (10th percentile)'); axes[0].legend(fontsize=7, frameon=False, loc='lower left')
    fig.suptitle('Per-base filter thresholds: modkit 0.6.4 estimate as run (circles) vs modkit 0.5.0 (triangles); dotted = 0.998', fontsize=10)
    save(fig, a.out, 'fig4_thresholds')

    md += ['## Missing inputs', ''] + ([f'- {m}' for m in sorted(set(missing))] or ['none']) + ['']
    open(f'{a.out}/digest.md', 'w').write('\n'.join(md))
    print(f'wrote {a.out}/digest.md, tables/ and figures/; {len(set(missing))} missing inputs', file=sys.stderr)


if __name__ == '__main__':
    main()
