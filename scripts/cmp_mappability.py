#!/usr/bin/env python3
"""cmp_mappability.py: do calls in hard-to-map places perform worse against the orthogonal maps, and for which model?

    python cmp_mappability.py --mod m6A --model v6=A_d212_v6 --model v5.2=C_d212_v52 \\
        --sample 08_01_25_RNA004_HEK293 --sample 09_09_25_RNA004_HEK293 --sample 10_06_25_RNA004_HEK293_CP \\
        --all-dir <pileup064_nofilter> --q40-dir <pileup064_nofilter_q40> \\
        --strata 'GIAB low mappability=GRCh38_lowmappabilityall.bed.gz' --strata 'GIAB segdups=GRCh38_segdups.bed.gz' \\
        --ortho-dir /projects/Genometechlab/andrew/orthogonal --fasta GRCh38.fa --mask-bed clusters.bed --out-dir mapp/m6A

Two lenses, per site:
  unique-read fraction  valid reads in the MAPQ >= 40 pileup / valid reads in the all-primary pileup, same model and
                        sample (<all-dir>/<run>/<sample>.bed and <q40-dir>/<run>/<sample>.bed). Strata: unique
                        (>= 0.95), partly multimapped (0.5-0.95), mostly multimapped (< 0.5). This asks whether
                        multimapping reads corrupt the nanopore call.
  static regions        membership of each --strata BED (e.g. GIAB GRCh38 low mappability, segmental duplications).
                        These mark where short reads cannot place uniquely, i.e. where the short-read truth maps
                        are blind, whatever the nanopore reads do.
Sites: >= --min-valid valid reads (all-primary pileup); calls: percent modified >= --pct, made two ways,
  all reads  from the all-primary pileup (what the comparison used so far)
  MAPQ >= 40 from the filtered pileup, at sites with >= --min-valid filtered valid reads (the pipeline's new default)
Per sample, model, read set and stratum: covered positions, truth sites per 1,000 covered positions (how blind the
map is there), calls, precision (a lower bound everywhere, and lower still where the map is blind), recall, and the
median occupancy shift from dropping MAPQ < 40 reads at called sites. Because map blindness lowers every model's
precision alike, compare models WITHIN a stratum: gaps.tsv has first --model minus second --model per stratum.
Outputs in --out-dir: strata.tsv, gaps.tsv, mappability.pdf/png. Needs cmp_rules.py and cmp_ortho.py alongside.
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
import cmp_rules as cr  # noqa: E402

U_BINS = [-0.01, 0.5, 0.95, 1.01]
U_LABELS = ['mostly multimapped (<50% unique)', 'partly multimapped (50-95%)', 'unique (>=95%)']


def in_intervals(chrom, start0, bed):
    """True where 0-based start lies in any interval of bed (chrom, start, end)."""
    hit = np.zeros(len(chrom), dtype=bool)
    for c, iv in bed.groupby('chrom'):
        sel = np.flatnonzero(chrom == c)
        if not len(sel):
            continue
        s, e = iv['start'].to_numpy(), iv['end'].to_numpy()
        o = np.argsort(s); s, e = s[o], np.maximum.accumulate(e[o])
        z = start0[sel]
        j = np.searchsorted(s, z, side='right') - 1
        hit[sel] = (j >= 0) & (z < e[np.clip(j, 0, None)])
    return hit


def read_bed3(path):
    return pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'], comment='#',
                       dtype={'chrom': str})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--mod', required=True, choices=sorted(cr.TRUTH))
    ap.add_argument('--model', action='append', required=True, metavar='LABEL=RUNDIR',
                    help='two of these: the first minus the second goes in gaps.tsv')
    ap.add_argument('--sample', action='append', required=True)
    ap.add_argument('--all-dir', required=True)
    ap.add_argument('--q40-dir', required=True)
    ap.add_argument('--strata', action='append', default=[], metavar='NAME=BED')
    ap.add_argument('--ortho-dir')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--truth-tsv', help='testing: chrom, pos, strand')
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-valid', type=int, default=20)
    ap.add_argument('--pct', type=float, default=20)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    models = cr.pairs(a.model)
    strata = [(n, read_bed3(p)) for n, p in cr.pairs(a.strata)]
    st = cr.TRUTH[a.mod]
    if a.truth_tsv:
        truth = pd.read_csv(a.truth_tsv, sep='\t')
    else:
        listed, _ = co.load_sets(a.ortho_dir, '10ng_mRNA_input')
        listed = listed[listed['mod'] == a.mod]
        truth, _ = co.check_against_reference(listed, co.Fasta(a.fasta))
        truth = co.add_consensus(truth)
        truth = truth[truth['set'] == st]
    mask = read_bed3(a.mask_bed) if a.mask_bed else None
    ids = {}
    for line in open(a.fasta + '.fai'):
        ids[line.split('\t')[0]] = len(ids)
    tk = np.unique(co.keys(truth['chrom'], truth['pos'].to_numpy(), truth['strand'], ids))
    code = co.CODE[a.mod]

    rows = []
    for sample in a.sample:
        for label, run in models:
            d = cr.read_pileup(f'{a.all_dir}/{run}/{sample}.bed', code, a.min_valid, a.contig_regex, mask)
            q = cr.read_pileup(f'{a.q40_dir}/{run}/{sample}.bed', code, 1, a.contig_regex, mask)
            d['key'] = co.keys(d['chrom'], d['pos'].to_numpy(), d['strand'], ids)
            q['key'] = co.keys(q['chrom'], q['pos'].to_numpy(), q['strand'], ids)
            q = q.drop_duplicates('key').set_index('key')
            d['nvalid_q'] = q['nvalid'].reindex(d['key']).fillna(0).to_numpy()
            d['pct_q'] = q['pct'].reindex(d['key']).to_numpy()
            d['ufrac'] = (d['nvalid_q'] / d['nvalid']).clip(upper=1)
            d['in_truth'] = np.isin(d['key'], tk)
            d['u_stratum'] = pd.cut(d['ufrac'], U_BINS, labels=U_LABELS).astype(str)
            groupings = [('unique-read fraction', d['u_stratum'].to_numpy())]
            for name, bed in strata:
                inside = in_intervals(d['chrom'].to_numpy(), d['start'].to_numpy(), bed)
                groupings.append((name, np.where(inside, f'in {name}', f'outside {name}')))
            groupings.append(('all sites', np.full(len(d), 'all sites', dtype=object)))
            calls = {'all reads': (np.ones(len(d), bool), d['pct'].to_numpy() >= a.pct),
                     'MAPQ >= 40': (d['nvalid_q'].to_numpy() >= a.min_valid, d['pct_q'].fillna(-1).to_numpy() >= a.pct)}
            for reads, (universe, called) in calls.items():
                called = called & universe
                for gname, lab in groupings:
                    for s in pd.unique(lab):
                        sel = (lab == s) & universe
                        cs, ts = called & sel, d['in_truth'].to_numpy() & sel
                        shift = (d['pct'] - d['pct_q'])[cs]
                        rows.append({'mod': a.mod, 'sample': sample, 'model': label, 'reads': reads, 'grouping': gname,
                                     'stratum': s, 'covered': int(sel.sum()), 'truth_covered': int(ts.sum()),
                                     'truth_per_1000_covered': 1000 * ts.sum() / max(sel.sum(), 1),
                                     'called': int(cs.sum()),
                                     'precision': d['in_truth'][cs].mean() if cs.any() else np.nan,
                                     'recall': called[ts].mean() if ts.any() else np.nan,
                                     'median_pct_shift_at_calls': float(np.nanmedian(shift)) if cs.any() else np.nan})
            print(f'{sample} {label}: {len(d):,} sites; unique-read fraction strata: '
                  + ', '.join(f'{k.split(" (")[0]} {v:,}' for k, v in d['u_stratum'].value_counts().items()), file=sys.stderr)
    S = pd.DataFrame(rows)
    S.to_csv(f'{a.out_dir}/strata.tsv', sep='\t', index=False, float_format='%.4g')

    (m1, _), (m2, _) = models[0], models[1]
    k = ['sample', 'reads', 'grouping', 'stratum']
    g = S[S['model'] == m1].merge(S[S['model'] == m2], on=k, suffixes=(f'_{m1}', f'_{m2}'))
    g['precision_gap'] = g[f'precision_{m1}'] - g[f'precision_{m2}']
    g['recall_gap'] = g[f'recall_{m1}'] - g[f'recall_{m2}']
    G = g[k + ['precision_gap', 'recall_gap', f'called_{m1}', f'called_{m2}', f'truth_per_1000_covered_{m1}']]
    G.to_csv(f'{a.out_dir}/gaps.tsv', sep='\t', index=False, float_format='%.4g')

    groups = [x for x in S['grouping'].unique() if x != 'all sites']
    fig, axes = plt.subplots(len(groups), 2, figsize=(8.5, 2.6 * len(groups)), layout='constrained', squeeze=False)
    cols = {m1: '#d95f02', m2: '#1b9e77'}
    for i, gname in enumerate(groups):
        sub = S[(S['grouping'] == gname)]
        order = [s for s in (U_LABELS if gname == 'unique-read fraction' else sorted(sub['stratum'].unique())) if s in set(sub['stratum'])]
        ax = axes[i, 0]
        w = 0.2
        for j, (mdl, reads) in enumerate([(m1, 'all reads'), (m1, 'MAPQ >= 40'), (m2, 'all reads'), (m2, 'MAPQ >= 40')]):
            v = sub[(sub['model'] == mdl) & (sub['reads'] == reads)].groupby('stratum')['precision'].agg(['mean', 'min', 'max']).reindex(order)
            x = np.arange(len(order)) + (j - 1.5) * w
            ax.bar(x, v['mean'], w, color=cols[mdl], alpha=1 if reads == 'all reads' else 0.5, label=f'{mdl}, {reads}')
            ax.errorbar(x, v['mean'], yerr=[np.clip(v['mean'] - v['min'], 0, None), np.clip(v['max'] - v['mean'], 0, None)],
                        fmt='none', ecolor='k', lw=0.5)
        dens = sub[(sub['model'] == m1) & (sub['reads'] == 'all reads')].groupby('stratum')['truth_per_1000_covered'].mean().reindex(order)
        ax.set_xticks(np.arange(len(order))); ax.set_xticklabels([f'{s.split(" (")[0]}\n{t:.1f} truth/1k' for s, t in zip(order, dens)], fontsize=6.5)
        ax.set_ylabel('precision'); ax.set_title(f'{a.mod} against {st}: by {gname}', fontsize=8)
        if i == 0:
            ax.legend(frameon=False, fontsize=6)
        ax = axes[i, 1]
        gg = G[G['grouping'] == gname]
        for j, reads in enumerate(['all reads', 'MAPQ >= 40']):
            v = gg[gg['reads'] == reads].groupby('stratum')['precision_gap'].agg(['mean', 'min', 'max']).reindex(order) * 100
            x = np.arange(len(order)) + (j - 0.5) * 0.35
            ax.bar(x, v['mean'], 0.35, color='0.3' if reads == 'all reads' else '0.65', label=reads)
            ax.errorbar(x, v['mean'], yerr=[np.clip(v['mean'] - v['min'], 0, None), np.clip(v['max'] - v['mean'], 0, None)],
                        fmt='none', ecolor='k', lw=0.5)
        ax.axhline(0, color='k', lw=0.5)
        ax.set_xticks(np.arange(len(order))); ax.set_xticklabels([s.split(' (')[0] for s in order], fontsize=6.5)
        ax.set_ylabel(f'precision, {m1} minus {m2} (points)'); ax.set_title('Model gap within each stratum', fontsize=8)
        if i == 0:
            ax.legend(frameon=False, fontsize=6)
    fig.savefig(f'{a.out_dir}/mappability.pdf'); fig.savefig(f'{a.out_dir}/mappability.png', dpi=200); plt.close(fig)
    with pd.option_context('display.width', 220, 'display.max_rows', 100):
        t = S[S['reads'] == 'all reads'].groupby(['grouping', 'stratum', 'model'])[['covered', 'truth_per_1000_covered', 'called', 'precision', 'recall', 'median_pct_shift_at_calls']].mean()
        print(t.round(3).to_string(), file=sys.stderr)
        print(G.groupby(['reads', 'grouping', 'stratum'])[['precision_gap', 'recall_gap']].mean().round(3).to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()
