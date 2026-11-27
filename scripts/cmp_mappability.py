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

Options added 27 Nov (without them every output is byte-identical to the version before, md5 f2a86999):
  --partition NAME=BED  a labelled BED (chrom, start, end, label[, score, strand]; no strand or '.' = both strands)
                        that gives each site one class, e.g. transcript regions from make_transcript_partition.py.
                        Overlapping intervals on a strand are an error; sites in no interval are 'unlabelled'. A
                        '#labels<TAB>a,b,...' header line sets the figure's order.
  --truth-set SET       any cmp_ortho set of the modification (default: its main set, cmp_rules.TRUTH).
  --calls-all-dir, --calls-q40-dir
                        the pileups the calls come from (both or neither; default --all-dir and --q40-dir), e.g. the
                        fixed-0.70 pair for definition 2. The unique-read fraction always comes from --all-dir and
                        --q40-dir, which must then be the no-filter pair (fixed threshold 0, checked in their
                        <sample>.summary.txt, as is MAPQ >= 40 for the q40 dirs only); a site with no valid read in
                        the no-filter all-reads pileup goes in stratum 'no no-filter coverage'.
  --fail-as-canonical   for the calls pileups: depth N_valid + N_fail, percent modified 100 N_mod / (N_valid + N_fail)
                        (cmp_rules.read_pileup), for the --min-valid floor, both read sets and the occupancy shift.
  --calls-threshold X   stop unless every calls pileup's summary.txt shows fixed threshold X for A, C, G and T.
  --shared-strata       also group by a unique-read stratum shared by the two --model runs: per sample, the mean
                        of the two models' unique-read fractions (each from its own no-filter pair, at every site
                        either model covers), so both models' calls are split on the same positions ('not covered by
                        both models' where one model's all-reads pileup has no valid read there). The per-model
                        grouping stays. stratum_moves.tsv: per sample, sites in both models' universes (all reads),
                        how many have a different own stratum in the two models, and the share of moving and staying
                        sites inside each --strata BED (enrichment = moving share / staying share);
                        stratum_crosstab.tsv: own stratum of the first model x the second.
  With any of these, settings.tsv records the settings.
"""
import argparse
import gzip
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
import cmp_rules as cr  # noqa: E402

U_BINS = [-0.01, 0.5, 0.95, 1.01]
U_LABELS = ['mostly multimapped (<50% unique)', 'partly multimapped (50-95%)', 'unique (>=95%)']
NO_NF = 'no no-filter coverage'
SHARED_NONE = 'not covered by both models'


def ufrac_at(n_all, n_q, keys):
    """unique-read fraction at keys from one model's no-filter pair (Series of N_valid by key): MAPQ-40 / all valid
    reads, capped at 1; NaN where the all-reads pileup has no valid read."""
    na = n_all.reindex(keys).to_numpy(dtype=float)
    nq = n_q.reindex(keys).fillna(0).to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.clip(np.where(na > 0, nq / na, np.nan), None, 1)


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


def read_partition(path):
    """labelled BED -> ({(chrom, strand): (starts, ends, labels)}, display order or None). One class per position:
    overlapping intervals on a strand ('.' counts for both) are an error."""
    order = None
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('#labels\t'):
                order = line.rstrip('\n').split('\t', 1)[1].split(',')
    b = pd.read_csv(path, sep='\t', header=None, comment='#', dtype={0: str, 3: str})
    if b.shape[1] < 4:
        sys.exit(f'--partition {path}: needs a label in column 4')
    strand = b[5].astype(str).to_numpy() if b.shape[1] >= 6 else np.full(len(b), '.', dtype=object)
    bad = sorted(set(strand) - {'+', '-', '.'})
    if bad:
        sys.exit(f'--partition {path}: strand values {bad}')
    part = {}
    for s in '+-':
        x = b[(strand == s) | (strand == '.')]
        for c, iv in x.groupby(0):
            iv = iv.sort_values(1)
            st, en = iv[1].to_numpy(np.int64), iv[2].to_numpy(np.int64)
            if (st[1:] < en[:-1]).any():
                k = int(np.flatnonzero(st[1:] < en[:-1])[0])
                sys.exit(f'--partition {path}: overlapping intervals on {c} {s} at {st[k]}-{en[k]} and '
                         f'{st[k + 1]}-{en[k + 1]}: one class per site is required')
            part[(c, s)] = (st, en, iv[3].to_numpy(dtype=object))
    return part, order


def label_sites(chrom, strand, start0, part):
    """each site's class in a partition (0-based start, strand-aware); 'unlabelled' where no interval holds it."""
    lab = np.full(len(chrom), 'unlabelled', dtype=object)
    for (c, s), sel in pd.DataFrame({'c': chrom, 's': strand}).groupby(['c', 's']).indices.items():
        if (c, s) not in part:
            continue
        st, en, lb = part[(c, s)]
        z = start0[sel]
        j = np.searchsorted(st, z, side='right') - 1
        ok = (j >= 0) & (z < en[np.clip(j, 0, None)])
        lab[sel[ok]] = lb[j[ok]]
    return lab


def pileup_settings(bed):
    """(thresholds {base: value} or None if not fixed, MAPQ floor or None) from <prefix>.summary.txt of cmp_pileup.sh."""
    p = bed[:-4] + '.summary.txt' if bed.endswith('.bed') else bed + '.summary.txt'
    if not os.path.exists(p):
        sys.exit(f'{p} missing: cannot check how {bed} was made')
    txt = open(p).read()
    m = re.search(r'thresholds \(fixed[^)]*\): (.*?) in ', txt)
    th = dict(x.split(':') for x in m.group(1).split()) if m else None
    q = re.search(r'keeping primary alignments with MAPQ >= (\d+)', txt)
    return ({k: float(v) for k, v in th.items()} if th else None), (int(q.group(1)) if q else None)


def check_pileup(bed, want_th, want_mapq, what):
    th, mapq = pileup_settings(bed)
    if th is None or set(th) != set('ACGT') or any(abs(v - want_th) > 1e-9 for v in th.values()):
        sys.exit(f'{bed}: {what} needs fixed threshold {want_th:g} for A, C, G, T; its summary says {th}')
    if mapq != want_mapq:
        sys.exit(f'{bed}: {what} needs MAPQ filter {want_mapq}; its summary says {mapq}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--mod', required=True, choices=sorted(cr.TRUTH))
    ap.add_argument('--model', action='append', required=True, metavar='LABEL=RUNDIR',
                    help='two of these: the first minus the second goes in gaps.tsv')
    ap.add_argument('--sample', action='append', required=True)
    ap.add_argument('--all-dir', required=True)
    ap.add_argument('--q40-dir', required=True)
    ap.add_argument('--strata', action='append', default=[], metavar='NAME=BED')
    ap.add_argument('--partition', action='append', default=[], metavar='NAME=BED',
                    help='labelled BED, one class per site (strand-aware when column 6 is + or -)')
    ap.add_argument('--truth-set', help="cmp_ortho set to score against (default: the modification's main set)")
    ap.add_argument('--calls-all-dir', help='pileups the all-reads calls come from (default --all-dir)')
    ap.add_argument('--calls-q40-dir', help='pileups the MAPQ >= 40 calls come from (default --q40-dir)')
    ap.add_argument('--fail-as-canonical', action='store_true',
                    help='calls pileups: depth N_valid + N_fail, occupancy N_mod / (N_valid + N_fail)')
    ap.add_argument('--shared-strata', action='store_true',
                    help='also a unique-read stratum shared by the two models (mean of their fractions); writes '
                         'stratum_moves.tsv and stratum_crosstab.tsv')
    ap.add_argument('--calls-threshold', type=float,
                    help='check that the calls pileups were made at this fixed threshold for every base')
    ap.add_argument('--ortho-dir')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--truth-tsv', help='testing: chrom, pos, strand')
    ap.add_argument('--mask-bed')
    ap.add_argument('--contig-regex', default=r'_random|_alt|_fix|_decoy|chrUn_|EBV|HLA-')
    ap.add_argument('--min-valid', type=int, default=20)
    ap.add_argument('--pct', type=float, default=20)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    if bool(a.calls_all_dir) != bool(a.calls_q40_dir):
        ap.error('--calls-all-dir and --calls-q40-dir go together')
    new_calls = bool(a.calls_all_dir or a.fail_as_canonical or a.calls_threshold is not None)
    calls_all, calls_q40 = a.calls_all_dir or a.all_dir, a.calls_q40_dir or a.q40_dir
    os.makedirs(a.out_dir, exist_ok=True)
    models = cr.pairs(a.model)
    strata = [(n, read_bed3(p)) for n, p in cr.pairs(a.strata)]
    partitions = [(n, *read_partition(p)) for n, p in cr.pairs(a.partition)]
    st = a.truth_set or cr.TRUTH[a.mod]
    if a.truth_tsv:
        truth = pd.read_csv(a.truth_tsv, sep='\t')
    else:
        listed, _ = co.load_sets(a.ortho_dir, '10ng_mRNA_input')
        listed = listed[listed['mod'] == a.mod]
        truth, _ = co.check_against_reference(listed, co.Fasta(a.fasta))
        truth = co.add_consensus(truth)
        sets = sorted(truth['set'].unique())
        truth = truth[truth['set'] == st]
        if not len(truth):
            sys.exit(f'no {a.mod} sites in truth set {st}; sets: {", ".join(sets)}')
    mask = read_bed3(a.mask_bed) if a.mask_bed else None
    ids = {}
    for line in open(a.fasta + '.fai'):
        ids[line.split('\t')[0]] = len(ids)
    tk = np.unique(co.keys(truth['chrom'], truth['pos'].to_numpy(), truth['strand'], ids))
    code = co.CODE[a.mod]

    rows, moves, cross = [], [], []
    for sample in a.sample:
        prep, nf = {}, {}
        for label, run in models:
            if new_calls:
                for path, th, mq, what in [(f'{a.all_dir}/{run}/{sample}.bed', 0, None, 'the no-filter all-reads pileup'),
                                           (f'{a.q40_dir}/{run}/{sample}.bed', 0, 40, 'the no-filter MAPQ-40 pileup')]:
                    check_pileup(path, th, mq, what)
                if a.calls_threshold is not None:
                    check_pileup(f'{calls_all}/{run}/{sample}.bed', a.calls_threshold, None, 'the all-reads calls pileup')
                    check_pileup(f'{calls_q40}/{run}/{sample}.bed', a.calls_threshold, 40, 'the MAPQ-40 calls pileup')
            d = cr.read_pileup(f'{calls_all}/{run}/{sample}.bed', code, a.min_valid, a.contig_regex, mask,
                               a.fail_as_canonical)
            q = cr.read_pileup(f'{calls_q40}/{run}/{sample}.bed', code, 1, a.contig_regex, mask, a.fail_as_canonical)
            d['key'] = co.keys(d['chrom'], d['pos'].to_numpy(), d['strand'], ids)
            q['key'] = co.keys(q['chrom'], q['pos'].to_numpy(), q['strand'], ids)
            q = q.drop_duplicates('key').set_index('key')
            d['nvalid_q'] = q['nvalid'].reindex(d['key']).fillna(0).to_numpy()
            d['pct_q'] = q['pct'].reindex(d['key']).to_numpy()
            if not new_calls:
                d['ufrac'] = (d['nvalid_q'] / d['nvalid']).clip(upper=1)
            if new_calls or a.shared_strata:        # unique-read fraction from the no-filter pair, at any site
                n_all = cr.read_pileup(f'{a.all_dir}/{run}/{sample}.bed', code, 1, a.contig_regex, mask)
                n_q = cr.read_pileup(f'{a.q40_dir}/{run}/{sample}.bed', code, 1, a.contig_regex, mask)
                n_all = n_all.set_index(co.keys(n_all['chrom'], n_all['pos'].to_numpy(), n_all['strand'], ids))['nvalid']
                n_q = n_q.set_index(co.keys(n_q['chrom'], n_q['pos'].to_numpy(), n_q['strand'], ids))['nvalid']
                n_all, n_q = n_all[~n_all.index.duplicated()], n_q[~n_q.index.duplicated()]
                if new_calls:
                    d['ufrac'] = ufrac_at(n_all, n_q, d['key'])
                if a.shared_strata:
                    nf[label] = (n_all, n_q)
                del n_all, n_q
            d['in_truth'] = np.isin(d['key'], tk)
            d['u_stratum'] = pd.cut(d['ufrac'], U_BINS, labels=U_LABELS).astype(str)
            if new_calls:
                d.loc[d['ufrac'].isna(), 'u_stratum'] = NO_NF
            prep[label] = d
        if a.shared_strata:                         # one stratum per site: mean of the models' fractions
            for label, _ in models:
                d = prep[label]
                u = np.column_stack([ufrac_at(*nf[m], d['key']) for m, _ in models])
                us = pd.cut(pd.Series(u.mean(axis=1)), U_BINS, labels=U_LABELS).astype(str).to_numpy()
                us[np.isnan(u).any(axis=1)] = SHARED_NONE
                d['u_shared'] = us
            (m1_, _), (m2_, _) = models[0], models[1]
            x_, y_ = prep[m1_], prep[m2_]
            both, ix, iy = np.intersect1d(x_['key'].to_numpy(), y_['key'].to_numpy(), return_indices=True)
            sx, sy = x_['u_stratum'].to_numpy()[ix], y_['u_stratum'].to_numpy()[iy]
            ct = pd.crosstab(pd.Series(sx, name=f'stratum_{m1_}'), pd.Series(sy, name=f'stratum_{m2_}'))
            for (r1, r2), n in ct.stack().items():
                cross.append({'sample': sample, f'stratum_{m1_}': r1, f'stratum_{m2_}': r2, 'sites': int(n)})
            moving = sx != sy
            mv = {'sample': sample, 'sites_in_both_universes': len(both), 'moving': int(moving.sum()),
                  'moving_frac': moving.mean() if len(both) else np.nan}
            ch_, st_ = x_['chrom'].to_numpy()[ix], x_['start'].to_numpy()[ix]
            for name, bed in strata:
                inside = in_intervals(ch_, st_, bed)
                fm = inside[moving].mean() if moving.any() else np.nan
                fs = inside[~moving].mean() if (~moving).any() else np.nan
                mv.update({f'moving_in {name}': int(inside[moving].sum()), f'moving_frac_in {name}': fm,
                           f'staying_frac_in {name}': fs, f'enrichment {name}': fm / fs if fs else np.nan})
            moves.append(mv)
            del nf
        for label, run in models:
            d = prep.pop(label)
            groupings = [('unique-read fraction', d['u_stratum'].to_numpy())]
            if a.shared_strata:
                groupings.append(('unique-read fraction (shared)', d['u_shared'].to_numpy()))
            for name, bed in strata:
                inside = in_intervals(d['chrom'].to_numpy(), d['start'].to_numpy(), bed)
                groupings.append((name, np.where(inside, f'in {name}', f'outside {name}')))
            for name, part, _ in partitions:
                groupings.append((name, label_sites(d['chrom'].to_numpy(), d['strand'].to_numpy(),
                                                    d['start'].to_numpy(), part)))
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
    if a.shared_strata:
        pd.DataFrame(moves).to_csv(f'{a.out_dir}/stratum_moves.tsv', sep='\t', index=False, float_format='%.4g')
        pd.DataFrame(cross).to_csv(f'{a.out_dir}/stratum_crosstab.tsv', sep='\t', index=False)
    if new_calls or a.partition or a.truth_set or a.shared_strata:
        cfg = [('mod', a.mod), ('truth_set', 'from --truth-tsv' if a.truth_tsv else st),
               ('strata_all_dir', a.all_dir), ('strata_q40_dir', a.q40_dir), ('calls_all_dir', calls_all),
               ('calls_q40_dir', calls_q40), ('fail_as_canonical', a.fail_as_canonical),
               ('calls_threshold', '' if a.calls_threshold is None else a.calls_threshold),
               ('min_valid', a.min_valid), ('pct', a.pct), ('mask_bed', a.mask_bed or ''),
               ('contig_regex', a.contig_regex), ('shared_strata', a.shared_strata)] + [(f'strata:{n}', p) for n, p in cr.pairs(a.strata)] \
            + [(f'partition:{n}', p) for n, p in cr.pairs(a.partition)] + [(f'model:{n}', r) for n, r in models]
        pd.DataFrame(cfg, columns=['setting', 'value']).to_csv(f'{a.out_dir}/settings.tsv', sep='\t', index=False)

    (m1, _), (m2, _) = models[0], models[1]
    k = ['sample', 'reads', 'grouping', 'stratum']
    g = S[S['model'] == m1].merge(S[S['model'] == m2], on=k, suffixes=(f'_{m1}', f'_{m2}'))
    g['precision_gap'] = g[f'precision_{m1}'] - g[f'precision_{m2}']
    g['recall_gap'] = g[f'recall_{m1}'] - g[f'recall_{m2}']
    G = g[k + ['precision_gap', 'recall_gap', f'called_{m1}', f'called_{m2}', f'truth_per_1000_covered_{m1}']]
    G.to_csv(f'{a.out_dir}/gaps.tsv', sep='\t', index=False, float_format='%.4g')

    groups = [x for x in S['grouping'].unique() if x != 'all sites']
    orders = {'unique-read fraction': U_LABELS + [NO_NF], 'unique-read fraction (shared)': U_LABELS + [SHARED_NONE]}
    for name, _, order in partitions:
        if order:
            orders[name] = order + ['unlabelled']
    calls_desc = ''
    if new_calls:
        calls_desc = (f', calls at fixed {a.calls_threshold:g}' if a.calls_threshold is not None else ', calls') \
            + (' (fails canonical)' if a.fail_as_canonical else '')
    fig, axes = plt.subplots(len(groups), 2, figsize=(8.5, 2.6 * len(groups)), layout='constrained', squeeze=False)
    cols = {m1: '#d95f02', m2: '#1b9e77'}
    for i, gname in enumerate(groups):
        sub = S[(S['grouping'] == gname)]
        order = [s for s in (orders.get(gname) or sorted(sub['stratum'].unique())) if s in set(sub['stratum'])]
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
        ax.set_ylabel('precision'); ax.set_title(f'{a.mod} against {st}{calls_desc}: by {gname}', fontsize=8)
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
