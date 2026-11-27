#!/usr/bin/env python3
"""mock_mappability.py: tests of scripts/cmp_mappability.py (--partition, --truth-set, --calls-*-dir,
--fail-as-canonical, --calls-threshold) and scripts/make_transcript_partition.py on synthetic data.

    python tests/mock_mappability.py --work-dir mock/mapp [--old-script old/cmp_mappability.py]

1. default outputs byte-identical: with --old-script (the version before the new options, with its cmp_rules.py and
   cmp_ortho.py alongside), both versions on the same mock with --strata; strata.tsv, gaps.tsv, mappability.png and
   stderr must be identical, and mappability.pdf apart from its /CreationDate (if two old runs agree on the rest).
2. --partition on the no-filter calls: every strata.tsv row of the partition grouping recomputed independently
   (plain Python over the bedMethyl text, a per-position label table).
3. definition 2: --calls-all-dir/--calls-q40-dir on fixed-0.70 pileups with --fail-as-canonical --calls-threshold 0.7;
   every row of every grouping recomputed independently (depth N_valid + N_fail, occupancy from N_mod, unique-read
   fraction from the no-filter pair, 'no no-filter coverage' where the no-filter all-reads pileup has no read).
4. refusals: overlapping partition intervals, a wrong --calls-threshold, a calls pileup posing as the no-filter pair.
5. make_transcript_partition.py on a hand-made GTF (+ and - strand coding transcripts, a stop codon split across
   exons, a stop window crossing an exon junction, non-coding and untagged transcripts, a contig without genes)
   against labels written out by hand.
Exit status 1 if any check fails.
"""
import argparse
import filecmp
import gzip
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), 'scripts')
SAMPLES = ['S1', 'S2', 'S3']
RUNS = {'v6': 'runA', 'v5.2': 'runC'}
CONTIGS = {'chr1': 6000, 'chr2': 4000}
U_BINS = [(-0.01, 0.5, 'mostly multimapped (<50% unique)'), (0.5, 0.95, 'partly multimapped (50-95%)'),
          (0.95, 1.01, 'unique (>=95%)')]
NO_NF = 'no no-filter coverage'
FAILS = []


def ok(cond, msg):
    print(('PASS ' if cond else 'FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


def bedline(chrom, start, code, strand, nvalid, nmod, nfail=0):
    pct = 100 * nmod / nvalid if nvalid else 0
    return '\t'.join(map(str, [chrom, start, start + 1, code, nvalid, strand, start, start + 1, '255,0,0', nvalid,
                               f'{pct:.2f}', nmod, nvalid - nmod, 0, 0, nfail, 0, 0])) + '\n'


def summary(path, th, mapq):
    with open(path, 'w') as f:
        f.write('[2026-11-27 10:00:00] modkit 0.6.4; BAM x.bam; reference ref.fa; 16 threads\n')
        if mapq:
            f.write(f'[2026-11-27 10:00:00] keeping primary alignments with MAPQ >= {mapq}\n')
        f.write(f'[2026-11-27 10:00:01] thresholds (fixed, 10th percentile): A:{th} C:{th} G:{th} T:{th} in 0 min\n')


def make(w, seed=7):
    rng = np.random.default_rng(seed)
    os.makedirs(w, exist_ok=True)
    with open(f'{w}/ref.fa', 'w') as f, open(f'{w}/ref.fa.fai', 'w') as fai:
        off = 0
        for c, n in CONTIGS.items():
            s = ''.join(rng.choice(list('ACGT'), n))
            head, body = f'>{c}\n', ''.join(s[i:i + 60] + '\n' for i in range(0, n, 60))
            fai.write(f'{c}\t{n}\t{off + len(head)}\t60\t61\n'); f.write(head + body); off += len(head) + len(body)
    sites = [(c, int(p), s) for c, n in CONTIGS.items() for s in '+-'
             for p in rng.choice(n - 10, n // 4, replace=False) + 5]
    true = rng.random(len(sites)) < 0.06
    level = np.where(true, rng.uniform(0.25, 0.8, len(sites)), rng.uniform(0.0, 0.12, len(sites)))
    level[~true & (rng.random(len(sites)) < 0.03)] = 0.3                      # false positives
    depth = rng.lognormal(math.log(40), 0.8, len(sites))
    u = rng.choice([0.2, 0.7, 1.0], len(sites), p=[0.15, 0.15, 0.7])
    with open(f'{w}/truth.tsv', 'w') as f:
        f.write('chrom\tpos\tstrand\n')
        for (c, p, s), t, r in zip(sites, true, rng.random(len(sites))):
            if (t and r < 0.8) or (not t and r < 0.004):                      # incomplete, with a few misses
                f.write(f'{c}\t{p + 1}\t{s}\n')
    for label, run in RUNS.items():
        noise = 0.03 if label == 'v6' else 0.0
        for k, smp in enumerate(SAMPLES):
            dirs = {x: f'{w}/{x}/{run}' for x in ['nf_all', 'nf_q40', 'f070_all', 'f070_q40']}
            for x in dirs.values():
                os.makedirs(x, exist_ok=True)
            fh = {x: open(f'{d}/{smp}.bed', 'w') for x, d in dirs.items()}
            for i, (c, p, s) in enumerate(sites):
                n = rng.poisson(depth[i] * (0.8 + 0.2 * k))
                if n == 0:
                    continue
                pm = min(level[i] + noise * rng.random(), 1)
                m = rng.binomial(n, pm)
                nq = rng.binomial(n, u[i]) if rng.random() > 0.02 else 0
                mq = rng.binomial(m, nq / n) if n else 0
                fh['nf_all'].write(bedline(c, p, 'a', s, n, m))
                if nq:
                    fh['nf_q40'].write(bedline(c, p, 'a', s, nq, mq))
                fh['nf_all'].write(bedline(c, p, 'm', s, n, n))                # another code: must be ignored
                for x, (nn, mm) in {'f070_all': (n, m), 'f070_q40': (nq, mq)}.items():
                    if not nn:
                        continue
                    f = rng.binomial(nn, 0.2)                                   # calls below 0.70 fail
                    fm = rng.binomial(min(mm, f), 0.5)                          # some of them were modified calls
                    other = rng.poisson(0.3)                                    # failed calls on other-base reads
                    fh[x].write(bedline(c, p, 'a', s, nn - f, mm - fm, f + other))
            for j in range(10):                                                 # calls pileup only: no no-filter read
                c, p = 'chr2', 3995 + (j % 5)                                   # (sites end at 3994)
                s = '+' if j < 5 else '-'
                fh['f070_all'].write(bedline(c, p, 'a', s, 30, 12 if j % 2 else 2, 5))
            for x in fh.values():
                x.close()
            summary(f'{dirs["nf_all"]}/{smp}.summary.txt', 0, None)
            summary(f'{dirs["nf_q40"]}/{smp}.summary.txt', 0, 40)
            summary(f'{dirs["f070_all"]}/{smp}.summary.txt', 0.7, None)
            summary(f'{dirs["f070_q40"]}/{smp}.summary.txt', 0.7, 40)
    with open(f'{w}/lowmap.bed', 'w') as f:
        f.write('chr1\t1000\t1800\nchr1\t4000\t4100\nchr2\t0\t700\n')
    with gzip.open(f'{w}/regions.bed.gz', 'wt') as f:                          # stranded, chr2 partly unlabelled
        f.write('#mock partition\n#labels\tCDS,intron,intergenic\n')
        for c, s, a, b, lab in [('chr1', '+', 0, 2000, 'intergenic'), ('chr1', '+', 2000, 3500, 'CDS'),
                                ('chr1', '+', 3500, 6000, 'intron'), ('chr1', '-', 0, 4000, 'intron'),
                                ('chr1', '-', 4000, 6000, 'CDS'), ('chr2', '.', 0, 2500, 'intergenic')]:
            f.write(f'{c}\t{a}\t{b}\t{lab}\t0\t{s}\n')
    with open(f'{w}/regions_overlap.bed', 'w') as f:
        f.write('chr1\t0\t2000\tA\t0\t+\nchr1\t1500\t3000\tB\t0\t+\n')
    return w


def read_bed(path, fail_canonical=False):
    out = {}
    for line in open(path):
        x = line.rstrip('\n').split('\t')
        if x[3] != 'a':
            continue
        nvalid, nmod, nfail = int(x[9]), int(x[11]), int(x[15])
        if fail_canonical:
            n = nvalid + nfail
            out[(x[0], int(x[1]), x[5])] = (n, 100 * nmod / max(n, 1))
        else:
            out[(x[0], int(x[1]), x[5])] = (nvalid, float(x[10]))
    return out


def partition_table(path):
    lab = {}
    for line in gzip.open(path, 'rt') if path.endswith('.gz') else open(path):
        if line.startswith('#'):
            continue
        c, a, b, l, _, s = line.rstrip('\n').split('\t')[:6]
        for st in (['+', '-'] if s == '.' else [s]):
            for p in range(int(a), int(b)):
                lab[(c, st, p)] = l
    return lab


def expected_rows(w, calls_all, calls_q40, fail_canonical, new_calls, lowmap, part, min_valid=20, pct=20):
    truth = {(l.split('\t')[0], int(l.split('\t')[1]) - 1, l.rstrip('\n').split('\t')[2])
             for l in list(open(f'{w}/truth.tsv'))[1:]}
    iv = [tuple(l.split('\t')[:3]) for l in open(lowmap)]
    rows = {}
    for smp in SAMPLES:
        for label, run in RUNS.items():
            A = read_bed(f'{w}/{calls_all}/{run}/{smp}.bed', fail_canonical)
            Q = read_bed(f'{w}/{calls_q40}/{run}/{smp}.bed', fail_canonical)
            NA, NQ = read_bed(f'{w}/nf_all/{run}/{smp}.bed'), read_bed(f'{w}/nf_q40/{run}/{smp}.bed')
            sites = []
            for k, (n, p) in A.items():
                if n < min_valid:
                    continue
                nq, pq = Q.get(k, (0, float('nan')))
                if new_calls:
                    na = NA.get(k, (0, 0))[0]
                    uf = min(NQ.get(k, (0, 0))[0] / na, 1) if na > 0 else float('nan')
                else:
                    uf = min(nq / n, 1)
                us = NO_NF if uf != uf else next(l for lo, hi, l in U_BINS if lo < uf <= hi)
                inlow = any(k[0] == c and int(a) <= k[1] < int(b) for c, a, b in iv)
                g = {'unique-read fraction': us, 'low': 'in low' if inlow else 'outside low',
                     'region': part.get((k[0], k[2], k[1]), 'unlabelled'), 'all sites': 'all sites'}
                sites.append(dict(k=k, n=n, p=p, nq=nq, pq=pq, t=(k[0], k[1], k[2]) in truth, g=g))
            for reads in ['all reads', 'MAPQ >= 40']:
                for gname in ['unique-read fraction', 'low', 'region', 'all sites']:
                    for s in {x['g'][gname] for x in sites}:
                        sel = [x for x in sites if x['g'][gname] == s and (reads == 'all reads' or x['nq'] >= min_valid)]
                        called = [x for x in sel if (x['p'] if reads == 'all reads' else x['pq']) >= pct]
                        ts = [x for x in sel if x['t']]
                        sh = [x['p'] - x['pq'] for x in called if x['pq'] == x['pq']]
                        rows[(smp, label, reads, gname, s)] = dict(
                            covered=len(sel), truth_covered=len(ts), called=len(called),
                            precision=sum(x['t'] for x in called) / len(called) if called else float('nan'),
                            recall=sum(1 for x in ts if (x['p'] if reads == 'all reads' else x['pq']) >= pct) / len(ts)
                            if ts else float('nan'),
                            median_pct_shift_at_calls=float(np.median(sh)) if sh else float('nan'))
    return rows


def compare(tsv, rows, what):
    import csv
    got, bad, n = {}, 0, 0
    for r in csv.DictReader(open(tsv), delimiter='\t'):
        got[(r['sample'], r['model'], r['reads'], r['grouping'], r['stratum'])] = r
    for k, e in rows.items():
        if k[3] not in {g for _, _, _, g, _ in got} and k[3] in ('low', 'region'):
            continue
        r = got.get(k)
        if r is None:
            if e['covered'] > 0:
                bad += 1
                print(f'  missing row {k}')
            continue
        n += 1
        for f, v in e.items():
            x = float(r[f]) if r[f] != '' else float('nan')
            if not ((x != x and v != v) or (x == x and v == v and abs(x - v) <= 1e-3 * max(1, abs(v)))):
                bad += 1
                if bad < 6:
                    print(f'  {k} {f}: script {r[f]}, recomputed {v}')
    extra = [k for k in got if k not in rows and k[3] in {g for _, _, _, g, _ in rows}]
    ok(bad == 0 and n > 0 and not extra, f'{what}: {n} rows recomputed, {bad} mismatches, {len(extra)} unexpected rows')


def run(args, out=None):
    p = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    if out:
        open(out, 'w').write(p.stderr)
    return p


def base_args(script, w, out):
    return [script, '--mod', 'm6A', '--model', 'v6=runA', '--model', 'v5.2=runC'] + \
        sum([['--sample', s] for s in SAMPLES], []) + \
        ['--all-dir', f'{w}/nf_all', '--q40-dir', f'{w}/nf_q40', '--fasta', f'{w}/ref.fa',
         '--truth-tsv', f'{w}/truth.tsv', '--out-dir', out]


def build_test(w):
    g = [('T1', '+', True, [(101, 200), (301, 400), (501, 700)], [(151, 200), (301, 400), (501, 547)], [(548, 550)]),
         ('T3', '+', True, [(221, 260), (351, 380)], [], []),
         ('T2', '-', True, [(1001, 1100), (1201, 1300)], [(1054, 1100), (1201, 1250)], [(1051, 1053)]),
         ('T4', '+', False, [(2001, 2100)], [], []),
         ('T5', '-', True, [(101, 200)], [], []),
         ('T6', '+', True, [(3001, 3010), (3101, 3200)], [(3001, 3010), (3101, 3102)], [(3103, 3105)]),
         ('T7', '+', True, [(4001, 4050), (4101, 4150)], [(4001, 4049)], [(4050, 4050), (4101, 4102)])]
    with open(f'{w}/test.gtf', 'w') as f:
        f.write('##description: mock\n')
        for t, s, basic, ex, cds, stop in g:
            tag = ' tag "basic";' if basic else ' tag "TAGENE";'
            at = f'gene_id "G{t}"; transcript_id "{t}"; gene_type "x";{tag}'
            f.write(f'chr1\tT\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{s}\t.\t{at}\n')
            for kind, xs in [('exon', ex), ('CDS', cds), ('stop_codon', stop)]:
                for a, b in xs:
                    f.write(f'chr1\tT\t{kind}\t{a}\t{b}\t.\t{s}\t.\t{at}\n')
    with open(f'{w}/test.fa.fai', 'w') as f:
        f.write('chr1\t5000\t6\t60\t61\nchr2\t1000\t5200\t60\t61\n')
    p = run([f'{SCRIPTS}/make_transcript_partition.py', '--gtf', f'{w}/test.gtf', '--fai', f'{w}/test.fa.fai',
             '--window', '20', '--out', f'{w}/test_regions.bed.gz'])
    if p.returncode:
        ok(False, f'make_transcript_partition.py ran: {p.stderr[-500:]}')
        return
    S = 'stop codon +/-20 nt'
    plus = [(101, 150, "5' UTR"), (151, 200, 'CDS'), (201, 220, 'intron'), (221, 260, 'non-coding exon'),
            (261, 300, 'intron'), (301, 400, 'CDS'), (401, 500, 'intron'), (501, 527, 'CDS'), (528, 570, S), (571, 700, "3' UTR rest"),
            (3001, 3010, S), (3011, 3100, 'intron'), (3101, 3125, S), (3126, 3200, "3' UTR rest"),
            (4001, 4029, 'CDS'), (4030, 4050, S), (4051, 4100, 'intron'), (4101, 4122, S),
            (4123, 4150, "3' UTR rest")]
    minus = [(101, 200, 'non-coding exon'), (1001, 1030, "3' UTR rest"), (1031, 1073, S), (1074, 1100, 'CDS'),
             (1101, 1200, 'intron'), (1201, 1250, 'CDS'), (1251, 1300, "5' UTR")]
    exp = {}
    for c, n in [('chr1', 5000), ('chr2', 1000)]:
        for s in '+-':
            for p0 in range(n):
                exp[(c, s, p0)] = 'intergenic'
    for s, xs in [('+', plus), ('-', minus)]:
        for a, b, l in xs:
            for p1 in range(a, b + 1):
                exp[('chr1', s, p1 - 1)] = l
    got = partition_table(f'{w}/test_regions.bed.gz')
    diff = [k for k in exp if got.get(k) != exp[k]]
    ok(not diff and len(got) == len(exp),
       f'make_transcript_partition.py: {len(exp):,} positions x strands, {len(diff)} labelled differently'
       + (f', e.g. {diff[0]}: got {got.get(diff[0])}, expected {exp[diff[0]]}' if diff else ''))
    header = [l for l in gzip.open(f'{w}/test_regions.bed.gz', 'rt') if l.startswith('#labels')]
    ok(header == [f"#labels\t5' UTR,CDS,{S},3' UTR rest,non-coding exon,intron,intergenic\n"], 'partition #labels header')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--work-dir', required=True)
    ap.add_argument('--old-script', help='cmp_mappability.py before the new options (cmp_rules/cmp_ortho alongside)')
    a = ap.parse_args()
    w = make(a.work_dir)
    new = f'{SCRIPTS}/cmp_mappability.py'
    strata = ['--strata', f'low={w}/lowmap.bed']

    # 1. default outputs
    if a.old_script:
        po = run(base_args(a.old_script, w, f'{w}/old') + strata, f'{w}/old.stderr')
        po2 = run(base_args(a.old_script, w, f'{w}/old2') + strata, f'{w}/old2.stderr')
        pn = run(base_args(new, w, f'{w}/new') + strata, f'{w}/new.stderr')
        ok(po.returncode == 0 and pn.returncode == 0, f'old and new run (exit {po.returncode}, {pn.returncode}) '
           + (po.stderr[-300:] + pn.stderr[-300:] if po.returncode or pn.returncode else ''))
        for f in ['strata.tsv', 'gaps.tsv', 'mappability.png']:
            ok(filecmp.cmp(f'{w}/old/{f}', f'{w}/new/{f}', shallow=False), f'default {f} byte-identical to the old version')
        ok(filecmp.cmp(f'{w}/old.stderr', f'{w}/new.stderr', shallow=False), 'default stderr identical')
        pdf = {k: re.sub(rb'/CreationDate \(D:[^)]*\)', b'', open(f'{w}/{k}/mappability.pdf', 'rb').read())
               for k in ['old', 'old2', 'new']}
        if pdf['old'] == pdf['old2']:
            ok(pdf['old'] == pdf['new'], 'default mappability.pdf identical apart from its creation date')
        else:
            print('NOTE mappability.pdf differs between two runs of the old version beyond its date: not compared')
        ok(not os.path.exists(f'{w}/new/settings.tsv'), 'no settings.tsv without the new options')
    part = partition_table(f'{w}/regions.bed.gz')

    # 2. partition, no-filter calls
    p = run(base_args(new, w, f'{w}/part') + strata + ['--partition', f'region={w}/regions.bed.gz'])
    ok(p.returncode == 0, 'run with --partition' + (': ' + p.stderr[-400:] if p.returncode else ''))
    if p.returncode == 0:
        rows = expected_rows(w, 'nf_all', 'nf_q40', False, False, f'{w}/lowmap.bed', part)
        compare(f'{w}/part/strata.tsv', rows, '--partition (no filter)')
        if a.old_script:
            import csv
            old = [r for r in csv.DictReader(open(f'{w}/old/strata.tsv'), delimiter='\t')]
            newr = [r for r in csv.DictReader(open(f'{w}/part/strata.tsv'), delimiter='\t') if r['grouping'] != 'region']
            ok(old == newr, 'adding --partition leaves every other row unchanged')

    # 3. definition 2
    d2 = ['--calls-all-dir', f'{w}/f070_all', '--calls-q40-dir', f'{w}/f070_q40', '--fail-as-canonical',
          '--calls-threshold', '0.7']
    p = run(base_args(new, w, f'{w}/def2') + strata + ['--partition', f'region={w}/regions.bed.gz'] + d2)
    ok(p.returncode == 0, 'definition-2 run' + (': ' + p.stderr[-400:] if p.returncode else ''))
    if p.returncode == 0:
        rows = expected_rows(w, 'f070_all', 'f070_q40', True, True, f'{w}/lowmap.bed', part)
        compare(f'{w}/def2/strata.tsv', rows, 'definition 2 (fixed 0.70, fails canonical), all groupings')
        nn = sum(1 for k, v in rows.items() if k[4] == NO_NF and v['covered'])
        ok(nn > 0, f"definition 2: 'no no-filter coverage' rows present ({nn})")
        ok(os.path.exists(f'{w}/def2/settings.tsv'), 'settings.tsv written')

    # 4. refusals
    p = run(base_args(new, w, f'{w}/bad1') + ['--partition', f'r={w}/regions_overlap.bed'])
    ok(p.returncode != 0 and 'overlapping' in p.stderr, 'overlapping partition refused')
    p = run(base_args(new, w, f'{w}/bad2') + d2[:-1] + ['0.8'])
    ok(p.returncode != 0 and 'fixed threshold 0.8' in p.stderr, 'wrong --calls-threshold refused')
    bad = base_args(new, w, f'{w}/bad3')
    bad[bad.index('--all-dir') + 1] = f'{w}/f070_all'
    p = run(bad + d2)
    ok(p.returncode != 0 and 'no-filter all-reads' in p.stderr, 'a fixed-0.70 pileup as the strata pair refused')
    p = run(base_args(new, w, f'{w}/bad4') + ['--calls-all-dir', f'{w}/f070_all'])
    ok(p.returncode != 0, '--calls-all-dir without --calls-q40-dir refused')

    # 5. partition builder
    build_test(w)
    print(f'{len(FAILS)} failed' if FAILS else 'all checks passed')
    sys.exit(1 if FAILS else 0)


if __name__ == '__main__':
    main()
