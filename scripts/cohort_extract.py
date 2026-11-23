#!/usr/bin/env python3
"""Per-sample read-level QC for cohort_plots.py figure 11 (and its covariates in figure 10), from the
Stage 1 BAMs and the IsoQuant requantification of RNA_analysis.nf.

    cohort_extract.py all --results <results> --samplesheet <analysis_samplesheet.csv> [--ref GRCh38]
                          [--fraction 0.02] [--max-reads 100000] [--threads 8]
    cohort_extract.py bam <sample> <bam> --out <dir> [--fraction 0.02] [--sam-text FILE]
    cohort_extract.py isoquant <sample> <read_assignments.tsv.gz> <catalogue.gtf> --out <dir>

bam       samtools view -F 2308 -s <seed>.<fraction> of the BAM (primary alignments; samtools must be on
          PATH, e.g. the pipeline's mods image) -> <sample>.read_metrics.tsv.gz: read length, MAPQ, Dorado
          mean qscore (qs) and translocation speed from the move table (mv:B:c,stride,moves...): bases /
          (move entries x stride / sample rate; 4 kHz for RNA004). Dorado's ns/ts tags are not carried by
          Stage 1 (ts clashes with minimap2's), the move table is. Also <sample>.models.tsv from the @RG DS.
isoquant  IsoQuant read assignments (unique / unique_minor_difference reads) against the catalogue GTF ->
          <sample>.integrity.tsv: per transcript-length bin, reads, the fraction whose 5' end lies within
          --five-prime nt of the isoform's 5' end, and the median fraction of the isoform covered; and
          <sample>.read_categories.tsv: IsoQuant's Classification= of each read (or its assignment type).
all       both for every sample of the sheet, into <results>/sample_qc (the default --sample-qc of
          cohort_plots.py). BAM paths are relative to the sheet's directory, as Stage 1 writes them; read
          assignments are isoforms/<ref>/requant/<s>/<s>.read_assignments.tsv.gz and the catalogue
          isoforms/<ref>/cohort/cohort.<ref>.catalogue.gtf.

Run it on a compute node (it streams every BAM once): e.g.
    sbatch -p short -c 8 --mem 16G -t 4:00:00 --wrap "cohort_extract.py all --results R --samplesheet S --threads 8"
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

LEN_EDGES = [0, 500, 1000, 1500, 2000, 3000, 5000, 10 ** 9]
UNIQUE = {'unique', 'unique_minor_difference'}
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


# ---------------------------------------------------------------- BAM
def bam_metrics(sample, bam, out, fraction=0.02, max_reads=100000, sample_rate=4000.0, seed=42, sam_text=None):
    if sam_text:
        src = opener(sam_text)
        header = [ln for ln in opener(sam_text) if ln.startswith('@')]
        proc = None
    else:
        if not shutil_which('samtools'):
            sys.exit('cohort_extract.py: samtools is not on PATH (use the pipeline mods image)')
        header = subprocess.run(['samtools', 'view', '-H', bam], capture_output=True, text=True, check=True).stdout.splitlines(True)
        frac = f'{seed}.{str(fraction).split(".")[1]}' if fraction < 1 else str(seed)
        cmd = ['samtools', 'view', '-F', '2308'] + (['-s', frac] if fraction < 1 else []) + [bam]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        src = proc.stdout
    rows = []
    for line in src:
        if line.startswith('@'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 2308:
            continue
        rl = len(f[9]) if f[9] != '*' else 0
        q, speed = np.nan, np.nan
        for tag in f[11:]:
            if tag.startswith('qs:'):
                try:
                    q = float(tag.split(':', 2)[2])
                except ValueError:
                    pass
            elif tag.startswith('mv:B:c,'):
                body = tag[7:]
                k = body.find(',')
                if k > 0:
                    stride = int(body[:k])
                    rest = body[k + 1:]
                    entries = rest.count(',') + 1
                    moves = rest.count('1')
                    if entries and moves:
                        speed = moves / (entries * stride / sample_rate)
        rows.append((rl, int(f[4]), q, speed))
        if len(rows) >= max_reads:
            break
    if proc is not None:
        proc.stdout.close()
        proc.terminate()
        proc.wait()
    elif sam_text:
        src.close()
    os.makedirs(out, exist_ok=True)
    pd.DataFrame(rows, columns=['read_len', 'mapq', 'qscore', 'speed_bps']).to_csv(
        os.path.join(out, f'{sample}.read_metrics.tsv.gz'), sep='\t', index=False, float_format='%.4g', compression='gzip')
    models = []
    for h in header:
        if h.startswith('@RG'):
            d = dict(x.split(':', 1) for x in h.rstrip('\n').split('\t')[1:] if ':' in x)
            ds = d.get('DS', '')
            bm = re.search(r'basecall_model=(\S+)', ds)
            mm = re.search(r'modbase_models=(\S+)', ds)
            models.append((d.get('ID', '.'), bm.group(1) if bm else '.', mm.group(1) if mm else '.'))
    pd.DataFrame(models, columns=['read_group', 'basecall_model', 'modbase_models']).to_csv(
        os.path.join(out, f'{sample}.models.tsv'), sep='\t', index=False)
    return sample, len(rows)


def shutil_which(x):
    import shutil
    return shutil.which(x)


# ---------------------------------------------------------------- IsoQuant
def load_models(gtf):
    """transcript_id -> (strand, sorted exons [(start0, end)], cumulative offsets, length)."""
    rx = re.compile(r'transcript_id "([^"]+)"')
    ex = {}
    with opener(gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 8)
            if len(f) < 9 or f[2] != 'exon':
                continue
            m = rx.search(f[8])
            if m:
                ex.setdefault(m.group(1), [f[6], []])[1].append((int(f[3]) - 1, int(f[4])))
    out = {}
    for t, (st, e) in ex.items():
        e.sort()
        lens = np.array([b - a for a, b in e])
        out[t] = (st, e, np.concatenate([[0], np.cumsum(lens)[:-1]]), int(lens.sum()))
    return out


def tx_coord(model, g):
    """Transcript coordinate (0 = 5' end) of genomic position g, snapped to the nearest exon base."""
    st, ex, cum, L = model
    best, bd = None, None
    for (a, b), c0 in zip(ex, cum):
        if a <= g < b:
            best = c0 + (g - a)
            break
        d = a - g if g < a else g - (b - 1)
        if bd is None or d < bd:
            bd, best = d, (c0 if g < a else c0 + (b - a) - 1)
    return best if st == '+' else L - 1 - best


def isoquant_metrics(sample, assignments, gtf, out, five_prime=50, fraction=1.0, seed=42, models=None):
    models = models or load_models(gtf)
    names, skip = None, 0
    with opener(assignments) as fh:
        for n_, line in enumerate(fh):
            if not line.startswith('#'):
                break
            cand = line.lstrip('#').rstrip('\n').split('\t')
            if 'read_id' in cand and 'isoform_id' in cand:
                names = cand
            skip = n_ + 1
    if names is None:
        sys.exit(f'cohort_extract.py: {assignments}: no "#read_id ..." header line')
    use = [c for c in ('read_id', 'strand', 'isoform_id', 'assignment_type', 'exons', 'additional_info') if c in names]
    rng = np.random.default_rng(seed)
    bins = {i: [0, 0, []] for i in range(len(LEN_EDGES) - 1)}
    cats, seen = {}, set()
    rx_c = re.compile(r'Classification=([^;]+)')
    for chunk in pd.read_csv(assignments, sep='\t', header=None, names=names, usecols=use, skiprows=skip, dtype=str,
                             keep_default_na=False, chunksize=500000):
        if fraction < 1:
            chunk = chunk[rng.random(len(chunk)) < fraction]
        first = ~chunk['read_id'].duplicated() & ~chunk['read_id'].isin(seen)
        for rid, at, info in zip(chunk['read_id'][first], chunk['assignment_type'][first],
                                 chunk['additional_info'][first] if 'additional_info' in chunk else [''] * first.sum()):
            m = rx_c.search(info or '')
            k = m.group(1) if m else at
            cats[k] = cats.get(k, 0) + 1
        seen.update(chunk['read_id'][first].tolist())
        u = chunk[chunk['assignment_type'].isin(UNIQUE)]
        for t, sd, exs in zip(u['isoform_id'], u['strand'], u['exons']):
            mdl = models.get(t)
            if mdl is None or not exs or exs == '.':
                continue
            blocks = [tuple(int(x) for x in b.split('-')) for b in exs.split(',') if '-' in b]
            if not blocks:
                continue
            lo, hi = min(b[0] for b in blocks) - 1, max(b[1] for b in blocks) - 1
            g5, g3 = (lo, hi) if mdl[0] == '+' else (hi, lo)
            p5, p3 = tx_coord(mdl, g5), tx_coord(mdl, g3)
            L = mdl[3]
            j = int(np.searchsorted(LEN_EDGES, L, side='right') - 1)
            bins[j][0] += 1
            bins[j][1] += p5 <= five_prime
            if len(bins[j][2]) < 20000:
                bins[j][2].append(min(1.0, (abs(p3 - p5) + 1) / L))
    os.makedirs(out, exist_ok=True)
    rows = [(LEN_EDGES[i], LEN_EDGES[i + 1] if LEN_EDGES[i + 1] < 10 ** 9 else 50000, v[0],
             v[1] / v[0] if v[0] else float('nan'), float(np.median(v[2])) if v[2] else float('nan')) for i, v in bins.items()]
    pd.DataFrame(rows, columns=['length_lo', 'length_hi', 'reads', 'reached_5p', 'median_covered_fraction']).to_csv(
        os.path.join(out, f'{sample}.integrity.tsv'), sep='\t', index=False, float_format='%.4f')
    pd.DataFrame(sorted(cats.items(), key=lambda kv: -kv[1]), columns=['category', 'reads']).to_csv(
        os.path.join(out, f'{sample}.read_categories.tsv'), sep='\t', index=False)
    return sample, sum(v[0] for v in bins.values())


def _all_job(job):
    kind, args = job
    return (bam_metrics if kind == 'bam' else isoquant_metrics)(*args[0], **args[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('bam')
    b.add_argument('sample')
    b.add_argument('bam')
    b.add_argument('--out', required=True)
    b.add_argument('--sam-text', help='read SAM text from this file instead of running samtools (tests)')
    i = sub.add_parser('isoquant')
    i.add_argument('sample')
    i.add_argument('assignments')
    i.add_argument('gtf')
    i.add_argument('--out', required=True)
    a = sub.add_parser('all')
    a.add_argument('--results', required=True)
    a.add_argument('--samplesheet', required=True)
    a.add_argument('--ref', default='GRCh38')
    a.add_argument('--out', help='default: <results>/sample_qc')
    a.add_argument('--threads', type=int, default=4)
    for p in (b, a):
        p.add_argument('--fraction', type=float, default=0.02, help='BAM fraction sampled (samtools view -s)')
        p.add_argument('--max-reads', type=int, default=100000)
        p.add_argument('--sample-rate', type=float, default=4000.0, help='signal sample rate in Hz (RNA004: 4000)')
    for p in (i, a):
        p.add_argument('--five-prime', type=int, default=50, help="nt from the isoform's 5' end that still counts as complete")
        p.add_argument('--assign-fraction', type=float, default=1.0, help='fraction of assignment rows read')
    args = ap.parse_args()
    if args.cmd == 'bam':
        s, n = bam_metrics(args.sample, args.bam, args.out, args.fraction, args.max_reads, args.sample_rate, sam_text=args.sam_text)
        print(f'cohort_extract.py bam: {s}: {n:,} reads', file=sys.stderr)
    elif args.cmd == 'isoquant':
        s, n = isoquant_metrics(args.sample, args.assignments, args.gtf, args.out, args.five_prime, args.assign_fraction)
        print(f'cohort_extract.py isoquant: {s}: {n:,} uniquely assigned reads', file=sys.stderr)
    else:
        sheet = pd.read_csv(args.samplesheet, dtype=str, keep_default_na=False)
        if 'reference' in sheet.columns:
            sheet = sheet[sheet['reference'] == args.ref]
        base = os.path.dirname(os.path.abspath(args.samplesheet))
        out = args.out or os.path.join(args.results, 'sample_qc')
        gtf = os.path.join(args.results, 'isoforms', args.ref, 'cohort', f'cohort.{args.ref}.catalogue.gtf')
        jobs = []
        for s, bam in zip(sheet['sample_id'], sheet['bam'] if 'bam' in sheet.columns else [''] * len(sheet)):
            bp = bam if bam.startswith('/') else os.path.join(base, bam)
            if bam and os.path.exists(bp):
                jobs.append(('bam', ((s, bp, out, args.fraction, args.max_reads, args.sample_rate), {})))
            else:
                print(f'cohort_extract.py: {s}: no BAM at {bp}', file=sys.stderr)
            ra = os.path.join(args.results, 'isoforms', args.ref, 'requant', s, f'{s}.read_assignments.tsv.gz')
            if os.path.exists(ra) and os.path.exists(gtf):
                jobs.append(('iso', ((s, ra, gtf, out, args.five_prime, args.assign_fraction), {})))
            else:
                print(f'cohort_extract.py: {s}: no {ra} or {gtf}', file=sys.stderr)
        with ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
            for (kind, _), (s, n) in zip(jobs, pool.map(_all_job, jobs)):
                print(f'cohort_extract.py {kind}: {s}: {n:,} reads', file=sys.stderr)
        print(f'cohort_extract.py: outputs in {out}')


if __name__ == '__main__':
    main()
