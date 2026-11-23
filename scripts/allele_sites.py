#!/usr/bin/env python3
"""Allele-specific modification at chosen sites. In each donor heterozygous for a SNV near a site, the reads
are split by the base they carry at the SNV and the modification calls at the site are compared between the
two alleles of the same sample; donors are then combined (Cochran-Mantel-Haenszel, donors as strata), and
every donor's genotype is set against its cohort-matrix occupancy (the dosage view).

    allele_sites.py --results R --samplesheet S --site chr6:32659873:-:m6A [--site chr17:76736425:-:Psi ...]
                    [--sites-file sites.tsv] --out DIR [--window 300] [--min-mapq 40]
                    [--container /path/mods.sif] [--reference FA --errortable gIVT.tsv.gz]

--site is chrom:position(1-based):strand:modification (m6A, m5C, Psi, Ino, Am, Cm, Um, Gm or a modkit code);
--sites-file holds the same, one per line. Per-read calls come from `samtools view -u -q MAPQ -F 2308 <bam>
<region> | modkit extract calls - - --include-bed <positions> --mapped-only $(cat <s>.modkit_flags.txt)`,
i.e. each sample's own call thresholds, exactly as the pipeline's feature step reads them. A read's allele is
its base at the SNV (read orientation, turned back to the genomic strand); reads with any other base, or a
deletion there, are left out. With --reference and --errortable, a SNV inside the site's 9-mer (+-4 nt) gives
each allele its own gIVT false-positive rate (the reference 9-mer with the read's base substituted), and both
raw and corrected occupancy are reported; outside the 9-mer both alleles share the reference rate.

Interpretation: a difference between alleles of the same donor is free of donor-level confounders, but a SNV
inside the 9-mer can change what the model reads as well as what the enzyme does: confirm those with a
synthetic oligo pair. A SNV AT the site (offset 0) changes the base itself: the table then shows what every
model calls on each allele (e.g. m5C/Cm on C reads where the reference has a U).

Needs samtools and modkit (on PATH, or through --container, e.g. the pipeline's mods.sif), pandas, numpy,
scipy, and cohort_plots.py next to this script (for the FASTA reader, the gIVT lookup and the CMH test).
"""
import argparse
import os
import re
import shlex
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_plots import Fasta, GivtTable, MOD_INDEX, cmh_2xk, plain, read_matrix   # noqa: E402

ALIAS = {'m6a': 'a', 'm5c': 'm', 'psi': '17802', 'pseu': '17802', 'ino': '17596', 'inosine': '17596',
         'am': '69426', 'cm': '19228', 'um': '19227', 'gm': '19229'}
MOD_BASE = {'a': 'A', '17596': 'A', '69426': 'A', 'm': 'C', '19228': 'C', '17802': 'T', '19227': 'T', '19229': 'G'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
HET = {'0/1', '1/0', '0|1', '1|0'}
CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b'ACGT'):
    CODE[_b] = _i
SHIFTS = (2 * np.arange(8, -1, -1)).astype(np.int64)


def parse_site(spec):
    f = spec.strip().split(':')
    if len(f) != 4:
        sys.exit(f'allele_sites.py: --site wants chrom:pos:strand:mod, got {spec!r}')
    code = f[3] if f[3] in MOD_BASE else ALIAS.get(f[3].lower())
    if code is None or f[2] not in '+-':
        sys.exit(f'allele_sites.py: bad site {spec!r}')
    return dict(chrom=f[0], pos0=int(f[1].replace(',', '')) - 1, strand=f[2], code=code, label=spec.strip())


def load_variants(path):
    if not os.path.exists(path):
        return None
    v = pd.read_csv(path, sep='\t', dtype=str, usecols=lambda c: c in {'chrom', 'pos', 'ref', 'alt', 'filter', 'gt', 'vaf', 'class'})
    ok = (v['ref'].str.len() == 1) & (v['alt'].str.len() == 1)
    if 'class' in v:
        ok &= v['class'] != 'editing'
    v = v[ok].copy()
    v['p0'] = v['pos'].astype(int) - 1
    return {c: g.sort_values('p0').reset_index(drop=True) for c, g in v.groupby('chrom')}


def near(var, chrom, pos0, w):
    if var is None or chrom not in var:
        return pd.DataFrame(columns=['p0', 'ref', 'alt', 'gt', 'filter', 'vaf', 'class'])
    g = var[chrom]
    a, b = np.searchsorted(g['p0'].to_numpy(), [pos0 - w, pos0 + w + 1])
    return g.iloc[a:b]


def run_extract(args, bam, flags, site, snvs, tmp):
    lo = min([site['pos0']] + [int(p) for p in snvs['p0']]) - 50
    hi = max([site['pos0']] + [int(p) for p in snvs['p0']]) + 50
    region = f"{site['chrom']}:{max(1, lo)}-{hi}"
    bed = os.path.join(tmp, 'positions.bed')
    with open(bed, 'w') as fh:
        for p in sorted({site['pos0'], *[int(x) for x in snvs['p0']]}):
            fh.write(f"{site['chrom']}\t{p}\t{p + 1}\tpos\t0\t{site['strand']}\n")
    sam = shlex.split(args.samtools)
    mk = shlex.split(args.modkit)
    cmd = (' '.join(map(shlex.quote, sam + ['view', '-u', '-q', str(args.min_mapq), '-F', '2308', bam, region])) + ' | '
           + ' '.join(map(shlex.quote, mk + ['extract', 'calls', '-', '-', '--include-bed', bed, '--mapped-only']))
           + f' {flags} --threads {args.threads} --suppress-progress')
    r = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f'allele_sites.py: extraction failed ({r.returncode}):\n  {cmd}\n{r.stderr[-2000:]}')
    if not r.stdout.strip():
        return pd.DataFrame()
    from io import StringIO
    t = pd.read_csv(StringIO(r.stdout), sep='\t', dtype=str)
    need = {'read_id', 'ref_position', 'call_code', 'canonical_base'}
    if not need <= set(t.columns):
        sys.exit('allele_sites.py: modkit extract output lacks ' + ', '.join(sorted(need - set(t.columns)))
                 + '; its columns are: ' + ', '.join(t.columns))
    t['ref_position'] = pd.to_numeric(t['ref_position'], errors='coerce')
    t = t[t['ref_position'] >= 0]
    if 'ref_strand' in t.columns:
        t = t[t['ref_strand'] == site['strand']]
    t['fail'] = t['fail'].astype(str).str.lower().isin(['true', '1']) if 'fail' in t.columns else False
    return t


def allele_fp(table, fasta, site, snv_p0, read_base_genomic):
    """gIVT rate for the site's 9-mer carrying this allele (genomic base at the SNV substituted)."""
    if table is None or not fasta.ok:
        return np.nan
    ctx = fasta.context(np.array([site['chrom']], dtype=object), np.array([site['pos0']]),
                        np.array([site['strand']], dtype=object), flank=4)[0].copy()
    d = (snv_p0 - site['pos0']) * (1 if site['strand'] == '+' else -1)
    if abs(d) <= 4:
        rna = read_base_genomic if site['strand'] == '+' else COMP[read_base_genomic]
        ctx[4 + d] = ord(rna)
    v = CODE[ctx]
    if (v == 255).any():
        return np.nan
    code = int((v.astype(np.int64) << SHIFTS).sum())
    i = MOD_INDEX[site['code']]
    return float(table.fps[i, code]) if table.present[i, code] else np.nan


def thresholds_of(path):
    d = {}
    if os.path.exists(path):
        for tok in open(path).read().split():
            if ',' in tok:
                k, v = tok.split(',', 1)
                try:
                    d[k] = float(v)
                except ValueError:
                    pass
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True)
    ap.add_argument('--samplesheet', required=True)
    ap.add_argument('--ref', default='GRCh38')
    ap.add_argument('--site', action='append', default=[])
    ap.add_argument('--sites-file')
    ap.add_argument('--out', required=True)
    ap.add_argument('--window', type=int, default=300, help='nt either side of the site searched for heterozygous SNVs')
    ap.add_argument('--min-mapq', type=int, default=40)
    ap.add_argument('--min-reads', type=int, default=5, help='reads per allele for a test')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--container', help='run samtools and modkit inside this image (apptainer exec -B /scratch,/projects)')
    ap.add_argument('--samtools', default='samtools')
    ap.add_argument('--modkit', default='modkit')
    ap.add_argument('--reference', default='/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
    ap.add_argument('--errortable', default='/projects/Genometechlab/referenceFastaFiles/gIVT/GM12878_genomic_IVT_refmatch_9mer_all_threshold.dorado_1.0.tsv.gz')
    a = ap.parse_args()
    if a.container:
        pre = f'apptainer exec -B /scratch,/projects {shlex.quote(a.container)}'
        a.samtools, a.modkit = f'{pre} samtools', f'{pre} modkit'
    specs = list(a.site)
    if a.sites_file:
        specs += [ln.split()[0] for ln in open(a.sites_file) if ln.strip() and not ln.startswith('#')]
    if not specs:
        sys.exit('allele_sites.py: give --site or --sites-file')
    sites = [parse_site(s) for s in specs]
    os.makedirs(a.out, exist_ok=True)
    sheet = pd.read_csv(a.samplesheet, dtype=str, keep_default_na=False)
    if 'reference' in sheet.columns:
        sheet = sheet[sheet['reference'] == a.ref]
    base = os.path.dirname(os.path.abspath(a.samplesheet))
    bams = {s: (b if b.startswith('/') else os.path.join(base, b)) for s, b in zip(sheet['sample_id'], sheet['bam'])}
    S = list(bams)
    R = a.results
    var = {s: load_variants(f'{R}/variants/{a.ref}/{s}/{s}.variants.tsv.gz') for s in S}
    keys = ['chrom', 'start', 'end', 'mod', 'strand']
    coh = os.path.join(R, 'cohort', a.ref, f'cohort.{a.ref}.mods.')
    pm = read_matrix(coh + 'pct_fpcorr.tsv.gz', keys) if os.path.exists(coh + 'pct_fpcorr.tsv.gz') else None
    nvm = read_matrix(coh + 'n_valid.tsv.gz', keys) if os.path.exists(coh + 'n_valid.tsv.gz') else None
    fasta = Fasta(a.reference)
    tables = {}
    rows, dosage = [], []
    for site in sites:
        print(f"\n== {site['chrom']}:{site['pos0'] + 1:,} ({site['strand']}) {plain(site['code'])}")
        if pm is not None:
            hit = pm[(pm['chrom'] == site['chrom']) & (pm['start'].astype(int) == site['pos0']) &
                     (pm['mod'] == site['code']) & (pm['strand'] == site['strand'])]
        else:
            hit = []
        for s in S:                                                 # dosage view: genotype against occupancy
            g = near(var[s], site['chrom'], site['pos0'], 4)
            occ = float(hit[s].iat[0]) if len(hit) and s in hit else np.nan
            nv = float(nvm.loc[hit.index[0], s]) if (len(hit) and nvm is not None and s in nvm) else np.nan
            gts = '; '.join(f"{(r.p0 - site['pos0']) * (1 if site['strand'] == '+' else -1):+d} {r.ref}>{r.alt} {r.gt}"
                            for r in g.itertuples()) or '-'
            dosage.append((site['label'], s, occ, nv, gts))
        for s in S:
            hets = near(var[s], site['chrom'], site['pos0'], a.window)
            hets = hets[hets['gt'].isin(HET) & (hets.get('filter', 'PASS') == 'PASS')] if len(hets) else hets
            if not len(hets):
                continue
            flags_path = f'{R}/mods/{a.ref}/{s}/{s}.modkit_flags.txt'
            flags = open(flags_path).read().strip() if os.path.exists(flags_path) else ''
            thr = thresholds_of(f'{R}/mods/{a.ref}/{s}/{s}.givt_flags.txt')
            key = tuple(sorted(thr.items()))
            if key not in tables and a.errortable and os.path.exists(a.errortable) and fasta.ok:
                tables[key] = GivtTable(a.errortable, thr)
            table = tables.get(key)
            with tempfile.TemporaryDirectory() as tmp:
                t = run_extract(a, bams[s], flags, site, hets, tmp)
            if not len(t):
                print(f'  {s}: no reads at MAPQ >= {a.min_mapq}')
                continue
            at_site = t[t['ref_position'] == site['pos0']].drop_duplicates('read_id').set_index('read_id')
            for snv in hets.itertuples():
                at_snv = t[t['ref_position'] == snv.p0].drop_duplicates('read_id').set_index('read_id')
                rb = at_snv['canonical_base'].str.upper().str.replace('U', 'T')
                gb = rb if site['strand'] == '+' else rb.map(COMP)        # read base back on the genomic strand
                allele = pd.Series(np.where(gb == snv.ref, 'ref', np.where(gb == snv.alt, 'alt', 'other')), index=gb.index)
                d = (snv.p0 - site['pos0']) * (1 if site['strand'] == '+' else -1)
                if d == 0:                                          # the SNV is the site: every read has a row there
                    sb = at_site['canonical_base'].str.upper().str.replace('U', 'T')
                    sg = sb if site['strand'] == '+' else sb.map(COMP)
                    allele = pd.Series(np.where(sg == snv.ref, 'ref', np.where(sg == snv.alt, 'alt', 'other')), index=sg.index)
                joined = at_site.join(allele.rename('allele'), how='inner')
                res = {}
                for al, genomic in (('ref', snv.ref), ('alt', snv.alt)):
                    r_ = joined[joined['allele'] == al]
                    ok = r_[~r_['fail']] if 'fail' in r_ else r_
                    calls = ok['call_code'].value_counts().to_dict()
                    tgt_base = ok['canonical_base'].str.upper().str.replace('U', 'T') == MOD_BASE[site['code']]
                    n = int(tgt_base.sum())
                    m = int(((ok['call_code'] == site['code']) & tgt_base).sum())
                    fp = allele_fp(table, fasta, site, snv.p0, genomic)
                    raw = 100 * m / n if n else np.nan
                    corr = max(0.0, raw - 100 * fp) if np.isfinite(raw) and np.isfinite(fp) else np.nan
                    res[al] = (n, m, raw, fp, corr)
                    rows.append((site['label'], s, snv.p0 + 1, d, f'{snv.ref}>{snv.alt}', snv.gt, al, len(r_),
                                 int(r_['fail'].sum()) if 'fail' in r_ else 0, n, m, raw, fp, corr,
                                 ' '.join(f'{k}:{v}' for k, v in sorted(calls.items(), key=lambda kv: -kv[1]))))
                (nr, mr, rr, fr, cr), (na_, ma, ra, fa, ca) = res['ref'], res['alt']
                if d == 0:                                          # different base on each allele: show every call
                    summ = {al: ' '.join(f'{k}:{v}' for k, v in joined[joined['allele'] == al]['call_code']
                                         .value_counts().items()) or '-' for al in ('ref', 'alt')}
                    print(f"  {s}  SNV AT the site {snv.ref}>{snv.alt} {snv.gt}:  ref-allele reads "
                          f"{int((joined['allele'] == 'ref').sum())} ({plain(site['code'])} {rr:5.1f}%; calls {summ['ref']})  |  "
                          f"alt-allele reads {int((joined['allele'] == 'alt').sum())} (calls {summ['alt']})")
                    continue
                p = (stats.fisher_exact([[mr, nr - mr], [ma, na_ - ma]])[1]
                     if min(nr, na_) >= a.min_reads else np.nan)
                kmer = 'in 9-mer' if abs(d) <= 4 else 'outside 9-mer'
                print(f"  {s}  SNV {snv.p0 + 1:,} ({d:+d} nt, {kmer}) {snv.ref}>{snv.alt} {snv.gt}:  "
                      f"ref {mr}/{nr} = {rr:5.1f}%  alt {ma}/{na_} = {ra:5.1f}%"
                      + (f"  [gIVT-corrected ref {cr:5.1f}%, alt {ca:5.1f}%]" if np.isfinite(cr) or np.isfinite(ca) else '')
                      + (f"  Fisher p = {p:.2g}" if np.isfinite(p) else f"  (< {a.min_reads} reads on an allele)"))
    cols = ['site', 'sample', 'snv_pos', 'snv_offset_rna', 'snv', 'gt', 'allele', 'reads', 'reads_fail_at_site',
            'n_valid_target_base', 'n_mod', 'occupancy_raw', 'gIVT_fp', 'occupancy_corrected', 'calls_at_site']
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(os.path.join(a.out, 'allele_counts.tsv'), sep='\t', index=False, float_format='%.4g')
    dz = pd.DataFrame(dosage, columns=['site', 'sample', 'cohort_occupancy', 'cohort_n_valid', 'snvs_within_4nt_rna_offset'])
    dz.to_csv(os.path.join(a.out, 'genotype_dosage.tsv'), sep='\t', index=False, float_format='%.4g')
    if len(df):                                                     # combined over heterozygous donors, per SNV
        print('\ncombined over heterozygous donors (CMH, donors as strata):')
        for (site, pos), g in df.groupby(['site', 'snv_pos']):
            strata, deltas = [], []
            for s, h in g.groupby('sample'):
                r_, l_ = h[h['allele'] == 'ref'], h[h['allele'] == 'alt']
                if len(r_) and len(l_) and min(r_['n_valid_target_base'].iat[0], l_['n_valid_target_base'].iat[0]) >= a.min_reads:
                    strata.append((np.array([r_['n_mod'].iat[0], l_['n_mod'].iat[0]]),
                                   np.array([r_['n_valid_target_base'].iat[0], l_['n_valid_target_base'].iat[0]])))
                    deltas.append(l_['occupancy_raw'].iat[0] - r_['occupancy_raw'].iat[0])
            if strata:
                q, df_, p = cmh_2xk(strata) if len(strata) else (np.nan, 0, np.nan)
                print(f"  {site} SNV {pos:,}: {len(strata)} donor(s), mean alt - ref {np.mean(deltas):+.1f} points, CMH p = {p:.2g}")
    print(f"\ntables: {a.out}/allele_counts.tsv (per donor and allele), {a.out}/genotype_dosage.tsv (every donor)")


if __name__ == '__main__':
    main()
