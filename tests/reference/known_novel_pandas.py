#!/usr/bin/env python3
"""Known vs novel for Clair3-RNA calls (step 6), by exact allele match after normalisation.

    known_novel.py --calls S.calls.tsv --gnomad S.gnomad.tsv --mito S.mito.tsv --clinvar S.clinvar.tsv \\
        --mito-contig chrM --dp-split 10 --prefix S

Inputs are bcftools query tables (no header), the references already restricted to the
call positions and everything split to one ALT per row and left-aligned:
  calls    CHROM POS REF ALT FILTER QUAL GT GQ DP AF        Clair3-RNA
  gnomad   CHROM POS REF ALT FILTER AF_joint AF_exomes AF_genomes nhomalt_joint   v4.1 joint
  mito     CHROM POS REF ALT FILTER AF_hom AF_het            v3.1 chrM
  clinvar  CHROM POS REF ALT CLNSIG CLNREVSTAT CLNDN GENEINFO

Class of each call, in this order:
  refcall  Clair3-RNA RefCall (no variant called)
  lowqual  FILTER LowQual
  editing  FILTER RNAEditing: a REDIportal site. Split off before novelty, as SALRR does.
  known    PASS, and the same CHROM, POS, REF and ALT are in gnomAD and passed QC there:
           nuclear (joint) FILTER is not BOTH_FILTERED, i.e. passed in exomes or genomes;
           chrM FILTER is PASS
  novel    PASS and not known
A gnomAD record at the same position with another ALT, or one that failed QC everywhere,
does not make a call known. ClinVar significance is attached wherever CHROM/POS/REF/ALT
match, whatever the class.

Writes:
  <prefix>.variants.tsv.gz      one row per call: the call, its class, gnomAD frequencies and
                                filter, ClinVar fields
  <prefix>.variants.summary.tsv per class: calls, SNVs, indels, transitions, transversions,
                                Ti/Tv, A>G+T>C share of SNVs, ClinVar matches; novel also split
                                by depth (DP >= --dp-split vs below)
  <prefix>.annot.tsv, <prefix>.annot.hdr
                                CHROM POS REF ALT + INFO values and their header lines, for
                                bcftools annotate to write the classes into the VCF
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

KEY = ['chrom', 'pos', 'ref', 'alt']
CALL_COLS = KEY + ['filter', 'qual', 'gt', 'gq', 'dp', 'vaf']
GNOMAD_COLS = KEY + ['gnomad_filter', 'gnomad_af', 'gnomad_af_exomes', 'gnomad_af_genomes', 'gnomad_nhomalt']
MITO_COLS = KEY + ['mito_filter', 'gnomad_af_hom', 'gnomad_af_het']
CLINVAR_COLS = KEY + ['clinvar_sig', 'clinvar_revstat', 'clinvar_dn', 'clinvar_gene']
CLASSES = ['known', 'novel', 'editing', 'lowqual', 'refcall']
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

ANNOT = [  # (INFO id, column, Number, Type, Description)
    ('VARIANT_CLASS', 'class', '1', 'String', 'known, novel, editing or lowqual (bin/known_novel.py)'),
    ('gnomAD_AF', 'gnomad_af', '1', 'Float', 'gnomAD v4.1 joint allele frequency (chrM: v3.1 AF_hom)'),
    ('gnomAD_AF_exomes', 'gnomad_af_exomes', '1', 'Float', 'gnomAD v4.1 exomes allele frequency'),
    ('gnomAD_AF_genomes', 'gnomad_af_genomes', '1', 'Float', 'gnomAD v4.1 genomes allele frequency'),
    ('gnomAD_FILTER', 'gnomad_filter', '1', 'String', 'gnomAD FILTER of the matching record'),
    ('CLNSIG', 'clinvar_sig', '.', 'String', 'ClinVar clinical significance of the matching record'),
    ('CLNREVSTAT', 'clinvar_revstat', '.', 'String', 'ClinVar review status of the matching record'),
]
MASK_ANNOT = ('MASKED', 'masked', '1', 'String', 'IG, TR or HLA region (somatic or alignment-distorted variation), else .')


def call_cols(path):
    """The calls table may carry an 11th column, masked (IG / TR / HLA / .)."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            if fh.readline().count('\t') + 1 == len(CALL_COLS) + 1:
                return CALL_COLS + ['masked']
    return CALL_COLS


def read(path, cols):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path, sep='\t', header=None, dtype=str, keep_default_na=False, na_filter=False)
    if df.shape[1] != len(cols):
        sys.exit(f"known_novel.py: {path} has {df.shape[1]} columns, expected {len(cols)} ({' '.join(cols)})")
    df.columns = cols
    df['pos'] = df['pos'].astype(np.int64)
    return df.drop_duplicates(subset=KEY, keep='first')


def classify(calls, gnomad, mito, mito_contig):
    nuclear_ok = gnomad[gnomad['gnomad_filter'] != 'BOTH_FILTERED']
    mito_ok = mito[mito['mito_filter'] == 'PASS']
    c = calls.merge(gnomad, on=KEY, how='left').merge(mito, on=KEY, how='left')
    is_mito = c['chrom'] == mito_contig
    known_nuc = ~is_mito & c[KEY].apply(tuple, axis=1).isin(set(map(tuple, nuclear_ok[KEY].values)))
    known_mito = is_mito & c[KEY].apply(tuple, axis=1).isin(set(map(tuple, mito_ok[KEY].values)))
    # chrM rows show the mito filter and homoplasmic AF in the shared columns
    c['gnomad_filter'] = np.where(is_mito, c['mito_filter'], c['gnomad_filter'])
    c['gnomad_af'] = np.where(is_mito, c['gnomad_af_hom'], c['gnomad_af'])
    c = c.drop(columns=['mito_filter'])
    flt = c['filter']
    c['class'] = np.select(
        [flt == 'RefCall', flt == 'LowQual', flt == 'RNAEditing', (flt == 'PASS') & (known_nuc | known_mito), flt == 'PASS'],
        ['refcall', 'lowqual', 'editing', 'known', 'novel'], default='other')
    for col in GNOMAD_COLS[4:] + MITO_COLS[5:]:
        c[col] = c[col].fillna('.').replace('', '.')
    return c


def substitution_stats(d):
    snv = d[(d['ref'].str.len() == 1) & (d['alt'].str.len() == 1)]
    pairs = list(zip(snv['ref'], snv['alt']))
    ti = sum(p in TRANSITIONS for p in pairs)
    tv = len(pairs) - ti
    agtc = sum(p in {('A', 'G'), ('T', 'C')} for p in pairs)
    return dict(calls=len(d), snv=len(snv), indel=len(d) - len(snv), transitions=ti, transversions=tv,
                ti_tv=f"{ti / tv:.2f}" if tv else 'NA', ag_tc_share=f"{agtc / len(snv):.3f}" if len(snv) else 'NA',
                in_clinvar=int((d['clinvar_sig'] != '.').sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--calls', required=True)
    ap.add_argument('--gnomad', required=True)
    ap.add_argument('--mito', required=True)
    ap.add_argument('--clinvar', required=True)
    ap.add_argument('--mito-contig', default='chrM')
    ap.add_argument('--dp-split', type=int, default=10, help='summary splits novel calls at DP >= this')
    ap.add_argument('--prefix', required=True)
    args = ap.parse_args()

    calls = read(args.calls, call_cols(args.calls))
    has_mask = 'masked' in calls.columns
    gnomad, mito = read(args.gnomad, GNOMAD_COLS), read(args.mito, MITO_COLS)
    clinvar = read(args.clinvar, CLINVAR_COLS)
    c = classify(calls, gnomad, mito, args.mito_contig)
    c = c.merge(clinvar, on=KEY, how='left')
    for col in CLINVAR_COLS[4:]:
        c[col] = c[col].fillna('.').replace('', '.')
    other = c[c['class'] == 'other']
    if len(other):
        print(f"known_novel.py: {len(other)} call(s) with FILTER values not handled: {sorted(other['filter'].unique())}",
              file=sys.stderr)

    p = args.prefix
    out_cols = CALL_COLS + (['masked'] if has_mask else []) + ['class', 'gnomad_af', 'gnomad_af_exomes', 'gnomad_af_genomes', 'gnomad_nhomalt',
                            'gnomad_af_het', 'gnomad_filter'] + CLINVAR_COLS[4:]
    c[out_cols].to_csv(f"{p}.variants.tsv.gz", sep='\t', index=False, compression='gzip')

    rows = []
    base = c[c['masked'] == '.'] if has_mask else c        # summary rows count unmasked calls
    for cls in CLASSES + (['other'] if len(other) else []):
        rows.append(dict(cls=cls, **substitution_stats(base[base['class'] == cls])))
    novel = base[base['class'] == 'novel']
    dp = pd.to_numeric(novel['dp'], errors='coerce')
    rows.append(dict(cls=f"novel_dp>={args.dp_split}", **substitution_stats(novel[dp >= args.dp_split])))
    rows.append(dict(cls=f"novel_dp<{args.dp_split}", **substitution_stats(novel[~(dp >= args.dp_split)])))
    if has_mask:
        for cls in CLASSES:
            rows.append(dict(cls=f"masked:{cls}", **substitution_stats(c[(c['masked'] != '.') & (c['class'] == cls)])))
    summary = pd.DataFrame(rows).rename(columns={'cls': 'class'})
    summary.to_csv(f"{p}.variants.summary.tsv", sep='\t', index=False)

    annot_spec = ANNOT + ([MASK_ANNOT] if has_mask else [])
    annot = c[KEY + [col for _, col, _, _, _ in annot_spec]].copy()
    annot = annot[annot['class'] != 'refcall']
    annot.to_csv(f"{p}.annot.tsv", sep='\t', index=False, header=False)
    with open(f"{p}.annot.hdr", 'w') as fh:
        for info_id, _, number, typ, desc in annot_spec:
            fh.write(f'##INFO=<ID={info_id},Number={number},Type={typ},Description="{desc}">\n')

    counts = c['class'].value_counts()
    print("known_novel.py: " + ', '.join(f"{k} {counts.get(k, 0)}" for k in CLASSES), file=sys.stderr)


if __name__ == '__main__':
    main()
