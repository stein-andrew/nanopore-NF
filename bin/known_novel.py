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
import sys

import polars as pl

from tabular import read_tsv, write_tsv

KEY = ['chrom', 'pos', 'ref', 'alt']
CALL_COLS = KEY + ['filter', 'qual', 'gt', 'gq', 'dp', 'vaf']
GNOMAD_COLS = KEY + ['gnomad_filter', 'gnomad_af', 'gnomad_af_exomes', 'gnomad_af_genomes', 'gnomad_nhomalt']
MITO_COLS = KEY + ['mito_filter', 'gnomad_af_hom', 'gnomad_af_het']
CLINVAR_COLS = KEY + ['clinvar_sig', 'clinvar_revstat', 'clinvar_dn', 'clinvar_gene']
CLASSES = ['known', 'novel', 'editing', 'lowqual', 'refcall']
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
OUT_COLS = CALL_COLS + ['class', 'gnomad_af', 'gnomad_af_exomes', 'gnomad_af_genomes', 'gnomad_nhomalt',
                        'gnomad_af_het', 'gnomad_filter'] + CLINVAR_COLS[4:]

ANNOT = [  # (INFO id, column, Number, Type, Description)
    ('VARIANT_CLASS', 'class', '1', 'String', 'known, novel, editing or lowqual (bin/known_novel.py)'),
    ('gnomAD_AF', 'gnomad_af', '1', 'Float', 'gnomAD v4.1 joint allele frequency (chrM: v3.1 AF_hom)'),
    ('gnomAD_AF_exomes', 'gnomad_af_exomes', '1', 'Float', 'gnomAD v4.1 exomes allele frequency'),
    ('gnomAD_AF_genomes', 'gnomad_af_genomes', '1', 'Float', 'gnomAD v4.1 genomes allele frequency'),
    ('gnomAD_FILTER', 'gnomad_filter', '1', 'String', 'gnomAD FILTER of the matching record'),
    ('CLNSIG', 'clinvar_sig', '.', 'String', 'ClinVar clinical significance of the matching record'),
    ('CLNREVSTAT', 'clinvar_revstat', '.', 'String', 'ClinVar review status of the matching record'),
]


def read(path, cols):
    df = read_tsv(path, cols).with_columns(pl.col('pos').cast(pl.Int64))
    return df.unique(subset=KEY, keep='first', maintain_order=True)


def classify(calls, gnomad, mito, clinvar, mito_contig):
    nuclear_ok = gnomad.filter(pl.col('gnomad_filter') != 'BOTH_FILTERED').select(KEY).with_columns(pl.lit(True).alias('_kn'))
    mito_ok = mito.filter(pl.col('mito_filter') == 'PASS').select(KEY).with_columns(pl.lit(True).alias('_km'))
    c = calls.with_row_index('_i') \
             .join(gnomad, on=KEY, how='left').join(mito, on=KEY, how='left') \
             .join(nuclear_ok, on=KEY, how='left').join(mito_ok, on=KEY, how='left') \
             .join(clinvar, on=KEY, how='left').sort('_i')
    is_mito = pl.col('chrom') == mito_contig
    known = (~is_mito & pl.col('_kn').fill_null(False)) | (is_mito & pl.col('_km').fill_null(False))
    flt = pl.col('filter')
    c = c.with_columns(
        # chrM rows show the mito filter and homoplasmic AF in the shared columns
        pl.when(is_mito).then(pl.col('mito_filter')).otherwise(pl.col('gnomad_filter')).alias('gnomad_filter'),
        pl.when(is_mito).then(pl.col('gnomad_af_hom')).otherwise(pl.col('gnomad_af')).alias('gnomad_af'),
        pl.when(flt == 'RefCall').then(pl.lit('refcall'))
          .when(flt == 'LowQual').then(pl.lit('lowqual'))
          .when(flt == 'RNAEditing').then(pl.lit('editing'))
          .when((flt == 'PASS') & known).then(pl.lit('known'))
          .when(flt == 'PASS').then(pl.lit('novel'))
          .otherwise(pl.lit('other')).alias('class'))
    fill = GNOMAD_COLS[4:] + ['gnomad_af_het'] + CLINVAR_COLS[4:]
    return c.with_columns([pl.col(x).fill_null('.') for x in fill])


def substitution_stats(d):
    refs, alts, sigs = d['ref'].to_list(), d['alt'].to_list(), d['clinvar_sig'].to_list()
    pairs = [(r, a) for r, a in zip(refs, alts) if len(r) == 1 and len(a) == 1]
    ti = sum(p in TRANSITIONS for p in pairs)
    tv = len(pairs) - ti
    agtc = sum(p in {('A', 'G'), ('T', 'C')} for p in pairs)
    return dict(calls=d.height, snv=len(pairs), indel=d.height - len(pairs), transitions=ti, transversions=tv,
                ti_tv=f"{ti / tv:.2f}" if tv else 'NA', ag_tc_share=f"{agtc / len(pairs):.3f}" if pairs else 'NA',
                in_clinvar=sum(s != '.' for s in sigs))


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

    c = classify(read(args.calls, CALL_COLS), read(args.gnomad, GNOMAD_COLS), read(args.mito, MITO_COLS),
                 read(args.clinvar, CLINVAR_COLS), args.mito_contig)
    other = c.filter(pl.col('class') == 'other')
    if other.height:
        print(f"known_novel.py: {other.height} call(s) with FILTER values not handled: {sorted(set(other['filter'].to_list()))}",
              file=sys.stderr)

    p = args.prefix
    write_tsv(c.select(OUT_COLS), f"{p}.variants.tsv.gz")

    rows = [dict(cls=cls, **substitution_stats(c.filter(pl.col('class') == cls)))
            for cls in CLASSES + (['other'] if other.height else [])]
    novel = c.filter(pl.col('class') == 'novel')
    deep = pl.col('dp').cast(pl.Float64, strict=False) >= args.dp_split
    rows.append(dict(cls=f"novel_dp>={args.dp_split}", **substitution_stats(novel.filter(deep.fill_null(False)))))
    rows.append(dict(cls=f"novel_dp<{args.dp_split}", **substitution_stats(novel.filter(~deep.fill_null(False)))))
    with open(f"{p}.variants.summary.tsv", 'w') as out:
        names = ['class', 'calls', 'snv', 'indel', 'transitions', 'transversions', 'ti_tv', 'ag_tc_share', 'in_clinvar']
        out.write('\t'.join(names) + '\n')
        for r in rows:
            out.write('\t'.join(str(r['cls'] if n == 'class' else r[n]) for n in names) + '\n')

    annot = c.filter(pl.col('class') != 'refcall').select(KEY + [col for _, col, _, _, _ in ANNOT])
    write_tsv(annot, f"{p}.annot.tsv", header=False)
    with open(f"{p}.annot.hdr", 'w') as fh:
        for info_id, _, number, typ, desc in ANNOT:
            fh.write(f'##INFO=<ID={info_id},Number={number},Type={typ},Description="{desc}">\n')

    counts = dict(c.group_by('class').agg(pl.len()).iter_rows())
    print("known_novel.py: " + ', '.join(f"{k} {counts.get(k, 0)}" for k in CLASSES), file=sys.stderr)


if __name__ == '__main__':
    main()
