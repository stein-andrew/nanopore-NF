#!/usr/bin/env python3
"""Isoform- and gene-level modification levels at the final-stage sites.

    feature_mods.py positions --candidates S.cov20_fpcorr.bed --candidate-pct 20 > positions.bed
    feature_mods.py levels --counts counts.tsv --candidates S.cov20_fpcorr.bed --candidate-pct 20 \\
        --catalogue cohort.GRCh38.catalogue.gtf --classification cohort.GRCh38.classify.classification.txt.gz \\
        --reference genome.fa --errortable table.tsv.gz --givt-flags S.givt_flags.txt \\
        --min-cov 20 --min-pct 20 --prefix S

Candidates are the sample-level sites whose FP-corrected percent (column 19 of the gIVT
stage) is >= --candidate-pct. With the default 20 they are exactly the final-stage sites;
lowering it also finds sites modified in one isoform but diluted below 20% sample-wide.

levels, from bin/mod_feature_counts.sh's counts, per (isoform or gene, site, mod code):
  n_valid   passing calls at the site from that feature's reads whose base is the mod code's
            base (canonical + any mod code on it), like modkit pileup's Nvalid_cov: reads
            with another base there are pileup's Ndiff, not valid coverage
  n_mod     of those, calls of this mod code
  percent   100 * n_mod / n_valid, then FP-corrected by fp_correct.py's own functions (same
            9mer, mod code and thresholds as the sample-level stage)
In the BEDs: n_valid >= --min-cov and corrected percent >= --min-pct, the final-stage rule
applied per feature.

First, the 'all' level (every read) must reproduce the pileup's Nvalid_cov and Nmod at the
candidate sites. If fewer than --min-agreement of them match exactly, extract calls and
pileup are not counting the same calls, and this stops instead of writing numbers.

Writes:
  <prefix>.isoform_mods.bed            final-rule rows per isoform, annotated ('#' header)
  <prefix>.gene_mods.bed               the same per gene
  <prefix>.feature_mods.tsv.gz         every isoform/gene row, before the final rule
  <prefix>.feature_mods.summary.tsv    candidates, pileup agreement, rows per level/mod code
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np
import polars as pl

from fp_correct import MOD_TO_BASE, correct, load_table
from tabular import codes_of, read_tsv, read_tsv_header, write_tsv

SITE = ['chrom', 'pos', 'strand']
TX_ID = re.compile(r'transcript_id "([^"]+)"')
GENE_ID = re.compile(r'gene_id "([^"]+)"')
NO_NAME = {'', '.', 'NA', 'nan'}
BED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'mod_code', 'feature_id', 'gene_id',
            'ref_gene_name', 'ref_tx_id', 'structural_category', 'n_valid', 'n_mod',
            'percent_modified', 'percent_fpcorr', 'sample_percent_fpcorr']
COUNT_NAMES = ['level', 'feature', 'chrom', 'pos', 'strand', 'base', 'code', 'n']


def read_candidates(path, pct):
    df = read_tsv(path, [f'c{i}' for i in range(19)])
    df = df.with_columns(pl.col('c18').cast(pl.Float64, strict=False).alias('sample_pct_fpcorr'))
    df = df.filter(pl.col('sample_pct_fpcorr') >= pct)
    return df.select(pl.col('c0').alias('chrom'), pl.col('c1').cast(pl.Int64).alias('pos'),
                     pl.col('c2').cast(pl.Int64).alias('end'), pl.col('c3').alias('mod'), pl.col('c5').alias('strand'),
                     pl.col('c9').cast(pl.Int64).alias('sample_n_valid'), pl.col('c10').cast(pl.Float64).alias('sample_pct'),
                     pl.col('c11').cast(pl.Int64).alias('sample_n_mod'), 'sample_pct_fpcorr')


def read_counts(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            width = fh.readline().count('\t') + 1
        if width != len(COUNT_NAMES):
            sys.exit(f"feature_mods.py: {path} has {width} columns, expected {len(COUNT_NAMES)} "
                     f"({', '.join(COUNT_NAMES)}): bin/mod_feature_counts.sh and bin/feature_mods.py are from different versions")
    df = read_tsv(path, COUNT_NAMES)
    return df.with_columns(pl.col('pos').cast(pl.Int64), pl.col('n').cast(pl.Int64))


def annotations(catalogue, classification):
    """isoform -> gene, isoform -> (ref_gene_name, ref_tx_id, structural_category), gene -> names."""
    tx_gene = {}
    with open(catalogue) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 9 and cols[2] == 'transcript':
                tx, g = TX_ID.search(cols[8]), GENE_ID.search(cols[8])
                if tx and g:
                    tx_gene[tx.group(1)] = g.group(1)
    cl = read_tsv_header(classification)
    for col in ('isoform_id', 'structural_category', 'ref_gene_name', 'ref_tx_id'):
        if col not in cl.columns:
            sys.exit(f"feature_mods.py: {classification} has no {col} column")
    iso_info = {iso: (name or '.', tx or '.', cat or '.') for iso, name, tx, cat in zip(
        cl['isoform_id'].to_list(), cl['ref_gene_name'].to_list(), cl['ref_tx_id'].to_list(),
        cl['structural_category'].to_list())}
    gene_names = defaultdict(set)
    for tx, g in tx_gene.items():
        name = iso_info.get(tx, ('.',))[0]
        if name not in NO_NAME:
            gene_names[g].add(name)
    return tx_gene, iso_info, {g: ','.join(sorted(n)) for g, n in gene_names.items()}


def thresholds_from(path):
    thresholds = defaultdict(lambda: 0.7)
    for pair in open(path).read().split():
        mod, value = pair.split(',')
        thresholds[mod] = float(value)
    return thresholds


def fp_corrected(rows, reference, errortable, givt_flags):
    """Corrected percent for each row, by fp_correct.py's own functions (NaN where none)."""
    if rows.height == 0:
        return np.array([])
    chrom_codes, chroms = codes_of(rows['chrom'])
    mod_codes, mods = codes_of(rows['mod'])
    sites = dict(chrom_codes=chrom_codes, chroms=[str(c) for c in chroms],
                 mod_codes=mod_codes, mods=[str(m) for m in mods],
                 start=rows['pos'].to_numpy().astype(np.int64), minus=(rows['strand'] == '-').to_numpy(),
                 pct=rows['percent'].to_numpy().astype(np.float64))
    index, fps, present = load_table(errortable, thresholds_from(givt_flags))
    adjusted, _reason = correct(sites, reference, index, fps, present)
    return adjusted


def fmt(x):
    return '' if x is None or x != x else f"{x:.2f}"


def write_bed(rows, path):
    with open(path, 'w') as out:
        out.write('#' + '\t'.join(BED_COLS) + '\n')
        for r in rows.iter_rows(named=True):
            out.write('\t'.join(map(str, [
                r['chrom'], r['pos'], r['pos'] + 1, f"{r['feature']}:{r['mod']}", min(int(r['n_valid']), 1000), r['strand'],
                r['mod'], r['feature'], r['gene_id'], r['ref_gene_name'], r['ref_tx_id'], r['structural_category'],
                r['n_valid'], r['n_mod'], fmt(r['percent']), fmt(r['percent_fpcorr']), fmt(r['sample_pct_fpcorr'])])) + '\n')


def cmd_positions(args):
    cand = read_candidates(args.candidates, args.candidate_pct)
    pos = cand.select('chrom', 'pos', 'end').unique().sort(['chrom', 'pos', 'end'])
    sys.stdout.write(pos.write_csv(separator='\t', include_header=False, quote_style='never'))
    print(f"feature_mods.py: {cand.height} candidate site/mod rows at {pos.height} positions", file=sys.stderr)


def cmd_levels(args):
    cand = read_candidates(args.candidates, args.candidate_pct)
    bases = pl.DataFrame({'mod': list(MOD_TO_BASE), 'base': list(MOD_TO_BASE.values())})
    cand = cand.join(bases, on='mod', how='left')        # candidates are gIVT-corrected, so every code maps
    counts = read_counts(args.counts)

    # valid coverage per site counts only calls on the mod's own base, as pileup does
    valid = counts.group_by(['level', 'feature'] + SITE + ['base']).agg(pl.col('n').sum().alias('n_valid'))
    rows = valid.join(cand, on=SITE + ['base'], how='inner')
    modn = counts.rename({'code': 'mod', 'n': 'n_mod'})
    rows = rows.join(modn, on=['level', 'feature'] + SITE + ['base', 'mod'], how='left') \
               .with_columns(pl.col('n_mod').fill_null(0))

    # every read together must give the pileup's own numbers at the candidate sites
    allr = cand.join(rows.filter(pl.col('level') == 'all').select(SITE + ['mod', 'n_valid', 'n_mod']),
                     on=SITE + ['mod'], how='left')
    match = ((pl.col('n_valid') == pl.col('sample_n_valid')) & (pl.col('n_mod') == pl.col('sample_n_mod'))).fill_null(False)
    allr = allr.with_columns(match.alias('match'))
    agree = int(allr['match'].sum()) if cand.height else 0
    agreement = agree / cand.height if cand.height else 1.0
    if agreement < args.min_agreement:
        bad = allr.filter(~pl.col('match')).head(5).select(SITE + ['mod', 'sample_n_valid', 'n_valid', 'sample_n_mod', 'n_mod'])
        sys.exit(f"feature_mods.py: extract calls reproduce the pileup at only {agree} of {cand.height} candidate sites "
                 f"({agreement:.1%} < {args.min_agreement:.0%}): they are not counting the same calls. First mismatches "
                 f"(pileup sample_n_valid/sample_n_mod vs extract n_valid/n_mod):\n"
                 + bad.write_csv(separator='\t', quote_style='never'))

    feat = rows.filter(pl.col('level').is_in(['isoform', 'gene'])) \
               .with_columns((100 * pl.col('n_mod') / pl.col('n_valid')).alias('percent'))
    adjusted = fp_corrected(feat, args.reference, args.errortable, args.givt_flags)
    feat = feat.with_columns(pl.Series('percent_fpcorr', adjusted, dtype=pl.Float64).fill_nan(None))
    feat = feat.with_columns(((pl.col('n_valid') >= args.min_cov) & (pl.col('percent_fpcorr') >= args.min_pct))
                             .fill_null(False).alias('final'))

    tx_gene, iso_info, gene_names = annotations(args.catalogue, args.classification)
    iso_df = pl.DataFrame({'feature': list(iso_info), 'iso_name': [v[0] for v in iso_info.values()],
                           'iso_tx': [v[1] for v in iso_info.values()], 'iso_cat': [v[2] for v in iso_info.values()]},
                          schema={'feature': pl.Utf8, 'iso_name': pl.Utf8, 'iso_tx': pl.Utf8, 'iso_cat': pl.Utf8})
    txg_df = pl.DataFrame({'feature': list(tx_gene), 'tx_gene': list(tx_gene.values())},
                          schema={'feature': pl.Utf8, 'tx_gene': pl.Utf8})
    gn_df = pl.DataFrame({'feature': list(gene_names), 'gene_name': [v or '.' for v in gene_names.values()]},
                         schema={'feature': pl.Utf8, 'gene_name': pl.Utf8})
    feat = feat.join(iso_df, on='feature', how='left').join(txg_df, on='feature', how='left') \
               .join(gn_df, on='feature', how='left')
    is_iso = pl.col('level') == 'isoform'
    feat = feat.with_columns(
        pl.when(is_iso).then(pl.col('tx_gene').fill_null('.')).otherwise(pl.col('feature')).alias('gene_id'),
        pl.when(is_iso).then(pl.col('iso_name').fill_null('.')).otherwise(pl.col('gene_name').fill_null('.')).alias('ref_gene_name'),
        pl.when(is_iso).then(pl.col('iso_tx').fill_null('.')).otherwise(pl.lit('.')).alias('ref_tx_id'),
        pl.when(is_iso).then(pl.col('iso_cat').fill_null('.')).otherwise(pl.lit('.')).alias('structural_category'),
    ).sort(['chrom', 'pos', 'strand', 'level', 'feature', 'mod'])

    p = args.prefix
    write_bed(feat.filter(is_iso & pl.col('final')), f"{p}.isoform_mods.bed")
    write_bed(feat.filter((pl.col('level') == 'gene') & pl.col('final')), f"{p}.gene_mods.bed")
    full = feat.select('level', 'feature', 'gene_id', 'ref_gene_name', 'ref_tx_id', 'structural_category', 'chrom', 'pos',
                       'strand', 'mod', 'n_valid', 'n_mod', 'percent', 'percent_fpcorr', 'sample_pct_fpcorr',
                       pl.col('final').cast(pl.Int64))
    write_tsv(full, f"{p}.feature_mods.tsv.gz", float_precision=4)

    with open(f"{p}.feature_mods.summary.tsv", 'w') as out:
        out.write("metric\tvalue\n")
        out.write(f"candidate_site_mod_rows\t{cand.height}\n")
        out.write(f"pileup_agreement\t{agree}/{cand.height}\n")
        for level in ('isoform', 'gene'):
            sub = feat.filter(pl.col('level') == level)
            out.write(f"{level}_rows\t{sub.height}\n")
            out.write(f"{level}_features\t{sub['feature'].n_unique()}\n")
            fin = sub.filter(pl.col('final')).group_by('mod').agg(pl.len().alias('n')).sort('mod')
            for mod, n in fin.iter_rows():
                out.write(f"{level}_final:{mod}\t{n}\n")
    print(f"feature_mods.py: pileup agreement {agree}/{cand.height}; final rows: isoform "
          f"{feat.filter(is_iso & pl.col('final')).height}, gene {feat.filter((pl.col('level') == 'gene') & pl.col('final')).height}",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('positions', 'levels'):
        s = sub.add_parser(name)
        s.add_argument('--candidates', required=True, help='sample-level gIVT stage: S.cov20_fpcorr.bed (19 columns)')
        s.add_argument('--candidate-pct', type=float, required=True, help='candidates: sample-level corrected percent >= this')
        if name == 'levels':
            s.add_argument('--counts', required=True, help='bin/mod_feature_counts.sh output')
            s.add_argument('--catalogue', required=True, help='cohort catalogue GTF (transcript -> gene)')
            s.add_argument('--classification', required=True, help='isomatch classification of the catalogue')
            s.add_argument('--reference', required=True, help='uncompressed reference FASTA with .fai')
            s.add_argument('--errortable', required=True, help='gIVT error table')
            s.add_argument('--givt-flags', required=True, help='S.givt_flags.txt (the sample-level thresholds)')
            s.add_argument('--min-cov', type=int, required=True)
            s.add_argument('--min-pct', type=float, required=True)
            s.add_argument('--min-agreement', type=float, default=0.95,
                           help='fraction of candidate sites where extract calls must reproduce the pileup (default 0.95)')
            s.add_argument('--prefix', required=True)
    args = ap.parse_args()
    cmd_positions(args) if args.cmd == 'positions' else cmd_levels(args)


if __name__ == '__main__':
    main()
