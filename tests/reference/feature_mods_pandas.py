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
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from fp_correct_pandas import MOD_TO_BASE, correct, load_table

SITE = ['chrom', 'pos', 'strand']
TX_ID = re.compile(r'transcript_id "([^"]+)"')
GENE_ID = re.compile(r'gene_id "([^"]+)"')
NO_NAME = {'', '.', 'NA', 'nan'}
BED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'mod_code', 'feature_id', 'gene_id',
            'ref_gene_name', 'ref_tx_id', 'structural_category', 'n_valid', 'n_mod',
            'percent_modified', 'percent_fpcorr', 'sample_percent_fpcorr']


def read_candidates(path, pct):
    cols = ['chrom', 'pos', 'end', 'mod', 'strand', 'sample_n_valid', 'sample_pct', 'sample_n_mod', 'sample_pct_fpcorr']
    if os.path.getsize(path) == 0:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path, sep='\t', header=None, dtype=str, keep_default_na=False, na_filter=False,
                     quoting=csv.QUOTE_NONE)
    if df.shape[1] != 19:
        sys.exit(f"feature_mods.py: {path} has {df.shape[1]} columns, expected 19 (bedMethyl + FP-corrected percent)")
    corrected = pd.to_numeric(df[18], errors='coerce')
    c, corrected = df[corrected >= pct], corrected[corrected >= pct]
    return pd.DataFrame({'chrom': c[0], 'pos': c[1].astype(np.int64), 'end': c[2].astype(np.int64),
                         'mod': c[3], 'strand': c[5], 'sample_n_valid': c[9].astype(np.int64),
                         'sample_pct': c[10].astype(np.float64), 'sample_n_mod': c[11].astype(np.int64),
                         'sample_pct_fpcorr': corrected}).reset_index(drop=True)


def read_counts(path):
    names = ['level', 'feature', 'chrom', 'pos', 'strand', 'base', 'code', 'n']
    if os.path.getsize(path) == 0:
        return pd.DataFrame(columns=names).astype({'pos': np.int64, 'n': np.int64})
    with open(path) as fh:
        width = fh.readline().count('\t') + 1
    if width != len(names):
        sys.exit(f"feature_mods.py: {path} has {width} columns, expected {len(names)} "
                 f"({', '.join(names)}): bin/mod_feature_counts.sh and bin/feature_mods.py are from different versions")
    return pd.read_csv(path, sep='\t', header=None, names=names, keep_default_na=False, na_filter=False,
                       dtype={'level': str, 'feature': str, 'chrom': str, 'pos': np.int64,
                              'strand': str, 'base': str, 'code': str, 'n': np.int64})


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
    opener = gzip.open if classification.endswith('.gz') else open
    with opener(classification, 'rt') as fh:
        cl = pd.read_csv(fh, sep='\t', dtype=str, keep_default_na=False, na_filter=False)
    for col in ('isoform_id', 'structural_category', 'ref_gene_name', 'ref_tx_id'):
        if col not in cl.columns:
            sys.exit(f"feature_mods.py: {classification} has no {col} column")
    iso_info = {r.isoform_id: (r.ref_gene_name or '.', r.ref_tx_id or '.', r.structural_category or '.')
                for r in cl.itertuples(index=False)}
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
    """Corrected percent for each row, by fp_correct.py's own functions."""
    if rows.empty:
        return np.array([])
    chrom_codes, chroms = pd.factorize(rows['chrom'])
    mod_codes, mods = pd.factorize(rows['mod'])
    sites = dict(chrom_codes=np.asarray(chrom_codes, dtype=np.int64), chroms=[str(c) for c in chroms],
                 mod_codes=np.asarray(mod_codes, dtype=np.int64), mods=[str(m) for m in mods],
                 start=rows['pos'].to_numpy(dtype=np.int64), minus=(rows['strand'] == '-').to_numpy(),
                 pct=rows['percent'].to_numpy(dtype=np.float64))
    index, fps, present = load_table(errortable, thresholds_from(givt_flags))
    adjusted, _reason = correct(sites, reference, index, fps, present)
    return adjusted


def fmt(x):
    return '' if x != x else f"{x:.2f}"


def write_bed(rows, path):
    with open(path, 'w') as out:
        out.write('#' + '\t'.join(BED_COLS) + '\n')
        for r in rows.itertuples(index=False):
            out.write('\t'.join(map(str, [
                r.chrom, r.pos, r.pos + 1, f"{r.feature}:{r.mod}", min(int(r.n_valid), 1000), r.strand,
                r.mod, r.feature, r.gene_id, r.ref_gene_name, r.ref_tx_id, r.structural_category,
                r.n_valid, r.n_mod, fmt(r.percent), fmt(r.percent_fpcorr), fmt(r.sample_pct_fpcorr)])) + '\n')


def cmd_positions(args):
    cand = read_candidates(args.candidates, args.candidate_pct)
    pos = cand[['chrom', 'pos', 'end']].drop_duplicates().sort_values(['chrom', 'pos'])
    pos.to_csv(sys.stdout, sep='\t', header=False, index=False)
    print(f"feature_mods.py: {len(cand)} candidate site/mod rows at {len(pos)} positions", file=sys.stderr)


def cmd_levels(args):
    cand = read_candidates(args.candidates, args.candidate_pct)
    cand['base'] = cand['mod'].map(MOD_TO_BASE)      # candidates are gIVT-corrected, so every code maps
    counts = read_counts(args.counts)

    # valid coverage per site counts only calls on the mod's own base, as pileup does
    valid = counts.groupby(['level', 'feature'] + SITE + ['base'], as_index=False)['n'].sum() \
                  .rename(columns={'n': 'n_valid'})
    rows = valid.merge(cand, on=SITE + ['base'])
    modn = counts.rename(columns={'code': 'mod', 'n': 'n_mod'})
    rows = rows.merge(modn, on=['level', 'feature'] + SITE + ['base', 'mod'], how='left')
    rows['n_mod'] = rows['n_mod'].fillna(0).astype(np.int64)

    # every read together must give the pileup's own numbers at the candidate sites
    allr = cand.merge(rows[rows['level'] == 'all'][SITE + ['mod', 'n_valid', 'n_mod']], on=SITE + ['mod'], how='left')
    agree = ((allr['n_valid'] == allr['sample_n_valid']) & (allr['n_mod'] == allr['sample_n_mod'])).sum()
    agreement = agree / len(cand) if len(cand) else 1.0
    if agreement < args.min_agreement:
        bad = allr[~((allr['n_valid'] == allr['sample_n_valid']) & (allr['n_mod'] == allr['sample_n_mod']))].head(5)
        sys.exit(f"feature_mods.py: extract calls reproduce the pileup at only {agree} of {len(cand)} candidate sites "
                 f"({agreement:.1%} < {args.min_agreement:.0%}): they are not counting the same calls. First mismatches "
                 f"(pileup sample_n_valid/sample_n_mod vs extract n_valid/n_mod):\n"
                 + bad[SITE + ['mod', 'sample_n_valid', 'n_valid', 'sample_n_mod', 'n_mod']].to_string(index=False))

    feat = rows[rows['level'].isin(['isoform', 'gene'])].copy()
    feat['percent'] = 100 * feat['n_mod'] / feat['n_valid']
    feat['percent_fpcorr'] = fp_corrected(feat, args.reference, args.errortable, args.givt_flags)
    feat['final'] = (feat['n_valid'] >= args.min_cov) & (feat['percent_fpcorr'] >= args.min_pct)

    tx_gene, iso_info, gene_names = annotations(args.catalogue, args.classification)
    iso = feat['level'] == 'isoform'
    info = feat['feature'].map(lambda f: iso_info.get(f, ('.', '.', '.')))
    feat['gene_id'] = np.where(iso, feat['feature'].map(lambda f: tx_gene.get(f, '.')), feat['feature'])
    feat['ref_gene_name'] = np.where(iso, info.str[0], feat['feature'].map(lambda g: gene_names.get(g, '.') or '.'))
    feat['ref_tx_id'] = np.where(iso, info.str[1], '.')
    feat['structural_category'] = np.where(iso, info.str[2], '.')
    feat = feat.sort_values(['chrom', 'pos', 'strand', 'level', 'feature', 'mod']).reset_index(drop=True)

    p = args.prefix
    write_bed(feat[(feat['level'] == 'isoform') & feat['final']], f"{p}.isoform_mods.bed")
    write_bed(feat[(feat['level'] == 'gene') & feat['final']], f"{p}.gene_mods.bed")
    full = feat[['level', 'feature', 'gene_id', 'ref_gene_name', 'ref_tx_id', 'structural_category', 'chrom', 'pos',
                 'strand', 'mod', 'n_valid', 'n_mod', 'percent', 'percent_fpcorr', 'sample_pct_fpcorr', 'final']].copy()
    full['final'] = full['final'].astype(int)
    full.to_csv(f"{p}.feature_mods.tsv.gz", sep='\t', index=False, float_format='%.4f', compression='gzip')

    with open(f"{p}.feature_mods.summary.tsv", 'w') as out:
        out.write("metric\tvalue\n")
        out.write(f"candidate_site_mod_rows\t{len(cand)}\n")
        out.write(f"pileup_agreement\t{agree}/{len(cand)}\n")
        for level in ('isoform', 'gene'):
            sub = feat[feat['level'] == level]
            out.write(f"{level}_rows\t{len(sub)}\n")
            out.write(f"{level}_features\t{sub['feature'].nunique()}\n")
            for mod, n in sub[sub['final']].groupby('mod').size().items():
                out.write(f"{level}_final:{mod}\t{n}\n")
    print(f"feature_mods.py: pileup agreement {agree}/{len(cand)}; final rows: isoform "
          f"{int(((feat['level'] == 'isoform') & feat['final']).sum())}, gene {int(((feat['level'] == 'gene') & feat['final']).sum())}",
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
