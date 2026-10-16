#!/usr/bin/env python3
"""Cohort high-confidence filter (SALRR): per-sample CPM from IsoQuant requantification
against the cohort catalogue; keep transcripts whose median CPM across samples > --min.

    cohort_cpm.py --gtf cohort.GRCh38.catalogue.gtf \
        --classification cohort.GRCh38.classify.classification.txt.gz \
        --counts ST1918=ST1918.transcript_counts.tsv --counts ST1922=ST1922.transcript_counts.tsv \
        --min 1 --prefix cohort.GRCh38

CPM is filter_models.py's, imported, not copied: count / total assigned reads in that
sample x 1e6, IsoQuant's __ summary rows excluded. A catalogue transcript missing from a
sample's counts is 0 there. A count row naming a transcript that is not in the catalogue
is an error: it means the requantification and the catalogue do not match.

Writes:
  <prefix>.counts.tsv                          transcript x sample read counts
  <prefix>.cpm.tsv                             transcript x sample CPM, median_cpm, kept
  <prefix>.high_confidence.gtf                 catalogue lines of kept transcripts (+ genes)
  <prefix>.high_confidence.classification.tsv  classification rows of kept transcripts
  <prefix>.summary.tsv                         assigned reads per sample; transcripts and
                                               structural categories, catalogue vs kept
"""
import argparse
import csv
import statistics
import sys
from collections import Counter

from filter_models import cpm, open_text, read_counts, read_gtf, write_subset


def parse_counts_args(items):
    samples = []
    for item in items:
        sample, sep, path = item.partition('=')
        if not sep or not sample or not path:
            sys.exit(f"--counts expects SAMPLE=PATH, got: {item}")
        if sample in (s for s, _ in samples):
            sys.exit(f"--counts: sample {sample} given twice")
        samples.append((sample, path))
    return samples


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', required=True, help='cohort catalogue GTF (with gene lines)')
    ap.add_argument('--classification', required=True, help='isomatch classify classification.txt[.gz]')
    ap.add_argument('--counts', action='append', required=True, metavar='SAMPLE=PATH',
                    help='IsoQuant transcript_counts.tsv from requantification; repeat per sample')
    ap.add_argument('--min', type=float, default=1.0, help='keep transcripts with median CPM > MIN (default 1)')
    ap.add_argument('--prefix', required=True)
    args = ap.parse_args()

    lines, tx_gene = read_gtf(args.gtf)
    catalogue = list(tx_gene)                     # catalogue order
    samples = sorted(parse_counts_args(args.counts))   # stable column order, whatever order tasks finished in

    counts, cpms, totals = {}, {}, {}
    for sample, path in samples:
        c = read_counts(path)
        unknown = sorted(set(c) - set(tx_gene))
        if unknown:
            sys.exit(f"{path}: {len(unknown)} transcript(s) not in the catalogue, e.g. {unknown[:3]}")
        counts[sample] = c
        cpms[sample], totals[sample] = cpm(c, path)

    names = [s for s, _ in samples]
    median = {tx: statistics.median(cpms[s].get(tx, 0.0) for s in names) for tx in catalogue}
    keep = {tx for tx in catalogue if median[tx] > args.min}

    with open(f"{args.prefix}.counts.tsv", 'w') as out:
        out.write('transcript_id\tgene_id\t' + '\t'.join(names) + '\n')
        for tx in catalogue:
            out.write(f"{tx}\t{tx_gene[tx] or ''}\t" + '\t'.join(f"{counts[s].get(tx, 0.0):.10g}" for s in names) + '\n')

    with open(f"{args.prefix}.cpm.tsv", 'w') as out:
        out.write('transcript_id\tgene_id\t' + '\t'.join(names) + '\tmedian_cpm\tkept\n')
        for tx in catalogue:
            row = '\t'.join(f"{cpms[s].get(tx, 0.0):.4f}" for s in names)
            out.write(f"{tx}\t{tx_gene[tx] or ''}\t{row}\t{median[tx]:.4f}\t{int(tx in keep)}\n")

    write_subset(lines, tx_gene, keep, f"{args.prefix}.high_confidence.gtf")

    category = {}
    with open_text(args.classification) as fh, open(f"{args.prefix}.high_confidence.classification.tsv", 'w') as out:
        reader = csv.DictReader(fh, delimiter='\t')
        for col in ('isoform_id', 'structural_category'):
            if col not in (reader.fieldnames or []):
                sys.exit(f"{args.classification}: no {col} column")
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in reader:
            category[row['isoform_id']] = row['structural_category']
            if row['isoform_id'] in keep:
                writer.writerow(row)
    unclassified = [tx for tx in catalogue if tx not in category]
    if unclassified:
        sys.exit(f"{args.classification}: {len(unclassified)} catalogue transcript(s) unclassified, e.g. {unclassified[:3]}")

    before = Counter(category[tx] for tx in catalogue)
    after = Counter(category[tx] for tx in keep)
    with open(f"{args.prefix}.summary.tsv", 'w') as out:
        out.write('metric\tcatalogue\thigh_confidence\n')
        out.write(f"transcripts\t{len(catalogue)}\t{len(keep)}\n")
        for cat in sorted(before):
            out.write(f"category:{cat}\t{before[cat]}\t{after.get(cat, 0)}\n")
        for s in names:
            out.write(f"assigned_reads:{s}\t{totals[s]:.10g}\t\n")

    print(f"{args.gtf}: kept {len(keep)} of {len(catalogue)} transcripts with median CPM > {args.min:g} "
          f"across {len(names)} sample(s)")


if __name__ == '__main__':
    main()
