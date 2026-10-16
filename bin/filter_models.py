#!/usr/bin/env python3
"""Keep transcript models above an abundance threshold (SALRR's per-sample CPM > 1 filter).

IsoQuant models, by CPM from IsoQuant's own model counts:
    filter_models.py --gtf S.transcript_models.gtf --counts S.discovered_transcript_counts.tsv \
        --min 1 --out S.isoquant.filtered.gtf --report S.isoquant.filtered.tsv

StringTie models, by an attribute on the transcript lines:
    filter_models.py --gtf S.stringtie.gtf --attr TPM \
        --min 1 --out S.stringtie.filtered.gtf --report S.stringtie.filtered.tsv

--counts: CPM = count / total assigned reads x 1e6, where "total assigned" is the sum over real
    features. IsoQuant's __ambiguous / __no_feature / __not_aligned summary rows are
    excluded: counting __no_feature as assigned would deflate every CPM by however many
    reads hit no model (20x in Connor's yeast run).
--attr: StringTie reports no read counts, so its per-transcript TPM stands in for CPM. SALRR
    states CPM > 1 for both tools without saying how StringTie's was derived.

A transcript is kept when its value is strictly greater than --min. Lines carrying a
transcript_id are kept or dropped with that transcript; gene lines are kept when any of
their transcripts is; comment lines pass through. Exits non-zero, rather than guessing, if
a model has no count, a transcript line lacks the attribute, or nothing was assigned.

read_gtf, read_counts, cpm and write_subset are shared with cohort_cpm.py, so the
per-sample and cohort filters cannot drift apart on what CPM means.
"""
import argparse
import gzip
import re
import sys

TX_ID = re.compile(r'transcript_id "([^"]+)"')
GENE_ID = re.compile(r'gene_id "([^"]+)"')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def attribute(attrs, key):
    """Value of `key "value"` in a GTF attribute column, or None."""
    m = re.search(r'(?:^|;)\s*' + re.escape(key) + r' "([^"]*)"', attrs)
    return m.group(1) if m else None


def read_gtf(path):
    """All lines, plus transcript_id -> gene_id for every transcript in the file."""
    with open_text(path) as fh:
        lines = fh.read().splitlines()
    tx_gene = {}
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 9:
            sys.exit(f"{path}: not a 9-column GTF line: {line[:80]}")
        tx = TX_ID.search(cols[8])
        if tx:
            gene = GENE_ID.search(cols[8])
            tx_gene[tx.group(1)] = gene.group(1) if gene else None
    return lines, tx_gene


def read_counts(path):
    """feature -> count from an IsoQuant feature_id/count TSV, summary (__) rows excluded."""
    counts = {}
    with open_text(path) as fh:
        fh.readline()                                   # feature_id <tab> count
        for line in fh:
            if not line.strip():
                continue
            feature, count = line.rstrip('\n').split('\t')[:2]
            if feature.startswith('__'):
                continue
            counts[feature] = float(count)
    return counts


def cpm(counts, path):
    """feature -> CPM over the total assigned reads, and that total."""
    total = sum(counts.values())
    if total <= 0:
        sys.exit(f"{path}: no assigned reads, so CPM is undefined")
    return {tx: c / total * 1e6 for tx, c in counts.items()}, total


def write_subset(lines, tx_gene, keep, out_path):
    """Write the lines of kept transcripts, the gene lines of their genes, and comments."""
    keep_genes = {tx_gene[tx] for tx in keep}
    with open(out_path, 'w') as out:
        for line in lines:
            if line.startswith('#') or not line.strip():
                out.write(line + '\n')
                continue
            cols = line.split('\t')
            tx = TX_ID.search(cols[8])
            if tx:
                if tx.group(1) in keep:
                    out.write(line + '\n')
            else:
                gene = GENE_ID.search(cols[8])
                if gene and gene.group(1) in keep_genes:
                    out.write(line + '\n')


def values_from_attr(lines, key, path):
    values = {}
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        if cols[2] != 'transcript':
            continue
        tx = TX_ID.search(cols[8]).group(1)
        v = attribute(cols[8], key)
        if v is None:
            sys.exit(f"{path}: transcript {tx} has no {key} attribute")
        values[tx] = float(v)
    return values


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', required=True, help='transcript models GTF')
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument('--counts', help='IsoQuant feature_id/count TSV for the models: filter on CPM')
    source.add_argument('--attr', help='transcript-line attribute to filter on, e.g. TPM for StringTie')
    ap.add_argument('--min', type=float, default=1.0, help='keep transcripts with value > MIN (default 1)')
    ap.add_argument('--out', required=True, help='filtered GTF')
    ap.add_argument('--report', required=True, help='per-transcript TSV: transcript_id, gene_id, value, kept')
    args = ap.parse_args()

    lines, tx_gene = read_gtf(args.gtf)
    if args.counts:
        values, total = cpm(read_counts(args.counts), args.counts)
        metric = 'CPM'
    else:
        values, total = values_from_attr(lines, args.attr, args.gtf), None
        metric = args.attr

    no_value = sorted(set(tx_gene) - set(values))
    if no_value:
        sys.exit(f"{len(no_value)} transcript(s) in {args.gtf} have no {metric}, e.g. {no_value[:3]}")

    keep = {tx for tx in tx_gene if values[tx] > args.min}
    write_subset(lines, tx_gene, keep, args.out)

    with open(args.report, 'w') as rep:
        rep.write(f"transcript_id\tgene_id\t{metric}\tkept\n")
        for tx in sorted(tx_gene, key=lambda t: -values[t]):
            rep.write(f"{tx}\t{tx_gene[tx] or ''}\t{values[tx]:.4f}\t{int(tx in keep)}\n")

    basis = f", total assigned reads {total:g}" if total is not None else ''
    print(f"{args.gtf}: kept {len(keep)} of {len(tx_gene)} transcripts with {metric} > {args.min:g}{basis}")


if __name__ == '__main__':
    main()
