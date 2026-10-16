#!/usr/bin/env python3
"""Make isomatch's catalogue GTF load in IsoQuant with --complete_genedb, by giving it the
shape GENCODE has (the shape that loads correctly):

    gtf_add_genes.py cohort.merged.gtf.gz > cohort.catalogue.gtf

  * a gene record for every gene_id, written directly before that gene's first transcript,
    spanning its transcripts (same chromosome and strand, min start to max end);
  * isomatch's ##ISOM header lines demoted to plain # comments.

Why: given isomatch's raw layout (##ISOM directive lines carrying key="value" pairs, and in
the first version of this script every gene line up front), IsoQuant's gffutils loader did
not take transcript_id/gene_id as IDs: every transcript came out as transcript_N and every
gene as having no transcripts, so all reads were intergenic. gffutils guesses the format
from the first lines it sees; these now look like GENCODE's.

A gene whose transcripts sit on different chromosomes or strands cannot be one gene
record, so that is an error rather than a guess; so is a GTF that already has gene lines.
"""
import sys

from filter_models import GENE_ID, open_text


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <transcripts+exons GTF[.gz]> > with_genes.gtf")
    path = sys.argv[1]

    with open_text(path) as fh:
        lines = fh.read().splitlines()

    genes = {}                                    # gene_id -> [chrom, source, strand, start, end]
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 9:
            sys.exit(f"{path}: not a 9-column GTF line: {line[:80]}")
        if cols[2] == 'gene':
            sys.exit(f"{path} already has gene lines; nothing to add")
        if cols[2] != 'transcript':
            continue
        gene = GENE_ID.search(cols[8])
        if not gene:
            sys.exit(f"{path}: transcript line without gene_id: {line[:120]}")
        gid = gene.group(1)
        chrom, source, start, end, strand = cols[0], cols[1], int(cols[3]), int(cols[4]), cols[6]
        if gid not in genes:
            genes[gid] = [chrom, source, strand, start, end]
            continue
        g = genes[gid]
        if (g[0], g[2]) != (chrom, strand):
            sys.exit(f"{path}: gene {gid} has transcripts on {g[0]}{g[2]} and {chrom}{strand}")
        g[3], g[4] = min(g[3], start), max(g[4], end)

    if not genes:
        sys.exit(f"{path}: no transcript lines")

    out = sys.stdout
    written = set()
    for line in lines:
        if not line.strip():
            continue
        if line.startswith('#'):
            out.write('# ' + line.lstrip('#').lstrip() + '\n')       # ##ISOM ... -> # ISOM ...
            continue
        cols = line.split('\t')
        if cols[2] == 'transcript':
            gid = GENE_ID.search(cols[8]).group(1)
            if gid not in written:
                chrom, source, strand, start, end = genes[gid]
                out.write(f'{chrom}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}";\n')
                written.add(gid)
        out.write(line + '\n')
    print(f"{path}: added {len(genes)} gene lines", file=sys.stderr)


if __name__ == '__main__':
    main()
