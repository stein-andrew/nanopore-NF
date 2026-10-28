#!/usr/bin/env bash
# Per-sample read composition: mapped reads, globin reads, mitochondrial reads.
#
#   read_qc.sh <bam> <globin.bed> <sample> <mito contig> > S.read_qc.tsv
#
# Whole-blood RNA can be mostly globin; every globin read is depth lost to everything
# else, and the fraction differs between samples, so it belongs with any between-sample
# comparison. Counts are primary alignments (-F 2308); a read overlapping several globin
# genes counts once.
set -euo pipefail
bam=$1 bed=$2 sample=$3 mito=$4
read -r mapped mt < <(samtools idxstats "$bam" | awk -v m="$mito" '{ s += $3 } $1 == m { mt = $3 } END { print s + 0, mt + 0 }')
# an annotation without globin genes gives an empty BED: 0, without asking samtools what an
# empty region list means
if [[ -s $bed ]]; then glob=$(samtools view -c -F 2308 -M -L "$bed" "$bam"); else glob=0; fi
awk -v OFS='\t' -v s="$sample" -v n="$mapped" -v g="$glob" -v m="$mt" 'BEGIN {
    print "sample", "mapped_reads", "globin_reads", "globin_pct", "mito_reads", "mito_pct", "non_globin_reads"
    printf "%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\n", s, n, g, n ? 100 * g / n : 0, m, n ? 100 * m / n : 0, n - g }'
