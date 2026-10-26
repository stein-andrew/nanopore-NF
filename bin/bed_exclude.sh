#!/usr/bin/env bash
# Drop BED rows (chrom, 0-based start, end in columns 1-3) that overlap any region in
# <regions.bed>; report how many were dropped on stderr.
#
#   bed_exclude.sh <regions.bed> <in.bed> > out.bed
set -euo pipefail
awk -F'\t' -v OFS='\t' '
    NR == FNR { n++; c[n] = $1; s[n] = $2; e[n] = $3; on[$1] = on[$1] " " n; next }
    { k = split(on[$1], idx, " ")
      for (j = 1; j <= k; j++) { i = idx[j]; if ($2 < e[i] && $3 > s[i]) { dropped++; next } }
      print }
    END { printf "bed_exclude.sh: dropped %d row(s) of %s inside %s\n", dropped + 0, FILENAME, ARGV[1] > "/dev/stderr" }
' "$1" "$2"
