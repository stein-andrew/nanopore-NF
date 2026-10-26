#!/usr/bin/env bash
# Add a masked column (the mask region's 4th column, IG / TR / HLA, or .) to a calls table
# whose first two columns are CHROM and 1-based POS.
#
#   mask_calls.sh <mask.bed> <calls.tsv> > calls.masked.tsv
#
# A call is masked when its position lies inside a mask interval (BED: start < POS <= end).
set -euo pipefail
awk -F'\t' -v OFS='\t' '
    NR == FNR { n++; c[n] = $1; s[n] = $2; e[n] = $3; t[n] = $4; on[$1] = on[$1] " " n; next }
    { m = "."; k = split(on[$1], idx, " ")
      for (j = 1; j <= k; j++) { i = idx[j]; if ($2 > s[i] && $2 <= e[i]) { m = t[i]; break } }
      print $0, m }
' "$1" "$2"
