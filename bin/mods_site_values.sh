#!/usr/bin/env bash
# The cohort's final-stage modification sites, and every sample's values at them.
#
#   mods_site_values.sh <min_pct> <sites.tsv> <site_values.tsv> SAMPLE=<S.cov20_fpcorr.bed> ...
#
# sites.tsv        chrom start end mod strand: every site whose gIVT-corrected percent
#                  (column 19) is >= min_pct in at least one sample; the input stage already
#                  holds only sites with coverage >= its minimum
# site_values.tsv  chrom start end mod strand sample n_valid n_mod pct pct_fpcorr, for each
#                  sample that has the site at all (column 10, 12, 11, 19); a sample without
#                  the site had too little coverage there
set -euo pipefail
min=$1 sites=$2 long=$3; shift 3
files=()
for pair in "$@"; do
    [[ $pair == *=* ]] || { echo "mods_site_values.sh: expected SAMPLE=PATH, got $pair" >&2; exit 1; }
    files+=("${pair#*=}")
done
awk -F'\t' -v OFS='\t' -v min="$min" '$19 != "" && $19 + 0 >= min + 0 { print $1, $2, $3, $4, $6 }' "${files[@]}" \
    | sort -u -k1,1 -k2,2n -k4,4 -k5,5 > "$sites"
: > "$long"
for pair in "$@"; do
    awk -F'\t' -v OFS='\t' -v s="${pair%%=*}" '
        NR == FNR { want[$1 FS $2 FS $3 FS $4 FS $5] = 1; next }
        ($1 FS $2 FS $3 FS $4 FS $6) in want { print $1, $2, $3, $4, $6, s, $10, $12, $11, $19 }
    ' "$sites" "${pair#*=}" >> "$long"
done
