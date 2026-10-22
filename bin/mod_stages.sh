#!/usr/bin/env bash
# Filter stages of a modkit pileup bedMethyl, before and after gIVT false-positive correction.
#
#   mod_stages.sh coverage <raw.bed> <min_cov> <min_pct> <prefix>
#       <prefix>.cov<C>.bed           Nvalid_cov (col 10) >= C
#       <prefix>.cov<C>_pct<P>.bed    ... and percent modified (col 11) >= P
#   mod_stages.sh fpcorr <fpcorr.bed> <min_pct> <out.bed>
#       rows whose FP-corrected percent (col 19, from IVT_fp_correction.py) is a number >= P.
#       gIVT leaves col 19 empty where it cannot correct (base mismatch, no 9mer, kmer not in
#       its table); those can't pass a corrected filter, so they are dropped here.
#   mod_stages.sh summary <out.tsv> <stage>=<bed> [<stage>=<bed> ...]
#       rows per stage per mod code (col 4); for 19-column files also how many are correctable
#
# Checks column counts (18 for modkit's bedMethyl, 19 after gIVT) instead of trusting them:
# a modkit release that changes the layout fails here, not as silently wrong filters.
set -euo pipefail

die() { echo "mod_stages.sh: $*" >&2; exit 1; }

# Column count of the first non-empty line; 0 for an empty file.
ncol() { awk -F'\t' 'NF { print NF; exit } END { if (NR == 0) print 0 }' "$1" | head -n1; }

expect_cols() {     # file n
    local n; n=$(ncol "$1")
    [[ $n == 0 || $n == "$2" ]] || die "$1 has $n tab-separated columns, expected $2"
}

number() { [[ $2 =~ ^[0-9]+(\.[0-9]+)?$ ]] || die "$1 must be a number, got '$2'"; }

case ${1:-} in
coverage)
    (( $# == 5 )) || die "usage: $0 coverage <raw.bed> <min_cov> <min_pct> <prefix>"
    raw=$2 cov=$3 pct=$4 prefix=$5
    number min_cov "$cov"; number min_pct "$pct"
    expect_cols "$raw" 18
    awk -F'\t' -v c="$cov" '$10 + 0 >= c' "$raw" > "$prefix.cov$cov.bed"
    awk -F'\t' -v p="$pct" '$11 + 0 >= p' "$prefix.cov$cov.bed" > "$prefix.cov${cov}_pct$pct.bed"
    ;;
fpcorr)
    (( $# == 4 )) || die "usage: $0 fpcorr <fpcorr.bed> <min_pct> <out.bed>"
    in=$2 pct=$3 out=$4
    number min_pct "$pct"
    expect_cols "$in" 19
    awk -F'\t' -v p="$pct" '$19 ~ /^[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/ && $19 + 0 >= p' "$in" > "$out"
    ;;
summary)
    (( $# >= 3 )) || die "usage: $0 summary <out.tsv> <stage>=<bed> [...]"
    out=$2; shift 2
    printf 'stage\tmod_code\trows\tcorrectable\n' > "$out"
    for item in "$@"; do
        stage=${item%%=*} bed=${item#*=}
        [[ $item == *=* && -n $stage && -f $bed ]] || die "summary: expected <stage>=<existing file>, got: $item"
        awk -F'\t' -v s="$stage" '
            NF { n[$4]++; if (NF >= 19 && $19 != "") c[$4]++ }
            END {
                for (m in n) printf "%s\t%s\t%d\t%s\n", s, m, n[m], (NF >= 19 ? c[m] + 0 : "NA")
                if (NR == 0) printf "%s\t(none)\t0\tNA\n", s
            }' "$bed" | sort -t$'\t' -k2,2 >> "$out"
    done
    ;;
*)
    die "usage: $0 {coverage|fpcorr|summary} ..."
    ;;
esac
