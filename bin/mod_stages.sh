#!/usr/bin/env bash
# Filter stages of a modkit pileup bedMethyl, before and after gIVT false-positive correction.
#
#   mod_stages.sh failcanon <raw.bed> <out.bed>
#       the pileup with calls that failed their threshold counted as canonical (the gIVT tables'
#       definition): Nvalid_cov (col 10, and the score, col 5) += Nfail, Ncanonical (col 13)
#       += Nfail, Nfail (col 16) = 0, percent modified (col 11) = 100 Nmod / Nvalid_cov with
#       two decimals, as modkit writes it; every other column as it was
#   mod_stages.sh coverage <raw.bed> <min_cov> <min_pct> <prefix>
#       <prefix>.cov<C>.bed           Nvalid_cov (col 10) >= C
#       <prefix>.cov<C>_pct<P>.bed    ... and percent modified (col 11) >= P
#   mod_stages.sh restage <fpcorr.bed> <min_cov> <min_pct> <prefix>
#       the same two stages again, from a 19-column file of fp_correct.py --filter_mismatch:
#       <prefix>.cov<C>.bed is its columns 1-18 (each kept input line as it was), so the
#       uncorrected stages hold exactly the rows the corrected ones were made from
#   mod_stages.sh fpcorr <fpcorr.bed> <min_pct> <out.bed>
#       rows whose FP-corrected percent (col 19, from IVT_fp_correction.py) is a number >= P.
#       gIVT leaves col 19 empty where it cannot correct (base mismatch, no 9mer, kmer not in
#       its table); those can't pass a corrected filter, so they are dropped here.
#   mod_stages.sh summary <out.tsv> <stage>=<bed> [<stage>=<bed> ...]
#       rows per stage per mod code (col 4); for 19-column files also how many are correctable
#
# Checks column counts (18 for modkit's bedMethyl, 19 after gIVT, 21 with fp_correct.py
# --wilson_columns) instead of trusting them: a modkit release that changes the layout fails
# here, not as silently wrong filters.
set -euo pipefail

die() { echo "mod_stages.sh: $*" >&2; exit 1; }

# Column count of the first non-empty line; 0 for an empty file.
ncol() { awk -F'\t' 'NF { print NF; exit } END { if (NR == 0) print 0 }' "$1" | head -n1; }

expect_cols() {     # file n [n ...]: any of the given counts
    local f=$1 n w; shift; n=$(ncol "$f")
    [[ $n == 0 ]] && return
    for w; do [[ $n == "$w" ]] && return; done
    die "$f has $n tab-separated columns, expected $*"
}

number() { [[ $2 =~ ^[0-9]+(\.[0-9]+)?$ ]] || die "$1 must be a number, got '$2'"; }

case ${1:-} in
failcanon)
    (( $# == 3 )) || die "usage: $0 failcanon <raw.bed> <out.bed>"
    in=$2 out=$3
    expect_cols "$in" 18
    awk -F'\t' -v OFS='\t' '{ f = $16 + 0; $5 += f; $10 += f; $13 += f; $16 = 0
                              $11 = sprintf("%.2f", $10 > 0 ? 100 * $12 / $10 : 0); print }' "$in" > "$out"
    ;;
coverage)
    (( $# == 5 )) || die "usage: $0 coverage <raw.bed> <min_cov> <min_pct> <prefix>"
    raw=$2 cov=$3 pct=$4 prefix=$5
    number min_cov "$cov"; number min_pct "$pct"
    expect_cols "$raw" 18
    awk -F'\t' -v c="$cov" '$10 + 0 >= c' "$raw" > "$prefix.cov$cov.bed"
    awk -F'\t' -v p="$pct" '$11 + 0 >= p' "$prefix.cov$cov.bed" > "$prefix.cov${cov}_pct$pct.bed"
    ;;
restage)
    (( $# == 5 )) || die "usage: $0 restage <fpcorr.bed> <min_cov> <min_pct> <prefix>"
    in=$2 cov=$3 pct=$4 prefix=$5
    number min_cov "$cov"; number min_pct "$pct"
    expect_cols "$in" 19 21
    cut -f1-18 "$in" > "$prefix.cov$cov.bed"
    awk -F'\t' -v p="$pct" '$11 + 0 >= p' "$prefix.cov$cov.bed" > "$prefix.cov${cov}_pct$pct.bed"
    ;;
fpcorr)
    (( $# == 4 )) || die "usage: $0 fpcorr <fpcorr.bed> <min_pct> <out.bed>"
    in=$2 pct=$3 out=$4
    number min_pct "$pct"
    expect_cols "$in" 19 21
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
    die "usage: $0 {failcanon|coverage|restage|fpcorr|summary} ..."
    ;;
esac
