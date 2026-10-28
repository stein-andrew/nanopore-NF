#!/usr/bin/env bash
# modkit pileup, written to the node's local disk, validated, then copied into the work dir
# and validated again.
#
#   validated_pileup.sh <bam> <out.bed> <threads> <max_depth> <flags file> <log file> [modkit args...]
#
# Valid means every line has bedMethyl's 18 columns and the line count equals the rows
# modkit reports writing. A pileup written straight to /scratch once lost ~11M rows in ~12k
# chunks, leaving cut-off lines, while modkit reported success (2026-10-27). The count goes to
# <out.bed>.rows so the next step can check the file it receives (bin/check_rows.sh).
set -euo pipefail
bam=$1 out=$2 threads=$3 depth=$4 flags=$5 log=$6; shift 6
die() { echo "validated_pileup.sh: $*" >&2; exit 1; }
tmp=$(mktemp -d "${TMPDIR:-/tmp}/pileup.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
echo "validated_pileup.sh: writing to $tmp ($(df -h --output=avail "$tmp" | tail -n 1 | tr -d ' ') free)" >&2

# shellcheck disable=SC2046
modkit pileup "$bam" "$tmp/out.bed" --threads "$threads" --max-depth "$depth" $(cat "$flags") \
    --log-filepath "$log" --suppress-progress "$@" 2> "$tmp/stderr" || { cat "$tmp/stderr" >&2; die "modkit pileup failed"; }
cat "$tmp/stderr" >&2
reported=$(grep -o 'processed [0-9]* rows' "$tmp/stderr" | grep -o '[0-9]*' | tail -n 1 || true)
[[ -n $reported ]] || die "modkit did not report how many rows it wrote, so the output cannot be checked"

count() { awk -F'\t' 'NF != 18 { b++ } END { print NR, b + 0 }' "$1"; }
read -r n bad < <(count "$tmp/out.bed")
[[ $n == "$reported" && $bad == 0 ]] || die "local output: $n lines ($bad not 18 columns), modkit reports $reported rows"
cp "$tmp/out.bed" "$out.partial"
read -r n2 bad2 < <(count "$out.partial")
[[ $n2 == "$n" && $bad2 == 0 ]] || die "copy in the work dir: $n2 lines ($bad2 not 18 columns), expected $n"
mv "$out.partial" "$out"
echo "$n" > "$out.rows"
echo "validated_pileup.sh: $out: $n rows, as modkit reported, all 18 columns" >&2
