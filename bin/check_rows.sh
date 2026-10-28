#!/usr/bin/env bash
# Fail unless a table arrived whole: exactly the line count recorded beside it when it was
# written, every line with the expected number of tab-separated columns.
#
#   check_rows.sh <file> <file.rows> <columns>
set -euo pipefail
f=$1 rows=$2 ncol=$3
want=$(tr -d '[:space:]' < "$rows")
[[ $want =~ ^[0-9]+$ ]] || { echo "check_rows.sh: $rows does not hold a line count ('$want')" >&2; exit 1; }
read -r n bad < <(awk -F'\t' -v c="$ncol" 'NF != c { b++ } END { print NR, b + 0 }' "$f")
if [[ $n != "$want" || $bad != 0 ]]; then
    echo "check_rows.sh: $f has $n lines ($bad without $ncol columns), but $want were written: damaged on the way" >&2
    exit 1
fi
echo "check_rows.sh: $f: $n lines, $ncol columns each, as written" >&2
