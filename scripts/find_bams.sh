#!/usr/bin/env bash
# Inventory of basecalled BAMs under one or more folders, a file list for rsync, and an
# RNA_analysis.nf samplesheet.
#
#   scripts/find_bams.sh <out prefix> <folder> [folder...]
#
# Finds every *.bam of at least MIN_SIZE that has an index beside it (.bam.bai, .bai or .csi:
# final sorted BAMs; Nextflow intermediates usually have none). A file reachable by several
# paths (outputs are published as hard links, so a BAM and its work-dir original are one
# file) is listed once, by a path outside any work/ directory if it has one. Each header gives:
# sample (@RG SM), basecall and modification models (@RG DS), dorado version (@PG), a
# fingerprint of the reference (every @SQ name and length), and mapped reads (from the index).
#
#   <prefix>.inventory.tsv    one row per BAM, all of the above
#   <prefix>.paths.txt        the RNA BAMs and their indexes, for rsync --files-from
#   <prefix>.samplesheet.csv  the RNA BAMs, for RNA_analysis.nf
#
# and prints what must hold before they are analysed together: one reference, one basecall
# model and dorado version (gIVT's error table is specific to dorado 1.0 calls, and mixed
# basecallers confound between-sample differences), and no sample twice.
#
# Environment: SAMTOOLS (default "samtools"; e.g. "apptainer exec --bind /scratch img.sif samtools"),
# MIN_SIZE (default 1G), REFERENCE_NAME (samplesheet reference, default GRCh38),
# SAMPLE_REGEX (sample from the file name when @RG SM is missing, default ST[0-9]+).
set -euo pipefail
die() { echo "find_bams.sh: $*" >&2; exit 1; }
(( $# >= 2 )) || die "usage: $0 <out prefix> <folder> [folder...]"
prefix=$1; shift
ST=${SAMTOOLS:-samtools}
MIN=${MIN_SIZE:-1G}
REF=${REFERENCE_NAME:-GRCh38}
RX=${SAMPLE_REGEX:-ST[0-9]+}
# shellcheck disable=SC2086
$ST --version > /dev/null 2>&1 || die "samtools not runnable as '$ST' (set SAMTOOLS)"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

# 1. candidates: fd if installed (parallel; -u also looks in hidden and git-ignored folders)
if command -v fd > /dev/null; then
    echo "searching with fd: $*" >&2
    fd --unrestricted --type f --extension bam --size "+${MIN,,}" --absolute-path . "$@" 2> "$tmp/errors" > "$tmp/candidates" || true
else
    echo "fd not found; searching with find: $*" >&2
    find "$@" -type f -name '*.bam' -size "+$MIN" 2> "$tmp/errors" > "$tmp/candidates" || true
fi
[[ -s $tmp/errors ]] && echo "   ($(wc -l < "$tmp/errors") folders could not be read; first: $(head -n 1 "$tmp/errors"))" >&2
echo "   $(wc -l < "$tmp/candidates") BAM files of at least $MIN" >&2

# 2. one path per file (device:inode), preferring a path outside work/
tr '\n' '\0' < "$tmp/candidates" | xargs -0 -r stat --printf '%d:%i\t%n\n' \
    | awk -F'\t' '{ k = $1; p = $2
                    if (!(k in best) || (best[k] ~ /\/work\// && p !~ /\/work\//)) best[k] = p }
                  END { for (k in best) print best[k] }' | sort > "$tmp/unique"
echo "   $(wc -l < "$tmp/unique") distinct files" >&2

# 3. indexed ones, with their headers
printf 'path\tindex\tsize_gb\tmodified\tsample\tdata_type\tbasecall_model\tmodbase_models\tdorado\tref_fingerprint\tn_contigs\tmapped_reads\n' > "$prefix.inventory.tsv"
noindex=0
while IFS= read -r p; do
    idx=""
    for c in "$p.bai" "${p%.bam}.bai" "$p.csi"; do [[ -e $c ]] && { idx=$c; break; }; done
    if [[ -z $idx ]]; then noindex=$((noindex + 1)); continue; fi
    # shellcheck disable=SC2086
    h=$($ST view -H "$p" 2> /dev/null) || h=""
    field() { grep "^@RG" <<< "$h" | grep -o "$1[^[:space:]]*" | cut -d"$2" -f2- | sort -u | paste -sd'|' -; }
    sm=$(field 'SM:' ':'); model=$(field 'basecall_model=' '='); mods=$(field 'modbase_models=' '=')
    dorado=$(grep '^@PG' <<< "$h" | grep -E 'ID:basecaller|PN:dorado' | grep -o 'VN:[^[:space:]]*' | cut -d: -f2 | sort -u | paste -sd'|' -)
    fp=$(grep '^@SQ' <<< "$h" | grep -o 'SN:[^[:space:]]*\|LN:[^[:space:]]*' | paste -d' ' - - | md5sum | cut -c1-10)
    nsq=$(grep -c '^@SQ' <<< "$h" || true)
    # shellcheck disable=SC2086
    mapped=$($ST idxstats "$p" 2> /dev/null | awk '{ s += $3 } END { print s + 0 }')
    if [[ -z $sm || $sm == *'|'* ]]; then
        sm=$(basename "$p" | grep -oE "$RX" | head -n 1 || true); sm=${sm:-$(basename "$p" .bam)}
    fi
    case ${model,,} in rna*) type=RNA ;; dna*) type=DNA ;; *) type=. ;; esac
    printf '%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$p" "$idx" "$(stat -c %s "$p" | awk '{ print $1 / 1e9 }')" \
        "$(stat -c %y "$p" | cut -d. -f1)" "$sm" "$type" "${model:-.}" "${mods:-.}" "${dorado:-.}" "$fp" "$nsq" "$mapped" >> "$prefix.inventory.tsv"
done < "$tmp/unique"

# 4. the RNA BAMs: rsync list and samplesheet
awk -F'\t' 'NR > 1 && $6 == "RNA" { print $1; print $2 }' "$prefix.inventory.tsv" > "$prefix.paths.txt"
{ echo "sample_id,reference,data_type,bam,bai,basecall_model,mod_models,dorado_version"
  awk -F'\t' -v OFS=, -v ref="$REF" 'NR > 1 && $6 == "RNA" { gsub(/,/, ";", $8); print $5, ref, "RNA", $1, $2, $7, $8, $9 }' "$prefix.inventory.tsv"
} > "$prefix.samplesheet.csv"

# 5. what must hold before these are analysed together
inv=$prefix.inventory.tsv
echo "== $(($(wc -l < "$inv") - 1)) indexed BAMs ($noindex without an index skipped): $(awk -F'\t' 'NR > 1 { n[$6]++ } END { for (t in n) printf "%s %d  ", t, n[t] }' "$inv")"
echo "   RNA total: $(awk -F'\t' 'NR > 1 && $6 == "RNA" { s += $3 } END { printf "%.0f GB", s }' "$inv")"
for col in 10:reference 7:basecall_model 9:dorado; do
    n=${col%%:*}; name=${col#*:}
    k=$(awk -F'\t' -v c="$n" 'NR > 1 && $6 == "RNA" { print $c }' "$inv" | sort -u | wc -l)
    if (( k > 1 )); then echo "   CHECK: $k different ${name}s among the RNA BAMs:"; awk -F'\t' -v c="$n" 'NR > 1 && $6 == "RNA" { print "      " $c }' "$inv" | sort | uniq -c
    else echo "   one $name among the RNA BAMs: $(awk -F'\t' -v c="$n" 'NR > 1 && $6 == "RNA" { print $c; exit }' "$inv")"; fi
done
dups=$(awk -F'\t' 'NR > 1 && $6 == "RNA" { print $5 }' "$inv" | sort | uniq -d)
if [[ -n $dups ]]; then echo "   CHECK: samples with more than one RNA BAM (keep one each in the samplesheet): $(paste -sd' ' <<< "$dups")"
else echo "   every RNA sample has one BAM"; fi
echo "written: $prefix.inventory.tsv  $prefix.paths.txt  $prefix.samplesheet.csv"
