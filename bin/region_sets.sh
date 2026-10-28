#!/usr/bin/env bash
# Region sets from the annotation GTF, once per reference.
#
#   region_sets.sh <annotation.gtf> <globin gene names, comma-separated> <MHC region or ''> <prefix>
#
# <prefix>.globin.bed   spans of the named globin genes (gene_name), for the per-sample
#                       globin read fraction
# <prefix>.mask.bed     variant mask: every GENCODE gene whose gene_type is IG_* or TR_*
#                       (antibody and T-cell receptor segments: rearranged and hypermutated
#                       within each cell, so pooled blood RNA shows somatic "variants" there),
#                       plus the MHC region as HLA (germline, but too polymorphic for
#                       reference alignment to give clean calls). 4th column: IG / TR / HLA.
# <prefix>.region_sets.log  what was found, so a renamed or missing gene shows up.
set -euo pipefail
gtf=$1 genes=$2 mhc=$3 prefix=$4

awk -F'\t' -v OFS='\t' -v genes="$genes" -v gb="$prefix.globin.tmp" -v mk="$prefix.igtr.tmp" '
    BEGIN { n = split(genes, a, ","); for (i = 1; i <= n; i++) want[a[i]] = 1; printf "" > gb; printf "" > mk }
    $3 != "gene" { next }
    match($9, /gene_name "[^"]+"/) { name = substr($9, RSTART + 11, RLENGTH - 12); if (name in want) print $1, $4 - 1, $5, name > gb }
    match($9, /gene_type "(IG|TR)_[^"]*"/) { print $1, $4 - 1, $5, substr($9, RSTART + 11, 2) > mk }
' "$gtf"

sort -k1,1 -k2,2n "$prefix.globin.tmp" > "$prefix.globin.bed"
{
    cat "$prefix.igtr.tmp"
    if [[ -n $mhc ]]; then
        [[ $mhc =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]] || { echo "region_sets.sh: MHC region must be chrom:start-end, got '$mhc'" >&2; exit 1; }
        printf '%s\t%d\t%d\tHLA\n' "${BASH_REMATCH[1]}" $(( BASH_REMATCH[2] - 1 )) "${BASH_REMATCH[3]}"
    fi
} | sort -k1,1 -k2,2n > "$prefix.mask.bed"
rm -f "$prefix.globin.tmp" "$prefix.igtr.tmp"

found=$(cut -f4 "$prefix.globin.bed" | sort -u | paste -sd, -)
missing=$(tr ',' '\n' <<< "$genes" | grep -vxF -f <(cut -f4 "$prefix.globin.bed" | sort -u) | paste -sd, - || true)
{
    echo "globin genes found: ${found:-none}"
    echo "globin genes not in the annotation: ${missing:-none}"
    echo "mask regions: $(cut -f4 "$prefix.mask.bed" | sort | uniq -c | awk '{printf "%s %s  ", $2, $1}')"
} > "$prefix.region_sets.log"
# No globin gene in the annotation (a single-chromosome test annotation, another species) is
# not an error: the read QC then reports 0 globin reads, and the log says why.
if [[ ! -s $prefix.globin.bed ]]; then
    echo "WARNING: none of the globin genes ($genes) are in $gtf: the read QC will report 0 globin reads" \
        | tee -a "$prefix.region_sets.log" >&2
fi
