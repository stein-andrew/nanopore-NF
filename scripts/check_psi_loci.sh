#!/usr/bin/env bash
# check_psi_loci.sh: where do the cohort's near-100% Psi sites come from?
#
#   1. clusters  Psi sites at >= 90% in >= MIN_SAMPLES samples, merged when
#                closer than GAP bp; clusters with >= MIN_SITES sites kept
#   2. share     for each modification, the fraction of all >= 90% calls
#                (sample x site) that fall inside those clusters
#   3. genes     GENCODE genes overlapping each cluster
#   4. reads     one sample's primary reads in each cluster (MAPQ, identity,
#                poly(A), length), with ACTB as a normal-gene baseline
#   5. rDNA      how much of each cluster's sequence aligns to RDNA
#                (chrUn_GL000220v1 carries a 45S unit; clusters on it are
#                the positive control)
#
# Usage: scripts/check_psi_loci.sh [SAMPLE]      (default ST1919)
# Needs samtools, minimap2 and gawk on PATH. Writes to OUT
# (default: the cohort folder's psi_loci_check/, so nothing lands in the repo).
# The >= 90 test is numeric-only: a plain awk ">= 90" passes "NA" as a string.
set -euo pipefail

SAMPLE=${1:-ST1919}
R=${R:-/scratch/stein.an/ont-pipeline/cohort-v1/results}
MATRIX=${MATRIX:-$R/cohort/GRCh38/cohort.GRCh38.mods.pct_fpcorr.tsv.gz}
SHEET=${SHEET:-/scratch/stein.an/ont-pipeline/cohort-v1/analysis_samplesheet.csv}
GTF=${GTF:-/projects/Genometechlab/referenceFastaFiles/human/gencode/gencode.v47.annotation.gtf}
FA=${FA:-/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa}
RDNA=${RDNA:-chrUn_GL000220v1}
MIN_SAMPLES=${MIN_SAMPLES:-10}
GAP=${GAP:-10000}
MIN_SITES=${MIN_SITES:-10}
MAX_READS=${MAX_READS:-20000}
OUT=${OUT:-$(dirname "$R")/psi_loci_check}

for t in samtools minimap2 gawk; do
    command -v "$t" > /dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done

if [[ -z ${BAM:-} ]]; then
    BAM=$(gawk -F',' -v s="$SAMPLE" '
        { gsub(/\r|"/, "") }
        NR == 1 { for (i = 1; i <= NF; i++) { h = tolower($i); if (h == "sample_id" || (h == "sample" && !cs)) cs = i; if (h == "bam") cb = i }; next }
        cs && cb && $cs == s { print $cb; exit }' "$SHEET")
fi
[[ -n $BAM && -s $BAM ]] || { echo "ERROR: no BAM for $SAMPLE in $SHEET (header: $(head -n 1 "$SHEET")); set BAM=" >&2; exit 1; }
mkdir -p "$OUT"
exec > >(tee "$OUT/report.txt")

echo "== 1. clusters: Psi >= 90% in >= $MIN_SAMPLES samples, merged within $GAP bp"
zcat "$MATRIX" \
  | gawk -F'\t' -v min="$MIN_SAMPLES" '
        NR > 1 && $4 == "17802" {
            n = 0
            for (i = 6; i <= NF; i++) if ($i ~ /^[0-9.]+$/ && $i + 0 >= 90) n++
            if (n >= min) print $1 "\t" $2
        }' \
  | sort -k1,1 -k2,2n \
  | gawk -F'\t' -v OFS='\t' -v gap="$GAP" -v minsites="$MIN_SITES" '
        function flush() { if (n >= minsites) print c, s, e + 1, n }
        NR == 1 || $1 != c || $2 - e > gap { if (NR > 1) flush(); c = $1; s = $2; n = 0 }
        { e = $2; n++ }
        END { if (NR) flush() }' > "$OUT/clusters.bed"
[[ -s $OUT/clusters.bed ]] || { echo "no clusters found"; exit 0; }
gawk -F'\t' 'BEGIN { printf "%-4s %-26s %11s %11s %8s %6s\n", "id", "chrom", "start", "end", "span_kb", "sites" }
    { printf "%-4s %-26s %11d %11d %8.1f %6d\n", "c" NR, $1, $2, $3, ($3 - $2) / 1000, $4 }' "$OUT/clusters.bed"

echo
echo "== 2. share of all >= 90% calls (sample x site) inside the clusters"
zcat "$MATRIX" | gawk -F'\t' '
    BEGIN {
        nk = split("a m 17802 17596 69426 19228 19227 19229", code, " ")
        split("m6A m5C Psi Ino Am Cm Um Gm", label, " ")
    }
    ARGIND == 1 { nc++; cc[nc] = $1; cs[nc] = $2; ce[nc] = $3; has[$1] = 1; next }
    FNR == 1 { next }
    {
        k = 0
        for (i = 6; i <= NF; i++) if ($i ~ /^[0-9.]+$/ && $i + 0 >= 90) k++
        if (!k) next
        tot[$4] += k
        if ($1 in has)
            for (j = 1; j <= nc; j++)
                if ($1 == cc[j] && $2 >= cs[j] && $2 < ce[j]) { inc[$4] += k; break }
    }
    END {
        printf "%-5s %12s %12s %7s\n", "mod", "calls_ge90", "in_clusters", "pct"
        for (j = 1; j <= nk; j++) if (code[j] in tot)
            printf "%-5s %12d %12d %6.1f%%\n", label[j], tot[code[j]], inc[code[j]], 100 * inc[code[j]] / tot[code[j]]
    }' "$OUT/clusters.bed" -

echo
echo "== 3. GENCODE genes overlapping each cluster"
gawk -F'\t' -v OFS='\t' -v actb="$OUT/actb.region" '
    ARGIND == 1 { nc++; cc[nc] = $1; cs[nc] = $2; ce[nc] = $3; has[$1] = 1; next }
    $3 != "gene" { next }
    {
        match($9, /gene_name "[^"]+"/); gn = substr($9, RSTART + 11, RLENGTH - 12)
        if (gn == "ACTB") print $1 ":" $4 "-" $5 > actb
        if (!($1 in has)) next
        match($9, /gene_type "[^"]+"/); gt = substr($9, RSTART + 11, RLENGTH - 12)
        for (j = 1; j <= nc; j++)
            if ($1 == cc[j] && $4 <= ce[j] && $5 > cs[j]) { print "c" j, gn, gt, $7; hit[j] = 1 }
    }
    END { for (j = 1; j <= nc; j++) if (!(j in hit)) print "c" j, "(no annotated gene)" }' \
    "$OUT/clusters.bed" "$GTF" | sort -k1,1V | column -t -s $'\t'

echo
echo "== 4. $SAMPLE primary reads, first $MAX_READS per region ($(basename "$BAM"))"
read_stats() {
    samtools view -F 0x904 "$BAM" "$2" | gawk -F'\t' -v label="$1" -v max="$MAX_READS" '
        function med(a, n,   b, k) {
            if (!n) return "NA"
            k = asort(a, b)
            return sprintf("%.1f", (k % 2) ? b[(k + 1) / 2] : (b[k / 2] + b[k / 2 + 1]) / 2)
        }
        NR > max { next }
        {
            n++; mq[n] = $5; if ($5 == 0) zero++
            if ($10 != "*") { nl++; rl[nl] = length($10) }
            M = I = D = 0; cig = $6
            while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
                len = substr(cig, RSTART, RLENGTH - 1) + 0; op = substr(cig, RSTART + RLENGTH - 1, 1)
                if (op ~ /[M=X]/) M += len; else if (op == "I") I += len; else if (op == "D") D += len
                cig = substr(cig, RSTART + RLENGTH)
            }
            for (f = 12; f <= NF; f++) {
                if ($f ~ /^NM:i:/ && M + I + D) { ni++; id[ni] = 100 * (M - (substr($f, 6) - I - D)) / (M + I + D) }
                else if ($f ~ /^pt:i:/) { np++; pt[np] = substr($f, 6) + 0 }
            }
        }
        END {
            if (!n) { printf "%-36s %7d\n", label, 0; exit }
            printf "%-36s %7d %8s %7.1f %10s %8s %8s\n", label, n, med(mq, n), 100 * zero / n, med(id, ni), med(pt, np), med(rl, nl)
        }'
}
printf "%-36s %7s %8s %7s %10s %8s %8s\n" region reads mapq_md "%mapq0" "ident_md%" pA_md len_md
[[ -s $OUT/actb.region ]] && read_stats "ACTB (baseline)" "$(cat "$OUT/actb.region")"
k=0
while IFS=$'\t' read -r c s e n; do
    ((++k))
    read_stats "c$k $c:$((s / 1000)) kb" "$c:$((s + 1))-$e"
done < "$OUT/clusters.bed"

echo
echo "== 5. cluster sequence aligned to $RDNA (minimap2 map-ont)"
samtools faidx "$FA" "$RDNA" > "$OUT/rdna.fa"
k=0
while IFS=$'\t' read -r c s e n; do
    ((++k))
    samtools faidx "$FA" "$c:$((s + 1))-$e" | sed "1s/.*/>c$k/"
done < "$OUT/clusters.bed" > "$OUT/clusters.fa"
minimap2 -c -x map-ont --secondary=no "$OUT/rdna.fa" "$OUT/clusters.fa" \
    2> "$OUT/minimap2.log" > "$OUT/clusters_vs_rdna.paf"
gawk -F'\t' '
    ARGIND == 1 { al[$1] += $4 - $3; m[$1] += $10; b[$1] += $11; next }
    {
        q = "c" FNR; len = $3 - $2
        if (q in al) printf "%-4s %-26s covered %5.1f%%   identity %5.1f%%\n", q, $1, 100 * (al[q] > len ? 1 : al[q] / len), 100 * m[q] / b[q]
        else printf "%-4s %-26s no alignment\n", q, $1
