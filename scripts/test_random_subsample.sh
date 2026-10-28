#!/usr/bin/env bash
#SBATCH --job-name=subsample-test
#SBATCH --partition=short
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=subsample-test_%j.out
#
# Does "the first N reads at a position" (modkit pileup --max-depth) give different site
# values from "a random N"? Builds a BAM in which each globin locus is randomly subsampled to
# at most CAP reads and every other read is kept, runs modkit pileup on it with the sample's
# own flags and the same --max-depth (which then never binds), and compares every site with
# the sample's existing capped pileup:
#   outside the globin loci   rows must be identical: same reads, cap never reached
#   inside                    per site: both percentages, the difference, and how many
#                             standard errors that is; summarised for deep (>= 1000 reads)
#                             and thin sites
#
#   sbatch scripts/test_random_subsample.sh <sample results dir> <sample id> <Stage 1 BAM> <out dir>
#
# Environment: CAP (1000000, = modkit_max_depth), PAD (50000: globin genes closer than this
# form one locus), SEED (42), SIF (the mods image), REUSE=1 (a validated random.raw.bed in the
# out dir is used as it is: only the comparison runs), REUSE=bam (random.bam is kept, the
# pileup is redone).
#
# The pileup is written to the node's local disk and validated before it is moved to the out
# dir: every line must have bedMethyl's 18 columns and the line count must equal the rows
# modkit reports. (A pileup written straight to /scratch once lost ~11M rows in ~12k chunks,
# leaving cut-off lines, while modkit reported success.)
set -euo pipefail
RES=$1 S=$2 BAM=$3 OUT=$4
CAP=${CAP:-1000000} PAD=${PAD:-50000} SEED=${SEED:-42}
SIF=${SIF:-/scratch/stein.an/ont-pipeline/containers/mods.sif}
CPUS=${SLURM_CPUS_PER_TASK:-8}
run() { apptainer exec --bind /projects/Genometechlab,/scratch "$SIF" "$@"; }

MODS=$RES/mods/GRCh38/$S
GLOBIN=$RES/qc/GRCh38/regions/GRCh38.globin.bed
for f in "$MODS/$S.raw.bed" "$MODS/$S.modkit_flags.txt" "$GLOBIN" "$BAM"; do
    [[ -s $f ]] || { echo "test_random_subsample.sh: missing $f" >&2; exit 1; }
done
mkdir -p "$OUT"; cd "$OUT"
echo "== $(date '+%F %T') sample $S, cap $CAP, seed $SEED"

# 1. globin loci: genes within PAD of each other merged (the HBB and HBA clusters)
sort -k1,1 -k2,2n "$GLOBIN" | awk -F'\t' -v OFS='\t' -v p="$PAD" '
    $1 == c && $2 <= e + p { if ($3 > e) e = $3; next }
    NR > 1 { print c, s, e }
    { c = $1; s = $2; e = $3 }
    END { if (NR) print c, s, e }' > loci.bed
echo "loci:"; sed 's/^/   /' loci.bed

if [[ ${REUSE:-0} == 1 && -s random.raw.bed && -e random.raw.bed.validated ]]; then
    echo "REUSE=1: using the validated random.raw.bed; comparing only"
else
if [[ ${REUSE:-0} == bam && -s random.bam && -s random.bam.bai ]]; then
    echo "REUSE=bam: keeping random.bam ($(run samtools view -c random.bam) reads); redoing the pileup"
else

# 2. reads off the loci as they are; each locus's reads randomly subsampled to <= CAP.
#    samtools -s keeps a read by a hash of its name and the seed: random and reproducible.
run samtools view -@ "$CPUS" -b -L loci.bed -U rest.bam -o loci_all.bam "$BAM"
run samtools index -@ "$CPUS" loci_all.bam
parts=(rest.bam); i=0
while read -r c s e; do
    i=$((i + 1)); reg="$c:$((s + 1))-$e"
    n=$(run samtools view -c loci_all.bam "$reg")
    frac=$(awk -v cap="$CAP" -v n="$n" 'BEGIN { printf "%.6f", (n > cap ? cap / n : 1) }')
    echo "locus $reg: $n reads, keeping fraction $frac"
    if [[ $frac == 1.000000 ]]; then
        run samtools view -b -o "locus$i.bam" loci_all.bam "$reg"
    else
        run samtools view -b -s "$SEED${frac#0}" -o "locus$i.bam" loci_all.bam "$reg"
    fi
    parts+=("locus$i.bam")
done < loci.bed
run samtools merge -@ "$CPUS" -f -o random.bam "${parts[@]}"
run samtools index -@ "$CPUS" random.bam
rm -f loci_all.bam loci_all.bam.bai rest.bam locus*.bam
echo "reads: original $(run samtools view -c "$BAM"), subsampled $(run samtools view -c random.bam)"

fi

# 3. pileup on the subsampled BAM, exactly as MODKIT_PILEUP runs it, written to local disk
rm -f random.raw.bed random.raw.bed.validated
LOCAL=$(mktemp -d "${TMPDIR:-/tmp}/subsample-test.XXXXXX")
echo "== $(date '+%F %T') pileup on random.bam -> $LOCAL ($(df -h --output=avail "$LOCAL" | tail -1) free there)"
run modkit pileup random.bam "$LOCAL/random.raw.bed" --threads "$CPUS" --max-depth "$CAP" \
    $(cat "$MODS/$S.modkit_flags.txt") --log-filepath random.pileup.log --suppress-progress 2> random.pileup.stderr
cat random.pileup.stderr
reported=$(grep -o 'processed [0-9]* rows' random.pileup.stderr | grep -o '[0-9]*' | tail -1)
read -r lines bad < <(awk -F'\t' 'NF != 18 { b++ } END { print NR, b + 0 }' "$LOCAL/random.raw.bed")
echo "pileup check: modkit reports ${reported:-?} rows; the file has $lines lines, $bad not 18 columns"
if [[ -z $reported || $lines != "$reported" || $bad != 0 ]]; then
    echo "PILEUP OUTPUT INVALID: left in $LOCAL for inspection; stopping" >&2; exit 1
fi
cp "$LOCAL/random.raw.bed" random.raw.bed.partial && mv random.raw.bed.partial random.raw.bed
read -r lines2 bad2 < <(awk -F'\t' 'NF != 18 { b++ } END { print NR, b + 0 }' random.raw.bed)
[[ $lines2 == "$lines" && $bad2 == 0 ]] || { echo "COPY TO $OUT INVALID ($lines2 lines, $bad2 bad); stopping" >&2; exit 1; }
rm -rf "$LOCAL"; touch random.raw.bed.validated
echo "== $(date '+%F %T') pileup done and validated ($lines rows)"
fi

# 4a. outside the loci the two pileups must agree site for site. Compared by site, not by file
#     order (modkit's threads need not write the same order twice); every differing row is
#     binned by its distance to the nearest locus, so an edge effect (reads overlapping a locus
#     were subsampled with it) shows up as differences close to the loci and nowhere else.
echo "== $(date '+%F %T') comparing outside the loci, site by site"
awk -F'\t' '
    function dist_of(chrom, pos,    i, d, best) {
        best = -1
        for (i = 1; i <= n; i++) if (chrom == c[i]) {
            if (pos >= s[i] && pos < e[i]) return -2
            d = pos < s[i] ? s[i] - pos : pos - e[i] + 1; if (best < 0 || d < best) best = d }
        return best }
    function bin(d) { return d == -1 ? "other chromosome" : (d < 100000 ? "< 100 kb" : (d < 1000000 ? "100 kb - 1 Mb" : "> 1 Mb")) }
    FNR == 1 { f++ }
    f == 1 { n++; c[n] = $1; s[n] = $2; e[n] = $3; next }
    { d = dist_of($1, $2); if (d == -2) next
      k = $1 FS $2 FS $4 FS $6; v = $10 FS $12
      if (f == 2) { A[k] = v; next }
      if (!(k in A)) { only_r++; B["only in random, " bin(d)]++; next }
      if (A[k] != v) { diff++; B["different counts, " bin(d)]++ } else same++
      delete A[k] }
    END { for (k in A) { split(k, q, FS); only_f++; B["only in first, " bin(dist_of(q[1], q[2]))]++ }
          far = 0; for (b in B) if (b !~ /< 100 kb/) far += B[b]
          printf "OUTSIDE THE LOCI: same %d, different counts %d, only in first %d, only in random %d\n", same, diff, only_f, only_r
          for (b in B) printf "   %-40s %d\n", b, B[b]
          print (far == 0 ? "OUTSIDE THE LOCI: every difference is within 100 kb of a locus (edge effect only)" \
                          : "OUTSIDE THE LOCI: " far " differing rows are 100 kb or more from any locus: NOT only an edge effect") }
' loci.bed "$MODS/$S.raw.bed" random.raw.bed

# 4b. inside: site by site, from the counts (Nmod col 12 / Nvalid col 10), not the two-decimal
#     percent column, which makes 0.01 vs 0.02 look like several standard errors at 1M reads
inside() { awk -F'\t' 'NR == FNR { n++; c[n] = $1; s[n] = $2; e[n] = $3; next }
                       { for (i = 1; i <= n; i++) if ($1 == c[i] && $2 >= s[i] && $2 < e[i]) { print; next } }' loci.bed "$1"; }
inside "$MODS/$S.raw.bed" > first.loci.bed
inside random.raw.bed > random.loci.bed
{
    printf 'class\tchrom\tstart\tmod\tstrand\tn_first\tpct_first\tn_random\tpct_random\tdifference\tz\n'
    awk -F'\t' '
        FNR == 1 { f++ }
        f == 1 { a[$1 FS $2 FS $4 FS $6] = $10 FS $12; next }
        { k = $1 FS $2 FS $4 FS $6; if (!(k in a)) next
          split(a[k], x, FS); n1 = x[1]; m1 = x[2]; n2 = $10; m2 = $12; if (n1 == 0 || n2 == 0) next
          p1 = 100 * m1 / n1; p2 = 100 * m2 / n2; pp = (m1 + m2) / (n1 + n2)
          se = 100 * sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2)); d = p2 - p1; z = se > 0 ? d / se : 0
          cls = n2 < 1000 ? "thin" : (pp >= 0.05 ? "deep_modified" : "deep_unmodified")
          printf "%s\t%s\t%s\t%s\t%s\t%d\t%.3f\t%d\t%.3f\t%+.3f\t%+.2f\n", cls, $1, $2, $4, $6, n1, p1, n2, p2, d, z
          N[cls]++; D[cls] += d; A[cls] += (d < 0 ? -d : d); if (z > 3 || z < -3) Z[cls]++; if (d > 0) P[cls]++; if (d < 0) M[cls]++ }
        END { for (c in N) printf "SUMMARY %-15s sites: %6d  mean diff %+.3f  mean |diff| %.3f  |z|>3: %d  random higher/lower: %d/%d\n", c, N[c], D[c] / N[c], A[c] / N[c], Z[c] + 0, P[c] + 0, M[c] + 0 > "/dev/stderr" }
    ' first.loci.bed random.loci.bed
} > compare_loci.tsv
echo "per-site table: $OUT/compare_loci.tsv"
echo "== $(date '+%F %T') done"
