#!/usr/bin/env bash
# cmp_pileup.sh: modkit 0.6.4 pileup of one aligned BAM for the model comparison (runs A, B, C):
#
#   1. per-sample thresholds, THRESH_MODE:
#        explicit      (recommended) the 10th percentile of the calls the model actually made: a random
#                      SAMPLE_FRAC of reads through modkit extract calls --no-filtering, keeping aligned,
#                      non-inferred calls. Bases dorado skipped (canonical by implication in '.' MM mode)
#                      are left out; with them, a base whose models rarely emit probabilities gets its
#                      threshold pushed to the top of the scale (seen for G and U with v5.2.0), and nearly
#                      every explicit call on it then fails.
#        fixed         THRESHOLDS="A:0.8 C:0.8 G:0.87 T:0.9" for every run alike
#        sample-probs  (v1 behaviour) modkit sample-probs on EVERY read (--no-sampling: exact, read-weighted,
#      deterministic), restricted to bases aligned to the reference (--ref + --motif X 0 for each base;
#      0.6.4 has no --only-mapped, and this is what pileup's own estimate uses), at the 10th percentile
#      (0.6.4 takes percentiles in percent). Written to <prefix>.thresholds.tsv; a threshold at the
#      top of the scale (>= 0.99) is flagged: in '.' MM mode skipped bases enter as canonical, so a
#      rarely called modification (e.g. 2'-O-Me G) can push its base's threshold to the ceiling.
#   2. random subsampling (bin/subsample_loci.sh, as in Stage 2): every 1 kb window with more than
#      0.98 x 65,535 overlapping reads keeps a random fraction of them, chosen by a hash of read
#      name and seed, so the same reads are kept in every run. Thresholds come from the original BAM.
#   3. modkit pileup --high-depth --max-depth 65535 with all eight RNA modifications
#      (--modified-bases + --ref: 0.6.4's fast workers). Without --high-depth/--max-depth, 0.6.4 caps
#      every position at its default of 8,000 reads, taken in coordinate order; after step 2 no
#      position reaches 65,535, so the cap never chooses reads.
#   4. the bedMethyl is written to node-local disk, checked (row count = modkit's report, 18 columns),
#      then copied, as bin/validated_pileup.sh does.
#
# Usage: cmp_pileup.sh <aligned.bam> <reference.fa> <out_prefix>
#   env: THREADS (16), IMAGE (dorado-2.1.2-dna.sif, carries modkit 0.6.4 + samtools), REPO
#        (~/ont-pipeline-v2, for bin/subsample_loci.sh), CAP (65535), SEED (42), PERCENTILE (10),
#        MIN_MAPQ (unset: every primary alignment, like Stage 2's modification branch)
# Outputs: <prefix>.thresholds.tsv, <prefix>.sample_probs/ (+ .txt, .log), <prefix>.loci.bed,
#          <prefix>.bed (+ .rows), <prefix>.pileup.log, <prefix>.summary.txt
set -euo pipefail

bam=${1:?usage: cmp_pileup.sh <aligned.bam> <reference.fa> <out_prefix>}
fasta=${2:?reference FASTA}
prefix=${3:?output prefix}
THREADS=${THREADS:-16}
THRESHOLDS_THREADS=${THRESHOLDS_THREADS:-$THREADS}
IMAGE=${IMAGE:-/scratch/stein.an/ont-pipeline/containers/dorado-2.1.2-dna.sif}
REPO=${REPO:-$HOME/ont-pipeline-v2}
CAP=${CAP:-65535}
SEED=${SEED:-42}
PERCENTILE=${PERCENTILE:-10}
THRESH_MODE=${THRESH_MODE:-explicit}
SAMPLE_FRAC=${SAMPLE_FRAC:-0.1}
[[ $THRESH_MODE =~ ^(explicit|fixed|sample-probs)$ ]] || { echo "cmp_pileup.sh: THRESH_MODE must be explicit, fixed or sample-probs" >&2; exit 1; }
MODS=(m6A inosine 2OmeA m5C 2OmeC pseU 2OmeU 2OmeG)

die() { echo "cmp_pileup.sh: $*" >&2; exit 1; }
say() { echo "[$(date '+%F %T')] $*" | tee -a "$prefix.summary.txt" >&2; }
run() { apptainer exec -B /scratch,/projects "$IMAGE" "$@"; }

[[ -s $bam ]] || die "no BAM $bam"
bam=$(readlink -f "$bam"); fasta=$(readlink -f "$fasta")
[[ -s $bam.bai || -s ${bam%.bam}.bai ]] || die "no index for $bam"
[[ -s $fasta && -s $fasta.fai ]] || die "need $fasta and $fasta.fai"
[[ -s $IMAGE ]] || die "no image $IMAGE"
[[ -x $REPO/bin/subsample_loci.sh || -f $REPO/bin/subsample_loci.sh ]] || die "no $REPO/bin/subsample_loci.sh"
(( CAP <= 65535 )) || die "CAP $CAP: modkit 0.6.4 counts in 16 bits, at most 65535"
mkdir -p "$(dirname "$prefix")"
prefix=$(cd "$(dirname "$prefix")" && pwd)/$(basename "$prefix")
: > "$prefix.summary.txt"
ver=$(run modkit --version | awk '{print $2}')
[[ $ver == 0.6.* ]] || die "image has modkit $ver, expected 0.6.x"
run sh -c 'command -v samtools > /dev/null && command -v awk > /dev/null' || die "image lacks samtools or awk"
say "modkit $ver; BAM $bam; reference $fasta; $THREADS threads"

# ---- 0. optional MAPQ filter (off by default, matching Stage 2)
if [[ -n ${MIN_MAPQ:-} ]]; then
    say "keeping primary alignments with MAPQ >= $MIN_MAPQ"
    run samtools view -@ "$THREADS" -b -q "$MIN_MAPQ" -F 2308 -o "$prefix.mapq.bam" "$bam"
    run samtools index -@ "$THREADS" "$prefix.mapq.bam"
    bam=$prefix.mapq.bam
fi

# ---- 1. thresholds
t0=$(date +%s)
if [[ $THRESH_MODE == fixed ]]; then
    [[ -n ${THRESHOLDS:-} ]] || die "THRESH_MODE=fixed needs THRESHOLDS=\"A:x C:x G:x T:x\""
    : > "$prefix.thresholds.tsv"
    for kv in $THRESHOLDS; do printf '%s\t%s\n' "${kv%%:*}" "${kv#*:}" >> "$prefix.thresholds.tsv"; done
elif [[ $THRESH_MODE == explicit ]]; then
    frac=${SAMPLE_FRAC#0}
    run bash -c "samtools view -u -@ 2 -s ${SEED}${frac} '$bam' \
        | modkit extract calls - - --mapped-only --no-filtering --threads $THRESHOLDS_THREADS \
            --log-filepath '$prefix.extract_calls.log' --suppress-progress" \
    | awk -F'\t' -v OFS='\t' '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i
                  if (!("call_prob" in col) || !("canonical_base" in col) || !("inferred" in col)) {
                      print "missing columns in modkit extract calls header: " $0 > "/dev/stderr"; exit 3 }
                  next }
        ("ref_position" in col) && $col["ref_position"] < 0 { next }
        ("within_alignment" in col) && $col["within_alignment"] == "false" { next }
        $col["inferred"] == "true" { inf[$col["canonical_base"]]++; next }
        { n[$col["canonical_base"] "\t" $col["call_prob"]]++; tot[$col["canonical_base"]]++ }
        END { for (k in n) print k, n[k]
              for (b in tot) printf "explicit calls %s: %d, inferred (left out): %d\n", b, tot[b], inf[b] > "/dev/stderr" }' \
        2> "$prefix.explicit_counts.txt" | sort -t$'\t' -k1,1 -k2,2g > "$prefix.explicit_hist.tsv" \
        || die "explicit-call thresholds failed; see $prefix.extract_calls.log and $prefix.explicit_counts.txt"
    awk -F'\t' -v OFS='\t' -v pct="$PERCENTILE" '
        { b = $1; cnt[b] += $3; p[b, ++m[b]] = $2; c[b, m[b]] = $3 }
        END { for (b in cnt) { need = cnt[b] * pct / 100; s = 0
                  for (i = 1; i <= m[b]; i++) { s += c[b, i]; if (s >= need) { print b, p[b, i]; break } } } }' \
        "$prefix.explicit_hist.tsv" | sort > "$prefix.thresholds.tsv"
    say "explicit-call sample ($SAMPLE_FRAC of reads, seed $SEED): $(grep -c . "$prefix.explicit_hist.tsv") probability bins; $(tr '\n' ';' < "$prefix.explicit_counts.txt")"
else
rm -rf "$prefix.sample_probs"
run modkit sample-probs "$bam" --ref "$fasta" \
    --motif A 0 --motif C 0 --motif G 0 --motif T 0 \
    --no-sampling --percentiles "$PERCENTILE" --threads "$THREADS" \
    --out-dir "$prefix.sample_probs" --prefix "$(basename "$prefix")" --hist --force \
    --log-filepath "$prefix.sample_probs.log" --suppress-progress > "$prefix.sample_probs.txt" \
    || die "modkit sample-probs failed; see $prefix.sample_probs.log"
table=$(ls "$prefix.sample_probs"/*thresholds*.tsv 2> /dev/null | head -1 || true)
src=${table:-$prefix.sample_probs.txt}
# parse any layout: tab-separated, space-aligned or box-drawn; header names the columns
sed 's/[│┃|║─━═┌┐└┘├┤┬┴┼╭╮╰╯]/ /g' "$src" | awk -v want="$PERCENTILE" -v OFS='\t' '
    tolower($0) ~ /base/ && tolower($0) ~ /threshold/ && !hdr {
        for (i = 1; i <= NF; i++) { h = tolower($i); if (h ~ /^(base|primary_base)$/) b = i; else if (h ~ /threshold/) t = i; else if (h ~ /percentile/) p = i }
        hdr = 1; next }
    hdr && b && t && $b ~ /^[ACGT]$/ {
        if (p && !($p + 0 == want + 0 || $p * 100 == want + 0)) next
        if (!($b in seen)) { seen[$b] = 1; print $b, $t } }' > "$prefix.thresholds.tsv"
fi
[[ $(wc -l < "$prefix.thresholds.tsv") -eq 4 ]] \
    || die "did not get four per-base thresholds ($THRESH_MODE mode; got: $(tr '\n' ' ' < "$prefix.thresholds.tsv"))"
flags=()
while IFS=$'\t' read -r base thr; do
    flags+=(--filter-threshold "$base:$thr")
    awk -v t="$thr" 'BEGIN { exit !(t + 0 >= 0.99) }' && say "WARNING: $base threshold $thr is at the top of the scale (see header note)"
done < "$prefix.thresholds.tsv"
say "thresholds ($THRESH_MODE, ${PERCENTILE}th percentile): $(tr '\t\n' ': ' < "$prefix.thresholds.tsv")in $(( ($(date +%s) - t0) / 60 )) min"

# ---- 2. random subsampling of loci deeper than the cap
t0=$(date +%s)
bai=$bam.bai; [[ -s $bai ]] || bai=${bam%.bam}.bai
sub=$prefix.subsample.d; rm -rf "$sub"; mkdir -p "$sub"      # its scratch BAMs go here, not beside the outputs
( cd "$sub" && run bash "$REPO/bin/subsample_loci.sh" run "$bam" "$bai" "$CAP" 1000 1000000 "$SEED" "$THREADS" s ) \
    2> "$prefix.subsample.err" || die "subsample_loci.sh failed; see $prefix.subsample.err"
mv -f "$sub/s.loci.bed" "$prefix.loci.bed"; mv -f "$sub/s.subsample.log" "$prefix.subsample.log" 2> /dev/null || true
if [[ -s $prefix.loci.bed ]]; then
    mv -f "$sub/s.pileup.bam" "$prefix.pileup.bam"; mv -f "$sub/s.pileup.bam.bai" "$prefix.pileup.bam.bai"
    pbam=$prefix.pileup.bam
    say "subsampled $(wc -l < "$prefix.loci.bed") loci deeper than 0.98 x $CAP reads (seed $SEED) in $(( ($(date +%s) - t0) / 60 )) min"
else
    pbam=$bam
    say "no locus deeper than 0.98 x $CAP reads: pileup reads the BAM itself"
fi
rm -rf "$sub"

# ---- 3. pileup, all eight modifications, cap above any remaining depth
t0=$(date +%s)
tmp=$(mktemp -d "${TMPDIR:-/tmp}/cmp_pileup.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
run modkit pileup "$pbam" "$tmp/out.bed" --ref "$fasta" --modified-bases "${MODS[@]}" \
    --high-depth --max-depth "$CAP" "${flags[@]}" --threads "$THREADS" \
    --log-filepath "$prefix.pileup.log" --suppress-progress 2> "$tmp/stderr" \
    || { cat "$tmp/stderr" >&2; die "modkit pileup failed"; }
cat "$tmp/stderr" >> "$prefix.pileup.log"
reported=$(grep -ho 'processed [0-9]* rows' "$tmp/stderr" "$prefix.pileup.log" | grep -o '[0-9]*' | tail -n 1 || true)
[[ -n $reported ]] || die "modkit did not report how many rows it wrote"

# ---- 4. validate, copy, validate again
count() { awk -F'\t' 'NF != 18 { b++ } END { print NR, b + 0 }' "$1"; }
read -r n bad < <(count "$tmp/out.bed")
[[ $n == "$reported" && $bad == 0 ]] || die "local output: $n lines ($bad not 18 columns), modkit reports $reported"
cp "$tmp/out.bed" "$prefix.bed.partial"
read -r n2 bad2 < <(count "$prefix.bed.partial")
[[ $n2 == "$n" && $bad2 == 0 ]] || die "copy: $n2 lines ($bad2 not 18 columns), expected $n"
mv "$prefix.bed.partial" "$prefix.bed"; echo "$n" > "$prefix.bed.rows"
grep -q 'using maximum of' "$prefix.pileup.log" || say "note: modkit did not log its depth cap; check $prefix.pileup.log"
say "pileup: $n rows in $(( ($(date +%s) - t0) / 60 )) min -> $prefix.bed"
