#!/usr/bin/env bash
# Random subsampling before modkit pileup, so that its --max-depth cap never decides which
# reads count. With a cap, modkit keeps the first reads at each position in alignment order,
# not a random set; on whole blood that biased the globin site values, by a few points at
# times even at ~450k reads (ST1918, 2026-10-27, scripts/test_random_subsample.sh). So each
# locus deeper than the cap is randomly subsampled beforehand to just under it, and every
# other read is kept as it is.
#
#   subsample_loci.sh run <bam> <bai> <cap> <window> <merge> <seed> <threads> <prefix>
#       -> <prefix>.loci.bed        chrom start end deepest_window_reads kept_fraction (may be empty)
#          <prefix>.pileup.bam/.bai the BAM modkit pileup reads (the original, linked, if no locus)
#          <prefix>.subsample.log
#   subsample_loci.sh loci <cap> <window> <merge> < "chrom<TAB>pos<TAB>cigar" (position-sorted)
#       -> chrom start end deepest_window_reads, one line per locus (what run uses)
#
# Depth is the number of reads overlapping each <window>-bp window, a read covering its whole
# aligned span including spliced gaps: modkit's depth limit counts those reads too. Windows
# deeper than the cap become loci, merged when within <merge> bp. A locus keeps a random
# 0.98 * cap / (its deepest window) of the reads that overlap it, chosen by samtools -s from a
# hash of read name and seed (random and reproducible), which brings its deepest window just
# under the cap. A read overlapping two loci belongs to the first, so none is duplicated.
# A dropped read is dropped wherever it aligns, so sites near a locus also see fewer reads:
# unbiased there, just less deep.
set -euo pipefail

loci() {
    local cap=$1 window=$2 merge=$3
    awk -F'\t' -v OFS='\t' -v W="$window" -v cap="$cap" -v M="$merge" '
        function emit(b, n,    bs, be) {
            bs = b * W; be = bs + W
            if (lc == cur && bs <= le + M) { if (be > le) le = be; if (n > lm) lm = n }
            else { if (lc != "") print lc, ls, le, lm; lc = cur; ls = bs; le = be; lm = n } }
        function flush_to(limit,    b) {
            for (b = lo; b < limit; b++) if (b in cnt) { if (cnt[b] > cap) emit(b, cnt[b]); delete cnt[b] }
            lo = limit }
        $1 != cur { flush_to(hi + 1); cur = $1; lo = 0; hi = -1 }
        { reflen = 0; nl = split($3, L, /[MIDNSHP=X]/); split($3, O, /[0-9]+/)
          for (i = 1; i < nl; i++) if (O[i + 1] ~ /^[MDN=X]$/) reflen += L[i]
          b0 = int(($2 - 1) / W); b1 = int(($2 + reflen - 2) / W); if (b1 < b0) b1 = b0
          if (b0 > lo) flush_to(b0)
          for (b = b0; b <= b1; b++) cnt[b]++
          if (b1 > hi) hi = b1; reads++ }
        END { flush_to(hi + 1); if (lc != "") print lc, ls, le, lm
              printf "subsample_loci.sh: %d reads counted in %d-bp windows\n", reads, W > "/dev/stderr" }'
}

run() {
    local bam=$1 bai=$2 cap=$3 window=$4 merge=$5 seed=$6 threads=$7 prefix=$8
    local log="$prefix.subsample.log"
    : > "$log"
    note() { echo "$*" | tee -a "$log" >&2; }
    if (( cap >= 2147483646 )); then
        : > "$prefix.loci.bed"
        note "cap $cap: modkit pileup is not capped, so nothing needs subsampling"
    else
        samtools view -@ "$threads" -F 2308 "$bam" | cut -f 3,4,6 | loci "$cap" "$window" "$merge" \
            | awk -F'\t' -v OFS='\t' -v cap="$cap" '{ printf "%s\t%s\t%s\t%s\t%.6f\n", $1, $2, $3, $4, 0.98 * cap / $4 }' \
            > "$prefix.loci.bed"
    fi
    if [[ ! -s $prefix.loci.bed ]]; then
        note "no window deeper than $cap reads: modkit pileup reads the original BAM"
        ln -sf "$(readlink -f "$bam")" "$prefix.pileup.bam"
        ln -sf "$(readlink -f "$bai")" "$prefix.pileup.bam.bai"
        return 0
    fi

    note "loci deeper than $cap reads (subsampled at random, seed $seed):"
    while IFS=$'\t' read -r c s e m f; do note "   $c:$((s + 1))-$e  deepest window $m reads, keeping fraction $f"; done < "$prefix.loci.bed"
    cut -f 1-3 "$prefix.loci.bed" > loci3.bed
    # every record either overlaps a locus (pending) or not (rest, kept as it is)
    samtools view -@ "$threads" -b -L loci3.bed -U rest.bam -o pending.bam "$bam"
    local parts=(rest.bam) i=0
    while IFS=$'\t' read -r c s e m f; do
        i=$((i + 1))
        printf '%s\t%s\t%s\n' "$c" "$s" "$e" > locus.bed
        # this locus takes the pending reads that overlap it; the rest wait for the next locus
        samtools view -@ "$threads" -b -L locus.bed -U next.bam -o "take$i.bam" pending.bam
        samtools view -@ "$threads" -b -s "$seed${f#0}" -o "part$i.bam" "take$i.bam"
        note "   locus $i: $(samtools view -c "take$i.bam") reads overlap it, $(samtools view -c "part$i.bam") kept"
        rm -f "take$i.bam"; mv next.bam pending.bam; parts+=("part$i.bam")
    done < "$prefix.loci.bed"
    local left; left=$(samtools view -c pending.bam)
    if (( left > 0 )); then note "   $left reads matched a locus overall but none singly: kept as they are"; parts+=(pending.bam); fi
    samtools merge -@ "$threads" -f -o "$prefix.pileup.bam" "${parts[@]}"
    samtools index -@ "$threads" "$prefix.pileup.bam"
    rm -f rest.bam pending.bam part*.bam locus.bed loci3.bed
    note "pileup BAM: $(samtools view -c "$prefix.pileup.bam") of $(samtools view -c "$bam") records"
}

case ${1:-} in
    loci) shift; loci "$@" ;;
    run)  shift; run "$@" ;;
    *)    echo "usage: $0 {run <bam> <bai> <cap> <window> <merge> <seed> <threads> <prefix> | loci <cap> <window> <merge>}" >&2; exit 1 ;;
esac
