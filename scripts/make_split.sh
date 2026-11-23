#!/usr/bin/env bash
# make_split.sh (v4): split each sample's pod5 files into chunks of about CHUNK_GB, each basecalled as
# its own one-GPU job, and write everything the three comparison runs (A, B, C) need, in $B/split:
#
#   chunks/<sample>/cNN/   symlinks to each chunk's pod5 files
#   samplesheet.csv        one row per chunk: sample_id <sample>, run_id <sample>_cNN; the pipeline merges
#                          each sample's chunks back into one BAM
#   split_<part>.config    one config per Explorer partition set, one run on each, so each run stays under
#                          its partitions' per-user QoS: short = gpu-short,sharing (Slurm starts a chunk in
#                          whichever frees first; 2 + 2 running), gpu (4 running), multigpu (4 running).
#                          One GPU per chunk, the GPU TYPE fixed per chunk (H200, or an 80 GB A100 for chunks
#                          numbered above H200_CHUNKS), so every read is basecalled on the same hardware in
#                          runs A, B and C. TIME must fit the shortest partition used (sharing: 1 h), which
#                          is why chunks are small. Includes the cohort config (node exclusions, retries)
#   {A,B,C}_*.yaml         params files: A and C load the dorado 2.1.2 batch-size benchmark file
#                          (deterministic batch size per GPU type); B pins --batchsize (dorado 1.3.0 can't
#                          read one)
#
# Samples come from SAMPLES, a file of "sample_id<TAB>pod5_folder" lines (several lines per sample are fine;
# folders are searched recursively), or for one sample from SAMPLE and POD5. Chunks are contiguous runs of
# files in path order, so the split is the same every time. Refuses to write anything if any pod5 file
# changed in the last STABLE_MIN minutes (a copy still running), or if the benchmark file lacks timings
# for a GPU type in use (that GPU would fall back to live, non-deterministic timing).
#
# Usage: SAMPLES=samples.tsv B=/scratch/stein.an/HEK293_basecall bash scripts/make_split.sh
set -euo pipefail

B=${B:-/scratch/stein.an/GM12878_basecall}
SAMPLES=${SAMPLES:-}
SAMPLE=${SAMPLE:-GM12878}
POD5=${POD5:-/scratch/stein.an/GM12878_R9RNA_pod5/pod5}
MD=${MD:-/projects/Genometechlab/andrew/custom_models}
BEN=${BEN:-$MD/batchsize_benchmarks}
COHORT_CONFIG=${COHORT_CONFIG:-/scratch/stein.an/ont-pipeline/cohort-v1/cohort.config}
CHUNK_GB=${CHUNK_GB:-12}               # ~8 min of v6 basecalling on an H200; fits TIME with modification models
H200_CHUNKS=${H200_CHUNKS:-1000}       # per sample: chunks 1..H200_CHUNKS on H200, any beyond on 80 GB A100s
B_BATCHSIZE=${B_BATCHSIZE:-960}        # run B: the v5.2.0 optimum on H200; fits an 80 GB A100
TIME=${TIME:-1h}                       # per chunk; sharing allows 1 h, gpu-short 2 h
MEMORY=${MEMORY:-32.GB}
FORKS=${FORKS:-4}                      # submitted chunks per run; within every partition's per-user limit
                                       # (sharing allows 4 submitted, the others 8), even if a job sent to
                                       # gpu-short,sharing counts against both
STABLE_MIN=${STABLE_MIN:-10}
SAMPLES_PER_BYTE=${SAMPLES_PER_BYTE:-1.087}
SPLIT=$B/split
V2=RNA004_SUP_v6.0.0_v2MODS
M6=rna004_sup@v6.0.0
M5=rna004_130bps_sup@v5.2.0

die() { echo "make_split.sh: $*" >&2; exit 1; }
[[ -f $COHORT_CONFIG ]] || die "no cohort config at $COHORT_CONFIG"
mkdir -p "$SPLIT"

# ---- 1. samples: sample_id <tab> folder
if [[ -n $SAMPLES ]]; then
    [[ -s $SAMPLES ]] || die "SAMPLES file $SAMPLES is empty or missing"
    grep -v '^[[:space:]]*\(#\|$\)' "$SAMPLES" | tr -s ' ' '\t' >| "$SPLIT/samples.tsv"
else
    printf '%s\t%s\n' "$SAMPLE" "$POD5" >| "$SPLIT/samples.tsv"
fi
: >| "$SPLIT/files.tsv"
while IFS=$'\t' read -r s dir _rest; do
    [[ -n $s && -n $dir ]] || die "bad line in samples: '$s $dir'"
    [[ $s =~ ^[A-Za-z0-9_.-]+$ ]] || die "sample id '$s' must be letters, digits, _ . - only"
    [[ -d $dir ]] || die "no folder $dir (sample $s)"
    recent=$(find "$dir" -name '*.pod5' -mmin "-$STABLE_MIN" | head -3)
    [[ -z $recent ]] || die "pod5 files under $dir changed in the last $STABLE_MIN min (copy still running?), e.g. $recent"
    find "$dir" -name '*.pod5' -printf "$s\t%s\t%p\n" >> "$SPLIT/files.tsv"
done < "$SPLIT/samples.tsv"
[[ -s $SPLIT/files.tsv ]] || die "no .pod5 files found"
dups=$(cut -f3 "$SPLIT/files.tsv" | sort | uniq -d | head -3)
[[ -z $dups ]] || die "pod5 files listed under more than one sample or folder, e.g. $dups"

# ---- 2. benchmark file for dorado 2.1.2: must cover each GPU type in use, for both models
shopt -s nullglob
parts=("$BEN"/d2.1.2_*.csv)
(( ${#parts[@]} )) || die "no $BEN/d2.1.2_*.csv benchmark files"
cat "${parts[@]}" >| "$BEN/bench_d212.csv"

# ---- 3. chunks per sample: n = round(sample bytes / CHUNK_GB), contiguous in path order
rm -rf "$SPLIT/chunks"
mkdir -p "$SPLIT/chunks"
sort -t$'\t' -k1,1 -k3,3 "$SPLIT/files.tsv" \
| awk -F'\t' -v OFS='\t' -v target="$CHUNK_GB" '
    { s[NR] = $1; b[NR] = $2; p[NR] = $3; tot[$1] += $2; nf[$1]++ }
    END {
        for (i = 1; i <= NR; i++) {
            x = s[i]
            if (!(x in n)) { n[x] = int(tot[x] / (target * 1e9) + 0.5); if (n[x] < 1) n[x] = 1; if (n[x] > nf[x]) n[x] = nf[x] }
            mid = cum[x] + b[i] / 2; cum[x] += b[i]
            c = int(mid * n[x] / tot[x]) + 1; if (c > n[x]) c = n[x]
            printf "%s\t%02d\t%s\t%s\n", x, c, b[i], p[i]
        }
    }' >| "$SPLIT/assignment.tsv"
while IFS=$'\t' read -r s c _size path; do
    mkdir -p "$SPLIT/chunks/$s/c$c"
    ln -s "$path" "$SPLIT/chunks/$s/c$c/"
done < "$SPLIT/assignment.tsv"

{
    echo "sample_id,run_id,pod5"
    cut -f1,2 "$SPLIT/assignment.tsv" | sort -u | while IFS=$'\t' read -r s c; do
        echo "$s,${s}_c$c,$SPLIT/chunks/$s/c$c"
    done
} >| "$SPLIT/samplesheet.csv"

# ---- 4. GPU type per chunk, fixed; check the benchmark file covers it
a100_runs=$(cut -f1,2 "$SPLIT/assignment.tsv" | sort -u | awk -F'\t' -v h="$H200_CHUNKS" '$2 + 0 > h { printf "%s'\''%s_c%s'\''", (k++ ? "," : ""), $1, $2 }')
need=("NVIDIA H200")
[[ -n $a100_runs ]] && need+=("A100")
for gpu in "${need[@]}"; do
    for m in "$M6" "$M5"; do
        awk -F, -v g="$gpu" -v m="$m" 'index($1, g) && $2 == m { found = 1 } END { exit !found }' "$BEN/bench_d212.csv" \
            || die "bench_d212.csv has no $gpu timings for $m (wait for that benchmark, or keep every chunk on H200)"
    done
done
if [[ -n $a100_runs ]] && [[ $(awk -F, 'index($1, "A100") { print $1 }' "$BEN/bench_d212.csv" | sort -u | wc -l) -gt 1 ]]; then
    die "benchmarks cover more than one A100 kind; keep A100 chunks on one kind"
fi
nf_time=${TIME//h/.h}; nf_time=${nf_time//min/.min}
for pair in short:gpu-short,sharing gpu:gpu multigpu:multigpu; do
    name=${pair%%:*}; parts=${pair#*:}
    cat >| "$SPLIT/split_$name.config" <<EOF
// Written by make_split.sh. One GPU per chunk; the GPU type is fixed per chunk, the same in every run.
// One run per partition set; maxForks keeps it under the per-user QoS limits of $parts.
includeConfig '$COHORT_CONFIG'

process {
    withLabel: gpu_dorado {
        queue          = '$parts'
        clusterOptions = { meta.run_id in [$a100_runs] ? '--gres=gpu:a100:1 --constraint=a100@80g' : '--gres=gpu:h200:1' }
        cpus           = 8
        memory         = $MEMORY
        time           = $nf_time
        maxForks       = $FORKS
    }
}
EOF
done

# ---- 5. params files
cat >| "$SPLIT/A_d212_v6.yaml" <<EOF
input: $SPLIT/samplesheet.csv
stage1_image: dorado-2.1.2-dna.sif
models_dir: $MD
dorado_model: $M6
dorado_args: '--batchsize-benchmarks-file $BEN/bench_d212.csv'
mod_models:
  - $V2/${M6}_inosine_m6A_2OmeA@v2
  - $V2/${M6}_m5C_2OmeC@v2
  - $V2/${M6}_pseU_2OmeU@v1
  - ${M6}_2OmeG@v1
EOF
cat >| "$SPLIT/B_d130_v52.yaml" <<EOF
input: $SPLIT/samplesheet.csv
stage1_image: dorado-1.3.0-rna.sif
dorado_model: $M5
dorado_args: '--batchsize $B_BATCHSIZE'
mod_models: [inosine_m6A_2OmeA@v1, pseU_2OmeU@v1, m5C_2OmeC@v1, 2OmeG@v1]
EOF
cat >| "$SPLIT/C_d212_v52.yaml" <<EOF
input: $SPLIT/samplesheet.csv
stage1_image: dorado-2.1.2-dna.sif
models_dir: $MD/v5.2.0
dorado_model: $M5
dorado_args: '--batchsize-benchmarks-file $BEN/bench_d212.csv'
mod_models: [inosine_m6A_2OmeA@v1, pseU_2OmeU@v1, m5C_2OmeC@v1, 2OmeG@v1]
EOF

# ---- 6. summary: per chunk, its size, GPU and estimated basecalling time (before modification models)
speed() { awk -F, -v g="$1" -v m="$2" 'index($1, g) && $2 == m && $4 > s { s = $4 } END { print s + 0 }' "$BEN/bench_d212.csv"; }
h6=$(speed "NVIDIA H200" "$M6"); a6=$(speed "A100" "$M6")
printf '%-24s %5s %8s  %-5s %16s\n' run_id files GB GPU 'v6.0.0 basecall'
awk -F'\t' '{ k = $1 "_c" $2; n[k]++; s[k] += $3 } END { for (k in n) print k "\t" n[k] "\t" s[k] }' "$SPLIT/assignment.tsv" | sort \
| while IFS=$'\t' read -r k n s; do
      c=${k##*_c}
      if (( 10#$c > H200_CHUNKS )); then g=a100; v=$a6; else g=h200; v=$h6; fi
      awk -v k="$k" -v n="$n" -v s="$s" -v g="$g" -v v="$v" -v spb="$SAMPLES_PER_BYTE" \
          'BEGIN { est = (v > 0) ? sprintf("%.1f h", s * spb / v / 3600) : "no timing"
                   printf "%-24s %5d %8.1f  %-5s %16s\n", k, n, s / 1e9, g, est }'
  done
awk -F'\t' '{ n[$1]++; s[$1] += $3 } END { for (x in n) printf "== %s: %d files, %.0f GB\n", x, n[x], s[x] / 1e9 }' "$SPLIT/assignment.tsv" | sort
echo "== $(($(wc -l < "$SPLIT/samplesheet.csv") - 1)) chunks; written to $SPLIT: samplesheet.csv, split_{short,gpu,multigpu}.config, A/B/C params files"
echo "   the modification models add GPU time on top: keep the largest estimate under a third of TIME ($TIME)"
