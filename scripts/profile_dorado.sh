#!/usr/bin/env bash
# profile_dorado.sh: for each dorado image, does the model load, and how fast is it on one GPU against
# all of them? Then: how long would the whole dataset, and each flow cell, take?
#
# Run inside a GPU allocation of the kind you will queue on, for example:
#   salloc -p multigpu --gres=gpu:h200:4 -c 32 --mem=128G -t 2:00:00
#   MODEL=... MODS=... srun --pty bash scripts/profile_dorado.sh
#
# Settings (environment):
#   MODEL    basecall model folder, relative to MODELS (required)
#   MODS     modification model folders, relative to MODELS, comma-separated (optional)
#   MODELS   folder holding the models (default /projects/Genometechlab/andrew/custom_models)
#   POD5     pod5 data, searched recursively; each subfolder is taken as one flow cell
#   IMG_DIR  folder with the images; IMAGES: image files in it, space-separated
#   FILES    pod5 files basecalled per test (default 10, ~40,000 reads at MinKNOW's 4,000 per file)
#   OUT      where BAMs, logs and report.txt go
# Throughput includes dorado's start-up (model load, batch-size search), so estimates lean long.
set -euo pipefail

MODEL=${MODEL:?set MODEL to the basecall model folder, relative to MODELS}
MODS=${MODS:-}
MODELS=${MODELS:-/projects/Genometechlab/andrew/custom_models}
POD5=${POD5:-/scratch/stein.an/GM12878_R9RNA_pod5}
IMG_DIR=${IMG_DIR:-/scratch/stein.an/ont-pipeline/containers}
IMAGES=${IMAGES:-"dorado-1.3.0-rna.sif dorado-2.1.2-dna.sif"}
FILES=${FILES:-10}
OUT=${OUT:-/scratch/stein.an/dorado_profile}

for t in apptainer nvidia-smi; do
    command -v "$t" > /dev/null || { echo "ERROR: $t not on PATH (run inside a GPU allocation)" >&2; exit 1; }
done
[[ -d $MODELS/$MODEL ]] || { echo "ERROR: no model folder $MODELS/$MODEL" >&2; exit 1; }
mods=""
if [[ -n $MODS ]]; then
    IFS=',' read -ra ml <<< "$MODS"
    for m in "${ml[@]}"; do
        [[ -d $MODELS/$m ]] || { echo "ERROR: no model folder $MODELS/$m" >&2; exit 1; }
        mods+="${mods:+,}$MODELS/$m"
    done
    mods="--modified-bases-models $mods"
fi
mkdir -p "$OUT"
exec > >(tee "$OUT/report.txt")

ngpu=$(nvidia-smi -L | wc -l)
echo "== $(hostname): $ngpu GPU(s), $(nproc) CPUs, $(date '+%F %T')"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader

mapfile -t all < <(find "$POD5" -name '*.pod5' | sort)
(( ${#all[@]} )) || { echo "ERROR: no .pod5 files under $POD5" >&2; exit 1; }
sub=$OUT/subset
rm -rf "$sub"
mkdir -p "$sub"
for f in "${all[@]:0:$FILES}"; do ln -s "$f" "$sub/"; done
sub_bytes=$(find -L "$sub" -name '*.pod5' -printf '%s\n' | awk '{ s += $1 } END { print s }')
tot_bytes=$(find "$POD5" -name '*.pod5' -printf '%s\n' | awk '{ s += $1 } END { print s }')
echo "== test set: $(ls "$sub" | wc -l) of ${#all[@]} pod5 files, $((sub_bytes / 1000000)) of $((tot_bytes / 1000000)) MB"
echo "== model: $MODEL${MODS:+ + $MODS}"
echo

devices="cuda:0"
(( ngpu > 1 )) && devices="cuda:0 cuda:all"
printf '%-24s %-9s %5s %7s %8s %11s %10s  %s\n' image device exit seconds reads samples/s hours_all 'GPU use (mean)'
best_bps=0
for img in $IMAGES; do
    if [[ ! -s $IMG_DIR/$img ]]; then
        printf '%-24s missing in %s\n' "$img" "$IMG_DIR"
        continue
    fi
    for dev in $devices; do
        tag="${img%.sif}.${dev/:/_}"
        nvidia-smi --query-gpu=index,utilization.gpu --format=csv,noheader,nounits -l 5 > "$OUT/$tag.gpu.csv" &
        mon=$!
        t0=$(date +%s)
        set +e
        # shellcheck disable=SC2086
        apptainer exec --nv -B "$MODELS" -B "$POD5" -B "$OUT" "$IMG_DIR/$img" \
            dorado basecaller "$MODELS/$MODEL" "$sub" --recursive --device "$dev" $mods \
            > "$OUT/$tag.bam" 2> "$OUT/$tag.log"
        rc=$?
        set -e
        secs=$(( $(date +%s) - t0 ))
        kill "$mon" 2> /dev/null || true
        wait "$mon" 2> /dev/null || true
        reads=$(grep -o 'reads basecalled: [0-9]*' "$OUT/$tag.log" | tail -1 | grep -o '[0-9]*$' || echo '?')
        sps=$(grep -o 'Samples/s: [0-9.e+]*' "$OUT/$tag.log" | tail -1 | awk '{ print $2 }' || echo '?')
        use=$(awk -F', *' '{ u[$1] += $2; n[$1]++ } END { for (g = 0; g in u; g++) printf "%s%d:%.0f%%", (g ? " " : ""), g, u[g] / n[g] }' "$OUT/$tag.gpu.csv")
        if (( rc == 0 && secs > 0 )); then
            bps=$(( sub_bytes / secs ))
            hours=$(awk -v t="$tot_bytes" -v b="$bps" 'BEGIN { printf "%.1f", t / b / 3600 }')
            (( bps > best_bps )) && best_bps=$bps best_tag=$tag
        else
            hours='-'
        fi
        printf '%-24s %-9s %5s %7s %8s %11s %10s  %s\n' "$img" "$dev" "$rc" "$secs" "$reads" "$sps" "$hours" "$use"
        if (( rc != 0 )); then
            echo "   last lines of $OUT/$tag.log:"
            tail -n 4 "$OUT/$tag.log" | sed 's/^/   | /'
        fi
    done
done

echo
if (( best_bps == 0 )); then
    echo "== no run finished; read the logs in $OUT"
    exit 1
fi
echo "== estimated hours per flow cell (subfolder of $POD5) at the fastest setting ($best_tag)"
find "$POD5" -name '*.pod5' -printf '%h\t%s\n' \
  | awk -F'\t' -v root="${POD5%/}/" -v bps="$best_bps" '
        { k = "."
          if (index($1 "/", root) == 1 && length($1) + 1 > length(root)) { split(substr($1, length(root) + 1), p, "/"); k = p[1] }
          s[k] += $2; n[k]++ }
        END { for (k in s) printf "  %-40s %6d files %9.1f GB %6.1f h\n", k, n[k], s[k] / 1e9, s[k] / bps / 3600 }' \
  | sort