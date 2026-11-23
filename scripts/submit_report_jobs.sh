#!/usr/bin/env bash
# submit_report_jobs.sh: submits every job still needed for the v6 vs v5.2 report. Sets its own paths, checks them,
# and stops at the first problem. Run from a shell where bioinf_tools is active:
#     mamba activate bioinf_tools && bash ~/ont-pipeline-v2/scripts/submit_report_jobs.sh
set -euo pipefail
SC=~/ont-pipeline-v2/scripts; B=/scratch/stein.an/basecall_cmp; PY=$(command -v python)
FA=/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
MASK=/scratch/stein.an/ont-pipeline/cohort-v1/psi_loci_check/clusters.bed
IMGD=/scratch/stein.an/ont-pipeline/containers; U=$B/tower_ubam; O=$B/sampleprobs_050; PL=$B/pileup064_lab050
P4=$B/pileup064_nofilter; ORTHO=/projects/Genometechlab/andrew/orthogonal
S_HEK="08_01_25_RNA004_HEK293 09_09_25_RNA004_HEK293 10_06_25_RNA004_HEK293_CP"; R="A_d212_v6 B_d130_v52 C_d212_v52"

# ---- checks before anything is submitted
"$PY" -c 'import importlib.util as u, sys; sys.exit(not all(u.find_spec(m) for m in ("pandas", "matplotlib", "openpyxl", "pyxlsb")))' || { echo "activate bioinf_tools first (python: $PY)"; exit 1; }
for f in "$SC/cmp_ortho.py" "$SC/cmp_plots.py" "$SC/cmp_pileup.sh" "$FA" "$FA.fai" "$MASK" "$IMGD/mods.sif" "$IMGD/dorado-2.1.2-dna.sif" "$ORTHO/MD5SUMS"; do
    [[ -s $f ]] || { echo "missing: $f"; exit 1; }
done
[[ $(md5sum < "$SC/cmp_ortho.py" | cut -d' ' -f1) == fef97ea2852dbdcb1aed08845de95157 ]] || { echo "cmp_ortho.py is not the latest version"; exit 1; }
for s in GM12878 $S_HEK; do for r in $R; do
    [[ -s $U/out_$r/bam/GRCh38/$s.GRCh38.bam ]] || { echo "missing BAM: $r $s"; exit 1; }
    [[ $s == GM12878 ]] || [[ -s $P4/$r/$s.bed ]] || { echo "missing threshold-free pileup: $r $s"; exit 1; }
done; done
for r in $R; do [[ -s $O/$r/thresholds.tsv ]] || { echo "missing GM12878 0.5.0 thresholds: $O/$r/thresholds.tsv"; exit 1; }; done
mkdir -p "$PL/logs" "$B/ortho/logs" "$B/plots/logs" "$B/sp064_diag"
sub() { sbatch --parsable "$@"; }
declare -A PILE                                                     # lab-arm pileup job id per run/sample
log() { printf '%-34s %s\n' "$1" "$2"; }

# ---- wrapper: modkit 0.6.4 pileup at the thresholds modkit 0.5.0 estimated (the lab's method)
cat > "$B/lab_pileup.sh" <<'IN'
#!/usr/bin/env bash
set -euo pipefail
th=$(awk 'NR > 1 && $1 ~ /^[ACGT]$/ { printf "%s:%s ", $1, $3 }' "$1")
[[ $(wc -w <<< "$th") -eq 4 ]] || { echo "lab_pileup.sh: expected 4 thresholds in $1, got '$th'" >&2; exit 1; }
echo "thresholds from $1: $th" >&2
THREADS=${THREADS:-8} THRESH_MODE=fixed THRESHOLDS="$th" bash ~/ont-pipeline-v2/scripts/cmp_pileup.sh "$2" "$3" "$4"
IN
chmod +x "$B/lab_pileup.sh"

echo "== 1. lab-rule arm: 0.5.0 thresholds -> 0.6.4 pileups -> plots and scoring"
for s in GM12878 $S_HEK; do
    ids=()
    for r in $R; do
        bam=$U/out_$r/bam/GRCh38/$s.GRCh38.bam; dep=()
        if [[ $s == GM12878 ]]; then
            t=$O/$r                                                   # estimated already
        else
            t=$O/$r/$s
            j=$(sub -p short -c 8 --mem=16G -t 1:00:00 -J sp050-$r-$s -o "$O/%x.%j.log" --wrap \
                "apptainer exec -B /scratch,/projects $IMGD/mods.sif modkit sample-probs $bam --sampling-frac 0.1 --seed 42 --percentiles 0.1 -t 8 --out-dir $t --force --log-filepath $t.log")
            dep=(--dependency=afterok:$j); log "sp050-$r-$s" "$j"
        fi
        j=$(sub ${dep[@]+"${dep[@]}"} -p short -c 8 --mem=32G -t 4:00:00 -J pilelab-$r-$s -o "$PL/logs/%x.%j.log" \
            --wrap "bash $B/lab_pileup.sh $t/thresholds.tsv $bam $FA $PL/$r/$s")
        ids+=("$j"); PILE[$r/$s]=$j; log "pilelab-$r-$s" "$j"
    done
    dep=$(IFS=:; echo "${ids[*]}")
    runs="--run A_d212_v6=$PL/A_d212_v6/$s.bed --run B_d130_v52=$PL/B_d130_v52/$s.bed --run C_d212_v52=$PL/C_d212_v52/$s.bed"
    log "plotslab-$s" "$(sub --dependency=afterok:$dep -p short -c 4 --mem=64G -t 2:00:00 -J plotslab-$s -o "$PL/logs/%x.%j.log" \
        --wrap "$PY $SC/cmp_plots.py --sample $s --mask-bed $MASK --out $B/plots/${s}_lab050 $runs")"
    [[ $s == GM12878 ]] || log "ortho-lab050-$s" "$(sub --dependency=afterok:$dep -p short -c 4 --mem=48G -t 2:00:00 \
        -J ortho-lab050-$s -o "$B/ortho/logs/%x.%j.log" \
        --wrap "$PY $SC/cmp_ortho.py --sample $s --ortho-dir $ORTHO --fasta $FA --mask-bed $MASK --out $B/ortho/${s}_lab050 $runs")"
done

echo "== 2. threshold-free scoring again, with confidence intervals"
for s in $S_HEK; do
    log "ortho-nofilter-$s" "$(sub -p short -c 4 --mem=48G -t 2:00:00 -J ortho-nofilter-$s -o "$B/ortho/logs/%x.%j.log" \
        --wrap "$PY $SC/cmp_ortho.py --sample $s --ortho-dir $ORTHO --fasta $FA --mask-bed $MASK --out $B/ortho/${s}_nofilter \
                --run A_d212_v6=$P4/A_d212_v6/$s.bed --run B_d130_v52=$P4/B_d130_v52/$s.bed --run C_d212_v52=$P4/C_d212_v52/$s.bed")"
done

echo "== 3. replicate comparisons per model, both arms (lab arm after its pileups)"
r=C_d212_v52
if [[ ! -s $B/plots/HEK293_replicates_${r}_nofilter/tables/all_runs_metrics.tsv ]]; then
    log "reps-nofilter-$r" "$(sub -p short -c 4 --mem=64G -t 2:00:00 -J reps-nofilter-$r -o "$B/plots/logs/%x.%j.log" \
        --wrap "$PY $SC/cmp_plots.py --sample 'HEK293 replicates, $r (nofilter)' --mask-bed $MASK --out $B/plots/HEK293_replicates_${r}_nofilter \
                --run rep_08_01=$P4/$r/08_01_25_RNA004_HEK293.bed --run rep_09_09=$P4/$r/09_09_25_RNA004_HEK293.bed --run rep_10_06_CP=$P4/$r/10_06_25_RNA004_HEK293_CP.bed")"
fi
for r in $R; do
    deps=$(for s in $S_HEK; do echo -n "${PILE[$r/$s]}:"; done); deps=${deps%:}
    log "reps-lab050-$r" "$(sub --dependency=afterok:$deps -p short -c 4 --mem=64G -t 2:00:00 -J reps-lab050-$r \
        -o "$B/plots/logs/%x.%j.log" --wrap "$PY $SC/cmp_plots.py --sample 'HEK293 replicates, $r (lab050)' --mask-bed $MASK \
        --out $B/plots/HEK293_replicates_${r}_lab050 --run rep_08_01=$PL/$r/08_01_25_RNA004_HEK293.bed \
        --run rep_09_09=$PL/$r/09_09_25_RNA004_HEK293.bed --run rep_10_06_CP=$PL/$r/10_06_25_RNA004_HEK293_CP.bed")"
done

echo "== 4. where the 0.998 ceiling comes from (modkit 0.6.4 with 0.5.0's flags; 0.6.4 on all data, no motif restriction)"
D=$B/sp064_diag; bam=$U/out_B_d130_v52/bam/GRCh38/GM12878.GRCh38.bam
log "sp064-like050" "$(sub -p short -c 8 --mem=16G -t 1:00:00 -J sp064-like050 -o "$D/%x.%j.log" --wrap \
    "apptainer exec -B /scratch,/projects $IMGD/dorado-2.1.2-dna.sif modkit sample-probs $bam --sampling-frac 0.1 --seed 42 --percentiles 10 -t 8 --out-dir $D/like050 --force")"
log "sp064-nosampling" "$(sub -p short -c 8 --mem=16G -t 1:00:00 -J sp064-nosampling -o "$D/%x.%j.log" --wrap \
    "apptainer exec -B /scratch,/projects $IMGD/dorado-2.1.2-dna.sif modkit sample-probs $bam --no-sampling --percentiles 10 -t 8 --out-dir $D/nosampling --force")"

echo; echo "submitted; $(squeue --me -h | wc -l) jobs in the queue. Track with: squeue --me -o '%.10i %.40j %.8T %.20R'"
