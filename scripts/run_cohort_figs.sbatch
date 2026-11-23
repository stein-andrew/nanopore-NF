#!/bin/bash
#SBATCH --job-name=cohort_figs
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
# Orthogonal sets -> BED, optional per-sample read QC, then every cohort figure.
#   sbatch run_cohort_figs.sbatch          (the log lands in the directory you submit from)
# Edit the block below first. Everything is written under /projects, not /scratch, so the purge
# cannot take the stream-pass cache or the converted BEDs.
set -euo pipefail

# ---------------------------------------------------------------- edit these
SCRIPTS=/home/stein.an/ont-pipeline-v2/scripts
RESULTS=/scratch/stein.an/ont-pipeline/cohort-v1/results
SHEET=/scratch/stein.an/ont-pipeline/cohort-v1/analysis_samplesheet.csv
OUT=/projects/Genometechlab/andrew/cohort_figs
ORTHO=/projects/Genometechlab/andrew/orthogonal
GROUP=""                     # samplesheet column for figure 8, e.g. "condition"; empty = no figure 8
LEVELS="ctrl case"           # the two levels of GROUP, reference first
REPEATS=""                   # optional RepeatMasker BED (inosine-in-Alu panel)
RUN_EXTRACT=auto             # auto = run cohort_extract.py unless sample_qc/ already has every sample; 1 = always; 0 = never
STYLE=/projects/Genometechlab/matplotlib_stylesheets/styles/genometechlab_main.mplstyle
FONTS=/projects/Genometechlab/matplotlib_stylesheets/fonts
# Python with pandas, numpy, scipy, matplotlib and openpyxl, e.g.:
set +u; source /projects/Genometechlab/andrew/miniforge3/etc/profile.d/conda.sh; conda activate bioinf_tools; set -u
# module load samtools        # if RUN_EXTRACT=1 and samtools is not already on PATH
# -----------------------------------------------------------------------------

THREADS=${SLURM_CPUS_PER_TASK:-8}

# ---------------------------------------------------------------- preflight (fails in seconds, not after extraction)
# A conda env whose scipy was built with GCC 13 needs the env's own libstdc++ (CXXABI_1.3.15); if the system's
# /lib64/libstdc++.so.6 is loaded first (e.g. by a pip-installed wheel imported earlier), scipy fails to import.
PYENV=$(python -c 'import sys; print(sys.prefix)')
if [[ -e "$PYENV/lib/libstdc++.so.6" ]]; then
    export LD_LIBRARY_PATH="$PYENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
if ! python -c 'import numpy, pandas, matplotlib.pyplot, scipy.stats, scipy.cluster.hierarchy, openpyxl' 2>"$OUT.preflight.err"; then
    cat "$OUT.preflight.err" >&2
    echo "Python environment $PYENV cannot import the plotting stack. If the error names CXXABI/GLIBCXX:" >&2
    echo "  mamba install -n $(basename "$PYENV") -c conda-forge 'libstdcxx-ng>=13' 'libgcc-ng>=13'" >&2
    echo "  (or reinstall pip-installed matplotlib/scipy from conda-forge: mamba install -c conda-forge matplotlib scipy)" >&2
    exit 1
fi
rm -f "$OUT.preflight.err"
for d in "$RESULTS" "$OUT"; do
    [[ "$d" == /scratch/* ]] && echo "WARNING: $d is on /scratch (purged by mtime): the stream cache, sample_qc and figures can disappear"
done
mkdir -p "$OUT"
echo "[$(date +%T)] job $SLURM_JOB_ID on $(hostname); results $RESULTS; out $OUT"
for f in cohort_plots.py cohort_extract.py orthogonal_to_bed.py; do
    [[ -s "$SCRIPTS/$f" ]] || { echo "missing $SCRIPTS/$f" >&2; exit 1; }
done
[[ -d "$RESULTS" && -s "$SHEET" ]] || { echo "set RESULTS and SHEET at the top of this script" >&2; exit 1; }

# ---------------------------------------------------------------- 0. pyxlsb for the GLORI2025 .xlsb files
if ! python -c 'import pyxlsb' 2>/dev/null; then
    pip install --user --quiet pyxlsb \
        || echo "WARNING: pyxlsb not installed (no internet on this node?): the two GLORI2025 .xlsb sets will be skipped;" \
                "run 'pip install --user pyxlsb' on a login node and resubmit"
fi

# ---------------------------------------------------------------- 1. per-sample read QC (figure 11)
if [[ "$RUN_EXTRACT" == auto ]]; then
    n_samples=$(python -c "import pandas as pd; print(len(pd.read_csv('$SHEET', dtype=str)))")
    n_done=$(ls "$RESULTS"/sample_qc/*.integrity.tsv 2>/dev/null | wc -l)
    if (( n_done >= n_samples )); then
        echo "[$(date +%T)] sample_qc has $n_done of $n_samples samples: cohort_extract.py skipped (RUN_EXTRACT=1 to redo)"
        RUN_EXTRACT=0
    else
        RUN_EXTRACT=1
    fi
fi
if [[ "$RUN_EXTRACT" == 1 ]]; then
    if command -v samtools >/dev/null; then
        echo "[$(date +%T)] cohort_extract.py"
        python "$SCRIPTS/cohort_extract.py" all --results "$RESULTS" --samplesheet "$SHEET" --threads "$THREADS" \
            2>&1 | tee "$OUT/cohort_extract.log"
    else
        echo "WARNING: samtools not on PATH: skipping cohort_extract.py (figure 11 panels A-E stay empty)"
    fi
fi

# ---------------------------------------------------------------- 2. orthogonal sets -> BED
echo "[$(date +%T)] orthogonal_to_bed.py --inspect (what it detects; nothing written)"
python "$SCRIPTS/orthogonal_to_bed.py" batch --dir "$ORTHO" --out "$ORTHO/bed" --inspect 2>&1 | tee "$OUT/orthogonal_inspect.log"
echo "[$(date +%T)] orthogonal_to_bed.py"
python "$SCRIPTS/orthogonal_to_bed.py" batch --dir "$ORTHO" --out "$ORTHO/bed" 2>&1 | tee "$OUT/orthogonal_convert.log"
KS=$(sed -n 's/^--known-sites //p' "$OUT/orthogonal_convert.log" || true)
if [[ -z "$KS" ]]; then
    echo "WARNING: no orthogonal set written (see $OUT/orthogonal_convert.log and $ORTHO/bed/report.tsv): no figure 16"
fi
echo "sets not written or skipped:"
grep -E '^\S.*: skipped: |NOT WRITTEN' "$OUT/orthogonal_convert.log" || true

# ---------------------------------------------------------------- 3. figures
ARGS=(--results "$RESULTS" --samplesheet "$SHEET" --out "$OUT" --threads "$THREADS" --style "$STYLE" --font-dir "$FONTS")
[[ -n "$GROUP" ]] && ARGS+=(--group "$GROUP" --levels $LEVELS)
[[ -n "$KS" ]] && ARGS+=(--known-sites $KS)
[[ -n "$REPEATS" ]] && ARGS+=(--repeats "$REPEATS")
echo "[$(date +%T)] cohort_plots.py ${ARGS[*]}"
python "$SCRIPTS/cohort_plots.py" "${ARGS[@]}" 2>&1 | tee "$OUT/cohort_plots.log"

echo "[$(date +%T)] done. Check $OUT/notes.txt for 'style:' and 'font: Helvetica', and any FAILED line:"
grep -E '^(style|font)|FAILED' "$OUT/notes.txt" || true
