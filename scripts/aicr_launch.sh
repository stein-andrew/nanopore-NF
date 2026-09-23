#!/usr/bin/env bash
# Run a pipeline on AICR. The login node kills the JVM within a second, so the
# Nextflow head job itself is submitted to the cpu partition.
#
#   scripts/aicr_launch.sh <launch_dir> basecall.nf -profile aicr,dna --input runs.csv --outdir results/run1
#
# <launch_dir> holds .nextflow/ and the head log; re-run the same command with
# -resume from the same launch_dir to continue. The head job's 24 h limit is the
# partition maximum; if a run outlives it, resubmit with -resume.
set -euo pipefail

LAUNCH=${1:?usage: $0 <launch_dir> <script.nf> [nextflow args...]}
shift
(( $# )) || { echo "usage: $0 <launch_dir> <script.nf> [nextflow args...]" >&2; exit 2; }

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
JAVA=${JAVA_HOME:-$HOME/.sdkman/candidates/java/current}
NEXTFLOW=${NEXTFLOW:-$(command -v nextflow || echo "$HOME/.local/bin/nextflow")}
[[ -x $JAVA/bin/java ]] || { echo "no java at $JAVA/bin/java -- set JAVA_HOME" >&2; exit 1; }
[[ -x $NEXTFLOW ]]      || { echo "no nextflow at $NEXTFLOW -- set NEXTFLOW" >&2; exit 1; }

script=$1; shift
[[ $script == /* ]] || script=$REPO/$script

mkdir -p "$LAUNCH"
cmd=$(printf '%q ' "$NEXTFLOW" run "$script" "$@")
sbatch -p cpu -c 2 --mem 8G -t 24:00:00 -J nf-head -D "$LAUNCH" -o "$LAUNCH/nf-head.%j.log" \
    --wrap "export JAVA_HOME=$JAVA PATH=$JAVA/bin:\$PATH; $cmd"
echo "log: $LAUNCH/nf-head.<jobid>.log"
