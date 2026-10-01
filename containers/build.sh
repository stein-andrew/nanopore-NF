#!/usr/bin/env bash
# Build one or more container images into $CONTAINER_DIR.
#
#   containers/build.sh dorado-2.1.2-dna dorado-1.3.0-rna
#   sbatch -p short -c 8 --mem 32G -t 4:00:00 containers/build.sh dorado-2.1.2-dna
#
# Submit from the repo root: sbatch runs a spooled copy of this script, so under
# SLURM the .def files are found via $SLURM_SUBMIT_DIR/containers.

set -euo pipefail

HERE=${SLURM_SUBMIT_DIR:+$SLURM_SUBMIT_DIR/containers}
[[ -n $HERE && -d $HERE ]] || HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CONTAINER_DIR=${CONTAINER_DIR:-/scratch/$USER/ont-pipeline/containers}
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-/scratch/$USER/tmp/apptainer}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/scratch/$USER/.apptainer_cache}
mkdir -p "$CONTAINER_DIR" "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

(( $# )) || { echo "usage: $0 <name> [name...]   (names of containers/*.def without .def)" >&2; exit 2; }

for name in "$@"; do
    def="$HERE/$name.def"
    sif="$CONTAINER_DIR/$name.sif"
    [[ -f $def ]] || { echo "no such definition: $def" >&2; exit 2; }
    if [[ -s $sif ]]; then
        echo "exists, skipping: $sif"
        continue
    fi
    echo "building $sif"
    # Serialize builds that share a cache: two concurrent pulls of the same base
    # image corrupt it ("conveyor failed to get ... unexpected end of JSON input").
    flock "$APPTAINER_CACHEDIR/.build.lock" apptainer build --fakeroot "$sif.partial" "$def"
    mv "$sif.partial" "$sif"
    apptainer test "$sif"
done
