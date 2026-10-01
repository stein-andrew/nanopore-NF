#!/usr/bin/env bash
# Push a Stage 1 outdir to another cluster with rsync, then verify it there.
#
#   scripts/transfer.sh <stage1_outdir> <user@host:/dest/dir> [--with-ubam]
#   scripts/transfer.sh results/run1 stein.an@<explorer-login-host>:/scratch/stein.an/ont-pipeline/results/run1
#
# The uBAMs are left behind unless --with-ubam is given: Stage 2 only reads the
# aligned BAMs, which already carry MM/ML (and mv/pt for RNA).
#
# Verification is `md5sum -c checksums.md5` at the destination, over the same
# relative paths analysis_samplesheet.csv uses.
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

SRC=${1:?usage: $0 <stage1_outdir> <user@host:/dest/dir> [--with-ubam]}
DEST=${2:?usage: $0 <stage1_outdir> <user@host:/dest/dir> [--with-ubam]}
WITH_UBAM=${3:-}

SRC=${SRC%/}
[[ -s $SRC/analysis_samplesheet.csv ]] || die "$SRC has no analysis_samplesheet.csv -- is it a finished Stage 1 outdir?"
[[ -s $SRC/checksums.md5 ]]            || die "$SRC has no checksums.md5"

EXCLUDE=(--exclude 'ubam/')
[[ $WITH_UBAM == --with-ubam ]] && EXCLUDE=()

if [[ $DEST == *:* ]]; then
    HOST=${DEST%%:*}
    DIR=${DEST#*:}
    ssh "$HOST" "mkdir -p '$DIR'"
else
    HOST=
    DIR=$DEST
    mkdir -p "$DIR"
fi

# -a keeps times, so a re-run skips files that already arrived intact.
# shellcheck disable=SC2086
rsync -a --partial --info=progress2 ${RSYNC_OPTS:-} "${EXCLUDE[@]}" "$SRC/" "$DEST/"

echo "verifying checksums at $DEST"
# Explicit if: set -e does not stop on a failure on the left of `&&`.
if [[ -n $HOST ]]; then
    verify=(ssh "$HOST" "cd '$DIR' && md5sum -c --quiet checksums.md5")
else
    verify=(bash -c "cd '$DIR' && md5sum -c --quiet checksums.md5")
fi
if ! "${verify[@]}"; then
    die "checksum mismatch at $DEST -- resend with: RSYNC_OPTS=--checksum $0 $*"
fi
echo "all BAMs verified"
echo "Stage 2 input: $DIR/analysis_samplesheet.csv"
