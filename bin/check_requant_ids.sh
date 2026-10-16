#!/usr/bin/env bash
# Fail unless every transcript IsoQuant counted is one of the catalogue's transcript_ids.
#
#   check_requant_ids.sh <catalogue.gtf> <S.transcript_counts.tsv>
#
# Catches IsoQuant loading the catalogue wrongly: in the first chr22 test, gffutils gave
# every transcript a placeholder name (transcript_1, transcript_10, ...) and attached none
# of them to genes, so every read came out intergenic and the counts were all zero. That
# failed only later, in COHORT_FILTER; this stops the requantification task itself.
set -euo pipefail

catalogue=${1:?usage: $0 <catalogue.gtf> <transcript_counts.tsv>}
counts=${2:?usage: $0 <catalogue.gtf> <transcript_counts.tsv>}

# POSIX awk only: the image's awk is Ubuntu's mawk.
awk -F'\t' '
    NR == FNR {
        if ($3 == "transcript" && match($9, /transcript_id "[^"]+"/))
            cat[substr($9, RSTART + 15, RLENGTH - 16)] = 1
        next
    }
    FNR > 1 && $1 !~ /^__/ { n++; if ($1 in cat) ok++ }
    END {
        if (n == 0 || ok < n) {
            printf "IsoQuant counted %d transcripts, %d of them catalogue IDs: it did not load the catalogue correctly\n", n, ok > "/dev/stderr"
            exit 1
        }
        printf "requant check: all %d counted transcripts are catalogue IDs\n", n
    }' "$catalogue" "$counts"
