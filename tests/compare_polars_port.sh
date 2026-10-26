#!/usr/bin/env bash
# The polars rewrite of bin/ checked against the pandas versions it replaced (tests/reference),
# run on the same inputs: the synthetic edge cases in tests/data, then real chr22 inputs found
# in the pipeline's work directory. Every output must be identical (read_features: as a set of
# rows; its pandas version did not sort them).
#
#   bash tests/compare_polars_port.sh [WORK_DIR]          (RESULTS=<published results dir> for the matrices)
#
# Needs the mods and variants images, both with polars and pandas. One line per output, "same"
# or "DIFFERENT" with the first differing lines; exits 1 if anything differs or fails.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=${1:-/scratch/stein.an/ont-pipeline/work}
CDIR=${CONTAINER_DIR:-/scratch/stein.an/ont-pipeline/containers}
MODS=$CDIR/mods.sif
VARIANTS=$CDIR/variants.sif
DATA=$REPO/tests/data
RESULTS=${RESULTS:-/scratch/stein.an/ont-pipeline/test-chr22/results}
OUT=$(mktemp -d "${TMPDIR:-/tmp}/polars_port.XXXXXX")
BINDS=/projects/Genometechlab,/scratch
fail=0

py() {      # image ref|new script args...   (in the current directory)
    local img=$1 which=$2 script=$3; shift 3
    local dir=$REPO/bin file=$script.py
    if [[ $which == ref ]]; then dir=$REPO/tests/reference; file=${script}_pandas.py; fi
    apptainer exec --bind "$BINDS" "$img" env PYTHONPATH="$dir" python3 "$dir/$file" "$@"
}

both() {    # name image script args...   ref runs in $OUT/name/ref (plus $REF_EXTRA), new in $OUT/name/new
    local name=$1 img=$2 script=$3; shift 3
    mkdir -p "$OUT/$name/ref" "$OUT/$name/new"
    # shellcheck disable=SC2086
    (cd "$OUT/$name/ref" && py "$img" ref "$script" "$@" ${REF_EXTRA:-} > stdout.txt 2> stderr.txt) \
        || { echo "  REF FAILED $name:"; tail -3 "$OUT/$name/ref/stderr.txt" | sed 's/^/      /'; fail=1; }
    (cd "$OUT/$name/new" && py "$img" new "$script" "$@" > stdout.txt 2> stderr.txt) \
        || { echo "  NEW FAILED $name:"; tail -5 "$OUT/$name/new/stderr.txt" | sed 's/^/      /'; fail=1; }
}

same() {    # label name file [sorted]   compares $OUT/name/ref/file with $OUT/name/new/file
    local label=$1 a=$OUT/$2/ref/$3 b=$OUT/$2/new/$3 how=${4:-}
    if [[ ! -e $a || ! -e $b ]]; then echo "  MISSING    $label ($3)"; fail=1; return 0; fi
    if [[ $how == sorted ]]; then
        if cmp -s <(zcat -f "$a" | sort) <(zcat -f "$b" | sort); then echo "  same       $label: $3"; return 0; fi
    elif cmp -s <(zcat -f "$a") <(zcat -f "$b"); then echo "  same       $label: $3"; return 0; fi
    echo "  DIFFERENT  $label: $3"
    diff <(zcat -f "$a") <(zcat -f "$b") | head -8 | sed 's/^/      /' || true    # diff exits 1 on a difference
    fail=1
}

find_dir() {  # newest work dir holding all the given files (the first one is searched for)
    local first=$1; shift
    find "$WORK" -mindepth 3 -maxdepth 3 -name "$first" -printf '%T@ %h\n' 2>/dev/null | sort -rn | while read -r _ d; do
        local ok=1 g
        for g in "$@"; do [[ -e $d/$g ]] || ok=0; done
        if (( ok )); then echo "$d"; break; fi
    done
}

one() { compgen -G "$1" | head -n1; }   # first file matching a glob

echo "outputs in $OUT"

echo "== fp_correct (edge cases, every flag combination)"
MT="a,0.7 m,0.85 17802,0.949 19229,0.99"
i=0
for flags in "" "--valid_only" "--filter_kmer" "--filter_mismatch" "-n 25"; do
    i=$((i + 1))
    # shellcheck disable=SC2086
    REF_EXTRA="--engine pandas" both "fp$i" "$MODS" fp_correct -m "$DATA/edge.bed" -ref "$DATA/ref2w.fa" \
        -e "$DATA/errtab2.tsv.gz" -mt $MT -o out.bed $flags
    same "fp_correct [${flags:-default}]" "fp$i" out.bed
done

echo "== known_novel (edge cases)"
both kn "$VARIANTS" known_novel --calls "$DATA/kn_calls.tsv" --gnomad "$DATA/kn_gnomad.tsv" --mito "$DATA/kn_mito.tsv" \
    --clinvar "$DATA/kn_clinvar.tsv" --prefix S
for f in S.variants.tsv.gz S.variants.summary.tsv S.annot.tsv S.annot.hdr; do same known_novel kn "$f"; done

echo "== feature_mods (edge cases)"
both fmp "$MODS" feature_mods positions --candidates "$DATA/fm_candidates.bed" --candidate-pct 20
same "feature_mods positions" fmp stdout.txt
both fm "$MODS" feature_mods levels --counts "$DATA/fm_counts.tsv" --candidates "$DATA/fm_candidates.bed" --candidate-pct 20 \
    --catalogue "$DATA/catalogue.gtf" --classification "$DATA/classification.txt.gz" --reference "$DATA/ref2w.fa" \
    --errortable "$DATA/errtab2.tsv.gz" --givt-flags "$DATA/givt_flags.txt" --min-cov 20 --min-pct 20 --prefix S
for f in S.isoform_mods.bed S.gene_mods.bed S.feature_mods.tsv.gz S.feature_mods.summary.tsv; do same feature_mods fm "$f"; done

echo "== polya_features (edge cases, without and with caps)"
both pa "$MODS" polya_features --reads "$DATA/pa_reads_pt.tsv" --features "$DATA/pa_read_features.tsv" \
    --catalogue "$DATA/catalogue.gtf" --classification "$DATA/classification.txt.gz" --prefix P
both pac "$MODS" polya_features --reads "$DATA/pa_reads_pt.tsv" --features "$DATA/pa_read_features.tsv" \
    --catalogue "$DATA/catalogue.gtf" --classification "$DATA/classification.txt.gz" --cap-nuclear 150 --cap-mito 250 --prefix P
for n in pa pac; do
    for f in P.polya_summary.tsv P.polya_hist.tsv P.polya_reads.tsv.gz P.polya_features.tsv.gz P.isoform_polya.bed P.gene_polya.bed; do
        same "polya_features ($n)" "$n" "$f"
    done
done

echo "== read_features (edge cases)"
both rf "$MODS" read_features "$DATA/rf_read_assignments.tsv.gz"
same read_features rf stdout.txt sorted

echo "== known_novel (edge cases, with the IG/TR/HLA column)"
both knm "$VARIANTS" known_novel --calls "$DATA/kn_calls_masked.tsv" --gnomad "$DATA/kn_gnomad.tsv" --mito "$DATA/kn_mito.tsv" \
    --clinvar "$DATA/kn_clinvar.tsv" --prefix S
for f in S.variants.tsv.gz S.variants.summary.tsv S.annot.tsv S.annot.hdr; do same "known_novel (masked)" knm "$f"; done

echo "== cohort_matrices (edge cases)"
both cmg "$MODS" cohort_matrices genes --counts "A=$DATA/cm_genes_A.tsv" --counts "B=$DATA/cm_genes_B.tsv" \
    --catalogue "$DATA/catalogue.gtf" --classification "$DATA/classification.txt.gz" --prefix C
for f in C.gene_counts.tsv.gz C.gene_cpm.tsv.gz; do same "cohort_matrices genes" cmg "$f"; done
both cmm "$MODS" cohort_matrices mods --sites "$DATA/cm_sites.tsv" --long "$DATA/cm_long.tsv" --samples A B --min-pct 20 --prefix C
for f in C.mods.n_valid.tsv.gz C.mods.n_mod.tsv.gz C.mods.pct_fpcorr.tsv.gz C.mods.summary.tsv; do same "cohort_matrices mods" cmm "$f"; done
both cmp "$MODS" cohort_matrices polya --features "A=$DATA/cm_polya_A.tsv.gz" --features "B=$DATA/cm_polya_B.tsv.gz" --prefix C
for f in C.polya.{isoform,gene}.{median_pt,n_tails}.tsv.gz; do same "cohort_matrices polya" cmp "$f"; done

echo "== cohort_matrices on the published chr22 results ($RESULTS)"
if [[ -d $RESULTS/isoforms/GRCh38/requant ]]; then
    ids=(); gargs=(); fpairs=(); pargs=()
    for d in "$RESULTS"/isoforms/GRCh38/requant/*/; do
        s=$(basename "$d"); ids+=("$s")
        gargs+=(--counts "$s=$d$s.gene_counts.tsv")
        fpairs+=("$s=$RESULTS/mods/GRCh38/$s/$s.cov20_fpcorr.bed")
        pargs+=(--features "$s=$RESULTS/polya/GRCh38/$s/$s.polya_features.tsv.gz")
    done
    C=$RESULTS/isoforms/GRCh38/cohort/cohort.GRCh38
    both cmg_real "$MODS" cohort_matrices genes "${gargs[@]}" --catalogue "$C.catalogue.gtf" \
        --classification "$C.classify.classification.txt.gz" --prefix C
    for f in C.gene_counts.tsv.gz C.gene_cpm.tsv.gz; do same "cohort_matrices genes (real)" cmg_real "$f"; done
    mkdir -p "$OUT/cm_real_in"
    bash "$REPO/bin/mods_site_values.sh" 20 "$OUT/cm_real_in/sites.tsv" "$OUT/cm_real_in/long.tsv" "${fpairs[@]}"
    echo "  (${#ids[@]} samples; $(wc -l < "$OUT/cm_real_in/sites.tsv") final sites in any sample)"
    both cmm_real "$MODS" cohort_matrices mods --sites "$OUT/cm_real_in/sites.tsv" --long "$OUT/cm_real_in/long.tsv" \
        --samples "${ids[@]}" --min-pct 20 --prefix C
    for f in C.mods.n_valid.tsv.gz C.mods.n_mod.tsv.gz C.mods.pct_fpcorr.tsv.gz C.mods.summary.tsv; do same "cohort_matrices mods (real)" cmm_real "$f"; done
    both cmp_real "$MODS" cohort_matrices polya "${pargs[@]}" --prefix C
    for f in C.polya.{isoform,gene}.{median_pt,n_tails}.tsv.gz; do same "cohort_matrices polya (real)" cmp_real "$f"; done
else echo "  skipped    no published results at $RESULTS"; fi

echo "== real chr22 inputs from $WORK (ST1918)"
d=$(find_dir ST1918.cov20.bed ST1918.givt_flags.txt || true)
if [[ -n $d ]]; then
    # shellcheck disable=SC2046
    REF_EXTRA="--engine pandas" both fp_real "$MODS" fp_correct -m "$d/ST1918.cov20.bed" -ref "$(one "$d/*.fa")" \
        -e "$(one "$d/*9mer*.tsv.gz")" -mt $(cat "$d/ST1918.givt_flags.txt") -o out.bed
    same "fp_correct (real)" fp_real out.bed
    if cmp -s "$OUT/fp_real/new/out.bed" "$d/ST1918.cov20_fpcorr.bed"; then echo "  same       fp_correct (real): new = the pipeline's own output"
    else echo "  DIFFERENT  fp_correct (real): new vs the pipeline's ST1918.cov20_fpcorr.bed"; fail=1; fi
else echo "  skipped    fp_correct: no MOD_STAGES work dir for ST1918"; fi

d=$(find_dir calls.tsv gnomad.tsv mito.tsv clinvar.tsv ST1918.clair3_rna.vcf.gz || true)
if [[ -n $d ]]; then
    both kn_real "$VARIANTS" known_novel --calls "$d/calls.tsv" --gnomad "$d/gnomad.tsv" --mito "$d/mito.tsv" \
        --clinvar "$d/clinvar.tsv" --prefix S
    for f in S.variants.tsv.gz S.variants.summary.tsv S.annot.tsv S.annot.hdr; do same "known_novel (real)" kn_real "$f"; done
else echo "  skipped    known_novel: no KNOWN_NOVEL work dir for ST1918"; fi

d=$(find_dir ST1918.counts.tsv ST1918.cov20_fpcorr.bed ST1918.givt_flags.txt || true)
if [[ -n $d ]]; then
    both fm_real "$MODS" feature_mods levels --counts "$d/ST1918.counts.tsv" --candidates "$d/ST1918.cov20_fpcorr.bed" \
        --candidate-pct 20 --catalogue "$(one "$d/*.catalogue.gtf")" --classification "$(one "$d/*.classification.txt.gz")" \
        --reference "$(one "$d/*.fa")" --errortable "$(one "$d/*9mer*.tsv.gz")" --givt-flags "$d/ST1918.givt_flags.txt" \
        --min-cov 20 --min-pct 20 --prefix S
    for f in S.isoform_mods.bed S.gene_mods.bed S.feature_mods.tsv.gz S.feature_mods.summary.tsv; do same "feature_mods (real)" fm_real "$f"; done
else echo "  skipped    feature_mods: no MOD_FEATURES work dir for ST1918"; fi

d=$(find_dir ST1918.reads_pt.tsv ST1918.read_features.tsv ST1918.read_assignments.tsv.gz || true)
if [[ -n $d ]]; then
    both pa_real "$MODS" polya_features --reads "$d/ST1918.reads_pt.tsv" --features "$d/ST1918.read_features.tsv" \
        --catalogue "$(one "$d/*.catalogue.gtf")" --classification "$(one "$d/*.classification.txt.gz")" --prefix P
    for f in P.polya_summary.tsv P.polya_hist.tsv P.polya_reads.tsv.gz P.polya_features.tsv.gz P.isoform_polya.bed P.gene_polya.bed; do
        same "polya_features (real)" pa_real "$f"
    done
    both rf_real "$MODS" read_features "$d/ST1918.read_assignments.tsv.gz"
    same "read_features (real)" rf_real stdout.txt sorted
else echo "  skipped    polya_features, read_features: no POLYA_FEATURES work dir for ST1918"; fi

if (( fail )); then echo "== DIFFERENCES or FAILURES above (outputs kept in $OUT)"; exit 1; fi
echo "== every output identical"
