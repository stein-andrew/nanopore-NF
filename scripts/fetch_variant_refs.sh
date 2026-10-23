#!/usr/bin/env bash
#SBATCH --job-name=varrefs
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=varrefs_%A_%a.out
#
# Genome-wide variant references for calling Clair3-RNA variants known or novel (step 6).
# Whole genome, every variant; sites only (no per-person genotypes) and stripped to the INFO
# fields the comparison uses. Everything is streamed, so nothing raw (~1 TB) is ever stored.
#
#   scripts/fetch_variant_refs.sh check                       # seconds: what each source keeps
#   sbatch --array=0-25%8 scripts/fetch_variant_refs.sh fetch
#   sbatch --dependency=afterok:<array job id> scripts/fetch_variant_refs.sh finish
#
# Array tasks: 0-23  gnomAD v4.1 joint (exomes + genomes), chr1-chr22, chrX, chrY
#              24    gnomAD v3.1 chrM (v4.1 has no chrM; mito frequencies are v3.1's)
#              25    ClinVar GRCh38, contigs renamed 1 -> chr1 ... MT -> chrM, all INFO kept
# finish: concatenates the 24 nuclear files into one indexed BCF, checks that no records were
# lost, removes the per-chromosome files and writes MANIFEST.tsv (sources, fields, counts).
#
# Kept INFO fields are matched on each file's own header (KEEP_RE), never guessed: a field
# name that isn't there would make `bcftools annotate -x ^...` silently keep nothing. A
# gnomAD file with no AF-type field left after matching is an error.
#
# Environment: VARREF_DIR (output, default below), VARREF_SIF (the variants image).
set -euo pipefail

OUT=${VARREF_DIR:-/projects/Genometechlab/referenceFastaFiles/variants}
SIF=${VARREF_SIF:-/scratch/stein.an/ont-pipeline/containers/variants.sif}
GNOMAD=https://storage.googleapis.com/gcp-public-data--gnomad/release
CLINVAR_URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
CHROMS=(chr{1..22} chrX chrY)
NUCLEAR_DIR=gnomad_v4.1_joint_by_chrom
NUCLEAR_BCF=gnomad.joint.v4.1.sites.min.bcf
MITO_BCF=gnomad.genomes.v3.1.sites.chrM.min.bcf
CLINVAR_BCF=clinvar.GRCh38.chr.bcf
KEEP_RE='^((AC|AN|AF|nhomalt)(_(joint|exomes|genomes))?|AC_(hom|het)|AF_(hom|het)|(exomes|genomes|joint)_filters)$'

die() { echo "fetch_variant_refs.sh: $*" >&2; exit 1; }
bcf() { apptainer exec --bind /projects/Genometechlab,/scratch "$SIF" bcftools "$@"; }

source_url() {          # task -> URL
    local i=$1
    if (( i >= 0 && i < 24 )); then echo "$GNOMAD/4.1/vcf/joint/gnomad.joint.v4.1.sites.${CHROMS[$i]}.vcf.bgz"
    elif (( i == 24 )); then echo "$GNOMAD/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
    elif (( i == 25 )); then echo "$CLINVAR_URL"
    else die "no task $i (tasks are 0-25)"; fi
}

out_path() {            # task -> output BCF
    local i=$1
    if (( i < 24 )); then echo "$OUT/$NUCLEAR_DIR/${CHROMS[$i]}.bcf"
    elif (( i == 24 )); then echo "$OUT/$MITO_BCF"
    else echo "$OUT/$CLINVAR_BCF"; fi
}

keep_fields() {         # VCF header on stdin -> INFO/A,INFO/B,... of the IDs matching KEEP_RE
    # empty output (not a failure) when nothing matches, so has_af can say so; with pipefail a
    # bare grep miss would end the script without a message
    sed -n 's/^##INFO=<ID=\([^,>]*\).*/\1/p' | { grep -E "$KEEP_RE" || true; } | sed 's/^/INFO\//' | paste -sd, -
}

has_af() {              # kept-field list -> true if any AF-type field survived
    tr ',' '\n' <<< "$1" | grep -qE '^INFO/AF'
}

rename_map() {          # ClinVar contig names -> the reference's chr names
    for c in {1..22} X Y; do printf '%s\tchr%s\n' "$c" "$c"; done
    printf 'MT\tchrM\n'
}

cmd_check() {
    [[ -s $SIF ]] || die "no image at $SIF (containers/build.sh variants)"
    for i in 21 24 25; do
        local url hdr keep
        url=$(source_url "$i")
        hdr=$(bcf view -h "$url") || die "cannot read the header of $url"
        echo "== $url"
        if (( i < 25 )); then
            keep=$(keep_fields <<< "$hdr")
            has_af "$keep" || die "no AF-type INFO field matched in $url; INFO IDs there: $(sed -n 's/^##INFO=<ID=\([^,>]*\).*/\1/p' <<< "$hdr" | head -40 | paste -sd' ' -)"
            echo "   INFO fields in file: $(grep -c '^##INFO=' <<< "$hdr")    kept: $keep"
            echo "   FILTER values: $(sed -n 's/^##FILTER=<ID=\([^,>]*\).*/\1/p' <<< "$hdr" | paste -sd' ' -)"
        else
            echo "   release: $(sed -n 's/^##fileDate=//p' <<< "$hdr")    contigs, first 5: $(sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' <<< "$hdr" | head -5 | paste -sd' ' -)"
            echo "   all INFO kept; contigs renamed: $(rename_map | head -2 | tr '\t\n' '> ') ... MT>chrM"
        fi
    done
    echo "== output: $OUT  (writable: $([[ -w $OUT || -w $(dirname "$OUT") ]] && echo yes || echo NO))"
}

cmd_fetch() {
    local i=${SLURM_ARRAY_TASK_ID:?run as an array task: sbatch --array=0-25%8 $0 fetch}
    [[ -s $SIF ]] || die "no image at $SIF"
    local url out part
    url=$(source_url "$i"); out=$(out_path "$i"); part="$out.partial"
    mkdir -p "$(dirname "$out")" "$OUT/.manifest"
    if [[ -s $out && -s $out.csi ]]; then echo "already done: $out"; return 0; fi

    echo "task $i: $url -> $out"
    if (( i < 25 )); then
        local keep
        keep=$(bcf view -h "$url" | keep_fields)
        has_af "$keep" || die "no AF-type INFO field matched in $url"
        echo "keeping: $keep"
        curl -sfL --retry 5 --retry-delay 30 "$url" \
          | bcf annotate --threads 2 -x "^$keep" -Ob -o "$part" -
    else
        local map; map=$(mktemp); rename_map > "$map"
        curl -sfL --retry 5 --retry-delay 30 "$url" \
          | bcf annotate --threads 2 --rename-chrs "$map" -Ob -o "$part" -
        rm -f "$map"
        keep='all'
    fi
    bcf index --threads 2 -f "$part"
    local n; n=$(bcf index -n "$part")
    (( n > 0 )) || die "$part has no records"
    mv "$part.csi" "$out.csi"; mv "$part" "$out"
    printf '%s\t%s\t%s\t%s\t%s\n' "$(basename "$out")" "$url" "$n" "$keep" "$(date -u +%FT%TZ)" > "$OUT/.manifest/$i.tsv"
    echo "done: $n records"
}

cmd_finish() {
    local files=() missing=() total=0 n
    for c in "${CHROMS[@]}"; do
        if [[ -s $OUT/$NUCLEAR_DIR/$c.bcf.csi ]]; then files+=("$OUT/$NUCLEAR_DIR/$c.bcf"); else missing+=("$c"); fi
    done
    (( ${#missing[@]} == 0 )) || die "not fetched yet: ${missing[*]}"
    for f in "${files[@]}"; do n=$(bcf index -n "$f"); total=$(( total + n )); done

    local part="$OUT/$NUCLEAR_BCF.partial"
    bcf concat --threads "${SLURM_CPUS_PER_TASK:-4}" -Ob -o "$part" "${files[@]}"
    bcf index --threads "${SLURM_CPUS_PER_TASK:-4}" -f "$part"
    n=$(bcf index -n "$part")
    (( n == total )) || die "concatenated $n records but the chromosomes hold $total: keeping the per-chromosome files"
    mv "$part.csi" "$OUT/$NUCLEAR_BCF.csi"; mv "$part" "$OUT/$NUCLEAR_BCF"
    rm -rf "${OUT:?}/$NUCLEAR_DIR"

    {
        printf 'file\tsource\trecords\tinfo_kept\tfetched\n'
        cat "$OUT"/.manifest/*.tsv | sort -V
        printf '%s\t%s\t%s\t%s\t%s\n' "$NUCLEAR_BCF" "concat of the 24 gnomAD v4.1 joint files" "$n" "as above" "$(date -u +%FT%TZ)"
    } > "$OUT/MANIFEST.tsv"
    echo "done: $OUT/$NUCLEAR_BCF ($n records); see $OUT/MANIFEST.tsv"
}

main() {
    case ${1:-} in
        check)  cmd_check ;;
        fetch)  cmd_fetch ;;
        finish) cmd_finish ;;
        *)      die "usage: $0 {check|fetch|finish}  (see the header of this script)" ;;
    esac
}

# run when executed (or by sbatch), not when sourced for testing
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then main "$@"; fi
