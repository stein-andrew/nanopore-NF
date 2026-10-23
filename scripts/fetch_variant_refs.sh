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
#   scripts/fetch_variant_refs.sh adopt                       # seconds: finish off old-style partials
#
# Array tasks: 0-23  gnomAD v4.1 joint (exomes + genomes), chr1-chr22, chrX, chrY
#              24    gnomAD v3.1 chrM (v4.1 has no chrM; mito frequencies are v3.1's)
#              25    ClinVar GRCh38, contigs renamed 1 -> chr1 ... MT -> chrM, all INFO kept. Its
#                    header declares no contigs, which bcftools needs to write BCF (and to
#                    rename), so the rename is done as text and the renamed contigs declared.
# finish: concatenates the 24 nuclear files into one indexed BCF, checks that no records were
# lost, removes the per-chromosome files and writes MANIFEST.tsv (sources, fields, counts).
# adopt: the first version named its temporary output <out>.partial; htslib finds an index by
# the data file's extension, so `bcftools index -n` on it failed after a complete stream.
# adopt renames such files to <out minus .bcf>.partial.bcf, counts them and finalises them,
# so nothing is fetched again. Temporary files are now named that way from the start.
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

clinvar_to_chr() {      # ClinVar VCF text on stdin -> chr-named VCF text with ##contig lines
    # Records on contigs other than 1-22, X, Y, MT are dropped and counted on stderr.
    awk -F'\t' -v OFS='\t' '
        BEGIN { for (i = 1; i <= 22; i++) m[i] = "chr" i; m["X"] = "chrX"; m["Y"] = "chrY"; m["MT"] = "chrM" }
        /^##contig=/ { next }
        /^##/        { print; next }
        /^#CHROM/    { for (i = 1; i <= 22; i++) print "##contig=<ID=chr" i ">"
                       print "##contig=<ID=chrX>"; print "##contig=<ID=chrY>"; print "##contig=<ID=chrM>"
                       print; next }
        ($1 in m)    { $1 = m[$1]; print; kept++; next }
                     { dropped[$1]++ }
        END          { printf "clinvar_to_chr: %d records kept\n", kept > "/dev/stderr"
                       for (c in dropped) printf "clinvar_to_chr: dropped %d records on contig %s\n", dropped[c], c > "/dev/stderr" }'
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
            echo "   all INFO kept; contigs renamed as text: 1>chr1 ... X>chrX Y>chrY MT>chrM (header declares none)"
        fi
    done
    echo "== output: $OUT  (writable: $([[ -w $OUT || -w $(dirname "$OUT") ]] && echo yes || echo NO))"
}

cmd_fetch() {
    local i=${SLURM_ARRAY_TASK_ID:?run as an array task: sbatch --array=0-25%8 $0 fetch}
    [[ -s $SIF ]] || die "no image at $SIF"
    local url out part
    url=$(source_url "$i"); out=$(out_path "$i"); part="${out%.bcf}.partial.bcf"
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
        curl -sfL --retry 5 --retry-delay 30 "$url" | gzip -dc | clinvar_to_chr \
          | bcf view --threads 2 -Ob -o "$part" -
        keep='all'
    fi
    bcf index --threads 2 -f "$part"
    finalise "$i" "$part" "$out" "$url" "$keep"
}

finalise() {            # task part out url kept-fields: count, move into place, manifest line
    local i=$1 part=$2 out=$3 url=$4 keep=$5 n
    n=$(bcf index -n "$part") || die "cannot count records in $part (index $part.csi)"
    (( n > 0 )) || die "$part has no records"
    mv "$part.csi" "$out.csi"; mv "$part" "$out"
    mkdir -p "$OUT/.manifest"
    printf '%s\t%s\t%s\t%s\t%s\n' "$(basename "$out")" "$url" "$n" "$keep" "$(date -u +%FT%TZ)" > "$OUT/.manifest/$i.tsv"
    write_manifest
    echo "done: $(basename "$out"): $n records"
}

write_manifest() {      # every finished file's line, regenerated after each one
    { printf 'file\tsource\trecords\tinfo_kept\tfinished\n'; cat "$OUT"/.manifest/*.tsv 2>/dev/null | sort -V; } > "$OUT/MANIFEST.tsv"
}

cmd_adopt() {
    [[ -s $SIF ]] || die "no image at $SIF"
    local i out old part url keep missing=()
    for i in $(seq 0 24); do
        out=$(out_path "$i"); old="$out.partial"; part="${out%.bcf}.partial.bcf"
        if [[ -s $out && -s $out.csi ]]; then echo "already final: $(basename "$out")"; continue; fi
        if [[ -s $old && -s $old.csi ]]; then mv "$old" "$part"; mv "$old.csi" "$part.csi"; fi
        if [[ ! ( -s $part && -s $part.csi ) ]]; then missing+=("$i"); echo "NOT FETCHED: $(basename "$out") (task $i)"; continue; fi
        url=$(source_url "$i"); keep=$(bcf view -h "$part" | keep_fields)
        finalise "$i" "$part" "$out" "$url" "$keep"
    done
    (( ${#missing[@]} == 0 )) || echo "re-fetch: sbatch --array=$(IFS=,; echo "${missing[*]}") $0 fetch"
}

cmd_finish() {
    local files=() missing=() total=0 n
    for c in "${CHROMS[@]}"; do
        if [[ -s $OUT/$NUCLEAR_DIR/$c.bcf.csi ]]; then files+=("$OUT/$NUCLEAR_DIR/$c.bcf"); else missing+=("$c"); fi
    done
    (( ${#missing[@]} == 0 )) || die "not fetched yet: ${missing[*]}"
    for f in "${files[@]}"; do n=$(bcf index -n "$f"); total=$(( total + n )); done

    local part="$OUT/${NUCLEAR_BCF%.bcf}.partial.bcf"
    bcf concat --threads "${SLURM_CPUS_PER_TASK:-4}" -Ob -o "$part" "${files[@]}"
    bcf index --threads "${SLURM_CPUS_PER_TASK:-4}" -f "$part"
    n=$(bcf index -n "$part")
    (( n == total )) || die "concatenated $n records but the chromosomes hold $total: keeping the per-chromosome files"
    mv "$part.csi" "$OUT/$NUCLEAR_BCF.csi"; mv "$part" "$OUT/$NUCLEAR_BCF"
    rm -rf "${OUT:?}/$NUCLEAR_DIR"

    printf '%s\t%s\t%s\t%s\t%s\n' "$NUCLEAR_BCF" "concat of the 24 gnomAD v4.1 joint chromosome files above" \
        "$n" "as above" "$(date -u +%FT%TZ)" > "$OUT/.manifest/nuclear.tsv"
    write_manifest
    echo "done: $OUT/$NUCLEAR_BCF ($n records); see $OUT/MANIFEST.tsv"
}

main() {
    case ${1:-} in
        check)  cmd_check ;;
        fetch)  cmd_fetch ;;
        finish) cmd_finish ;;
        adopt)  cmd_adopt ;;
        *)      die "usage: $0 {check|fetch|finish|adopt}  (see the header of this script)" ;;
    esac
}

# run when executed (or by sbatch), not when sourced for testing
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then main "$@"; fi
