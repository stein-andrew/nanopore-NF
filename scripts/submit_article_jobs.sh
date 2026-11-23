#!/usr/bin/env bash
# submit_article_jobs.sh: reruns with the expanded truth sets (no basecalling): orthogonal scoring, all-against-all
# truth-set concordance (exact and +-1 nt) and m6A class support, for both calling rules. Run with bioinf_tools active:
#     mamba activate bioinf_tools && bash ~/ont-pipeline-v2/scripts/submit_article_jobs.sh
set -euo pipefail
SC=~/ont-pipeline-v2/scripts; B=/scratch/stein.an/basecall_cmp; PY=$(command -v python)
FA=/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
MASK=/scratch/stein.an/ont-pipeline/cohort-v1/psi_loci_check/clusters.bed; ORTHO=/projects/Genometechlab/andrew/orthogonal
HEK="08_01_25_RNA004_HEK293 09_09_25_RNA004_HEK293 10_06_25_RNA004_HEK293_CP"
declare -A ARMDIR=([nofilter]=$B/pileup064_nofilter [lab050]=$B/pileup064_lab050)

"$PY" -c 'import importlib.util as u, sys; sys.exit(not all(u.find_spec(m) for m in ("pandas", "matplotlib", "openpyxl", "pyxlsb")))' \
    || { echo "activate bioinf_tools first (python: $PY)"; exit 1; }
declare -A MD5=([make_article_figures.py]=6394c90863fd99dd0b83702753e95ea4 [cmp_ortho.py]=c6f76c9d103958d1a6c5962a3b876399 [cmp_motif.py]=812906c230c138aaa452045836698ebf
                [cmp_truthsets.py]=ca20684ffebd22c3a9ab064933d1de5c)
for f in "${!MD5[@]}"; do
    [[ $(md5sum < "$SC/$f" | cut -d' ' -f1) == "${MD5[$f]}" ]] || { echo "$f is not the expected version (${MD5[$f]})"; exit 1; }
done
for f in m6ASACseq_GSE162356_HEK293_polyA.bed.gz miCLIP_GSE63753_hek293.abcam.CIMS.m6A.9536.hg38.bed.gz \
         miCLIP_GSE63753_hek293.sysy.CITS.m6A.12051.hg38.bed.gz UBSseq_GSE225614_HeLa-shControl_sites.tsv.gz \
         UBSseq_GSE225614_HeLa-shNSUN2_sites.tsv.gz DARTseq_GSM4024125_HEK293T_APOBEC1YTH_processed.hg38.bed.gz \
         DARTseq_GSM4024129_HEK293T_APOBEC1YTHmut_processed.hg38.bed.gz; do
    [[ -s $ORTHO/$f ]] || { echo "missing truth set: $ORTHO/$f (hg38 files come from liftover_truthsets.py)"; exit 1; }
done
for arm in nofilter lab050; do for s in GM12878 $HEK; do for r in A_d212_v6 B_d130_v52 C_d212_v52; do
    [[ -s ${ARMDIR[$arm]}/$r/$s.bed ]] || { echo "missing pileup: $arm $r $s"; exit 1; }
done; done; done
mkdir -p $B/ortho/logs $B/truthsets/logs $B/motif/logs
IDFILE=$(mktemp); trap 'rm -f "$IDFILE"' EXIT
sub() { local id; id=$(sbatch --parsable "$@"); echo "$id" >> "$IDFILE"; echo "$id"; }   # runs in \$(...): IDs go via a file
log() { printf '%-50s %s\n' "$1" "$2"; }
runs() { local P=$1 s=$2; echo "--run A_d212_v6=$P/A_d212_v6/$s.bed --run C_d212_v52=$P/C_d212_v52/$s.bed --run B_d130_v52=$P/B_d130_v52/$s.bed"; }

for arm in nofilter lab050; do
    P=${ARMDIR[$arm]}
    for s in $HEK; do
        log "ortho-$arm-$s" "$(sub -p short -c 4 --mem=48G -t 2:00:00 -J ortho-$arm-$s -o "$B/ortho/logs/%x.%j.log" --wrap \
            "$PY $SC/cmp_ortho.py --sample $s --ortho-dir $ORTHO --fasta $FA --mask-bed $MASK --out $B/ortho/${s}_$arm $(runs $P $s)")"
        for tol in 0 1; do
            log "truthsets-$arm-tol$tol-$s" "$(sub -p short -c 4 --mem=48G -t 2:00:00 -J truthsets-$arm-tol$tol-$s -o "$B/truthsets/logs/%x.%j.log" --wrap \
                "$PY $SC/cmp_truthsets.py --sample $s --ortho-dir $ORTHO --fasta $FA --mask-bed $MASK --label $arm --tolerance $tol \
                 --out $B/truthsets/${s}_${arm}_tol$tol $(runs $P $s)")"
        done
    done
    for s in GM12878 $HEK; do
        o=""; [[ $s == GM12878 ]] || o="--ortho-dir $ORTHO"
        log "motif-$arm-$s" "$(sub -p short -c 2 --mem=32G -t 1:00:00 -J motif-$arm-$s -o "$B/motif/logs/%x.%j.log" --wrap \
            "$PY $SC/cmp_motif.py --sample $s --v6 $P/A_d212_v6/$s.bed --v52 $P/C_d212_v52/$s.bed --fasta $FA --mask-bed $MASK $o --label $arm --out $B/motif/${s}_$arm")"
    done
done
# figures for the article, once every job above has finished (durable storage: /scratch is purged)
ART=/projects/Genometechlab/andrew/basecall_cmp_results/article/figures; mkdir -p "$ART"
dep=$(paste -sd: "$IDFILE"); njobs=$(wc -l < "$IDFILE")
fid=$(sbatch --parsable -p short -c 2 --mem=24G -t 0:45:00 -J article-figures -o "$B/motif/logs/%x.%j.log" --dependency=afterany:$dep \
      --wrap "$PY $SC/make_article_figures.py --base $B --out $ART")
log "article-figures (after all $njobs jobs)" "$fid"
echo; echo "submitted; $(squeue --me -h | wc -l) jobs in the queue; figures will be in $ART"
