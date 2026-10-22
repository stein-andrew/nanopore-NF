// RNA modifications: modkit pileup, then the filter stages and gIVT false-positive correction.
//
// MODKIT_PILEUP  modkit sample-probs estimates per-base call thresholds (what pileup would
//                estimate itself); bin/mod_thresholds.py writes them as flags for BOTH modkit
//                pileup and gIVT, because gIVT's false-positive rates are only right at the
//                threshold the calls were made at. params.mod_filter_thresholds fixes them
//                for every sample instead (sample-probs still runs; its estimate is recorded).
//                --only-mapped: pileup's own estimate skips soft-clipped bases too.
// MOD_STAGES     raw -> Nvalid_cov >= C -> (and) percent >= P; FP correction on the coverage
//                stage -> FP-corrected percent >= P; a per-stage, per-mod-code row count.
//                Correcting the coverage stage, not the raw pileup, gives the same corrected
//                stage from far fewer rows.
//
// The correction is bin/fp_correct.py: a vectorized drop-in for gIVT's IVT_fp_correction.py
// (github.com/genometechlab/gIVT, kept in bin/ for comparison) with the same arguments and
// corrected values, checked row for row by tests/compare_fpcorr.py. It reads the reference
// by memory map through the .fai instead of loading the genome, and parses with polars on
// task.cpus threads. Both processes run in the mods image (modkit 0.5.0 + Python + polars).

process MODKIT_PILEUP {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_align'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.raw.bed"), path("${prefix}.givt_flags.txt"), emit: pileup
    tuple val(meta), path("${prefix}.modkit_flags.txt"), path("${prefix}.givt_flags.txt"), emit: flags
    tuple val(meta),
          path("${prefix}.{raw.bed,thresholds.tsv,modkit_flags.txt,givt_flags.txt,sample_probs.txt,sample_probs.log,pileup.log}"),
          path("${prefix}.sample_probs/*", arity: '0..*'), emit: out

    script:
    prefix = meta.id
    def fixed = params.mod_filter_thresholds ? "--fixed '${params.mod_filter_thresholds}'" : ''
    """
    modkit sample-probs ${bam} \\
        --threads ${task.cpus} \\
        --only-mapped \\
        --percentiles ${params.modkit_filter_percentile} \\
        --num-reads ${params.modkit_sample_reads} \\
        --hist \\
        --out-dir ${prefix}.sample_probs \\
        --prefix ${prefix} \\
        --log-filepath ${prefix}.sample_probs.log \\
        --suppress-progress \\
        > ${prefix}.sample_probs.txt

    # The threshold table has been on stdout and in --out-dir in different modkit versions;
    # mod_thresholds.py fails rather than guess if neither parses.
    table=\$(ls ${prefix}.sample_probs/*thresholds*.tsv 2>/dev/null | head -n1 || true)
    mod_thresholds.py --sample-probs "\${table:-${prefix}.sample_probs.txt}" \\
        --percentile ${params.modkit_filter_percentile} --prefix ${prefix} ${fixed}

    modkit pileup ${bam} ${prefix}.raw.bed \\
        --threads ${task.cpus} \\
        --max-depth ${params.modkit_max_depth} \\
        \$(cat ${prefix}.modkit_flags.txt) \\
        --log-filepath ${prefix}.pileup.log \\
        --suppress-progress \\
        ${params.modkit_pileup_args}
    """

    // Stub test for pipeline integration
    stub:
    prefix = meta.id
    """
    mkdir -p ${prefix}.sample_probs
    touch ${prefix}.raw.bed ${prefix}.thresholds.tsv ${prefix}.modkit_flags.txt ${prefix}.sample_probs.txt ${prefix}.sample_probs.log ${prefix}.pileup.log
    echo 'a,0.9' > ${prefix}.givt_flags.txt
    """
}

process MOD_STAGES {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_small'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(meta), path(raw), path(givt_flags), path(fasta), path(fai)
    path(error_table)

    output:
    // the FP-corrected coverage stage goes on to the isoform/gene-level step (features.nf)
    tuple val(meta), path("${prefix}.cov${min_cov}_fpcorr.bed"), emit: fpcorr
    tuple val(meta),
          path("${prefix}.{cov${min_cov}.bed,cov${min_cov}_pct${min_pct}.bed,cov${min_cov}_fpcorr.bed,cov${min_cov}_fpcorr_pct${min_pct}.bed,mod_stages.tsv,fpcorr.log.gz}"), emit: out

    script:
    prefix  = meta.id
    min_cov = params.mods_min_coverage
    min_pct = params.mods_min_percent
    def cov = "${prefix}.cov${min_cov}"
    """
    mod_stages.sh coverage ${raw} ${min_cov} ${min_pct} ${prefix}

    # An empty coverage stage just stays empty. The log holds fp_correct.py's per-step timings
    # and its count of sites per reason (corrected / no 9mer / base mismatch / not in table).
    export POLARS_MAX_THREADS=${task.cpus}
    if [[ -s ${cov}.bed ]]; then
        fp_correct.py -m ${cov}.bed -ref ${fasta} -e ${error_table} \\
            -mt \$(cat ${givt_flags}) -o ${cov}_fpcorr.bed > ${prefix}.fpcorr.log 2>&1
    else
        : > ${cov}_fpcorr.bed
        echo "no rows with Nvalid_cov >= ${min_cov}: nothing to correct" > ${prefix}.fpcorr.log
    fi
    gzip -f ${prefix}.fpcorr.log

    mod_stages.sh fpcorr ${cov}_fpcorr.bed ${min_pct} ${cov}_fpcorr_pct${min_pct}.bed
    mod_stages.sh summary ${prefix}.mod_stages.tsv \\
        raw=${raw} \\
        cov${min_cov}=${cov}.bed \\
        cov${min_cov}_pct${min_pct}=${cov}_pct${min_pct}.bed \\
        cov${min_cov}_fpcorr=${cov}_fpcorr.bed \\
        cov${min_cov}_fpcorr_pct${min_pct}=${cov}_fpcorr_pct${min_pct}.bed
    """

    // Stub test for pipeline integration
    stub:
    prefix  = meta.id
    min_cov = params.mods_min_coverage
    min_pct = params.mods_min_percent
    def cov = "${prefix}.cov${min_cov}"
    """
    touch ${cov}.bed ${cov}_pct${min_pct}.bed ${cov}_fpcorr.bed ${cov}_fpcorr_pct${min_pct}.bed ${prefix}.mod_stages.tsv
    echo stub | gzip > ${prefix}.fpcorr.log.gz
    """
}
