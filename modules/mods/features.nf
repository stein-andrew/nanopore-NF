// Modification levels per isoform and per gene (step 5a).
//
// Candidate sites: the sample-level FP-corrected stage at >= params.feature_mods_candidate_pct
// (default 20: the final-stage sites). modkit extract calls gets every read's call at just
// those sites, at the sample's own thresholds (the same flags as its pileup), and
// bin/mod_feature_counts.sh counts them per site, per isoform and per gene of each read as
// they stream past, so the per-read table (sites x depth rows) never exists.
// bin/feature_mods.py then checks the all-reads counts reproduce the pileup (or stops),
// computes per-feature percentages, FP-corrects them with fp_correct.py's own code, and
// applies the final-stage rule (>= mods_min_coverage reads, >= mods_min_percent) per feature.
//
// Calls come from the Stage 1 BAM; which isoform and gene a read belongs to comes from the
// IsoQuant requantification (MAPQ >= 40 reads): isoform = unique assignment, as in the
// unique_only counts; gene = a single gene across the read's assignments.

process MOD_FEATURES {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(fpcorr), path(modkit_flags), path(givt_flags), path(assignments),
          path(catalogue), path(classification), path(fasta), path(fai), path(globin_bed)
    path(error_table)

    output:
    tuple val(meta),
          path("${prefix}.{isoform_mods.bed,gene_mods.bed,feature_mods.tsv.gz,feature_mods.summary.tsv,extract_calls.log}")

    script:
    prefix = meta.id
    def cand = params.feature_mods_candidate_pct
    """
    set -o pipefail      # a modkit failure mid-stream must fail the task, not leave partial counts

    read_features.py ${assignments} > ${prefix}.read_features.tsv
    # Globin loci are left out: in whole blood their depth exceeds modkit_max_depth, so the
    # pileup there is downsampled and can never agree with uncapped per-read counts.
    bed_exclude.sh ${globin_bed} ${fpcorr} > ${prefix}.candidates.bed
    feature_mods.py positions --candidates ${prefix}.candidates.bed --candidate-pct ${cand} > ${prefix}.positions.bed

    if [[ -s ${prefix}.positions.bed ]]; then
        modkit extract calls ${bam} - \\
            --include-bed ${prefix}.positions.bed \\
            --mapped-only \\
            \$(cat ${modkit_flags}) \\
            --threads ${task.cpus} \\
            --log-filepath ${prefix}.extract_calls.log \\
            --suppress-progress \\
          | mod_feature_counts.sh ${prefix}.read_features.tsv > ${prefix}.counts.tsv
    else
        : > ${prefix}.counts.tsv
        echo "no candidate sites: nothing to extract" > ${prefix}.extract_calls.log
    fi

    feature_mods.py levels \\
        --counts ${prefix}.counts.tsv \\
        --candidates ${prefix}.candidates.bed --candidate-pct ${cand} \\
        --catalogue ${catalogue} --classification ${classification} \\
        --reference ${fasta} --errortable ${error_table} --givt-flags ${givt_flags} \\
        --min-cov ${params.mods_min_coverage} --min-pct ${params.mods_min_percent} \\
        --min-agreement ${params.feature_mods_min_agreement} \\
        --prefix ${prefix}
    """

    // Stub test for pipeline integration
    stub:
    prefix = meta.id
    """
    touch ${prefix}.isoform_mods.bed ${prefix}.gene_mods.bed ${prefix}.feature_mods.summary.tsv ${prefix}.extract_calls.log
    echo stub | gzip > ${prefix}.feature_mods.tsv.gz
    """
}
