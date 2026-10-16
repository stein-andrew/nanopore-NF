// Cohort isoform curation, SALRR:
//   COHORT_CATALOGUE  merge every sample's consensus with isomatch into one catalogue,
//                     classify it against the annotation, and add gene lines so IsoQuant
//                     can requantify against it with --complete_genedb.
//   COHORT_FILTER     per-sample CPM from the requantification; keep transcripts whose
//                     median CPM across samples > params.cohort_min_median_cpm.
//
// isomatch labels its inputs S1..Sn and records each one's file name in the ##ISOM <SAMPLE>
// header; the file names start with the sample ID (ST1918.merged.gtf.gz), so that header is
// the S-number -> sample key for the track and present_absent tables.

process COHORT_CATALOGUE {
    tag "${ref_name}"
    label 'cpu_small'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(ref_name), path(gtfs), path(fasta), path(fai), path(ref_gtf)

    output:
    tuple val(ref_name), path("${prefix}.catalogue.gtf"), emit: catalogue
    tuple val(ref_name), path("${prefix}.classify.classification.txt.gz"), emit: classification
    tuple val(ref_name),
          path("${prefix}.{merged.gtf.gz,track.tsv.gz,present_absent.tsv.gz,merged_info.json,merged_params.json,catalogue.gtf}"),
          path("${prefix}.classify.*"), emit: out

    script:
    prefix = "cohort.${ref_name}"
    """
    isomatch merge --ref-fa ${fasta} -o ${prefix} ${params.isomatch_merge_args} ${gtfs}
    isomatch classify --ref-fa ${fasta} --ref-gtf ${ref_gtf} -o ${prefix}.classify ${params.isomatch_classify_args} \\
        ${prefix}.merged.gtf.gz
    gtf_add_genes.py ${prefix}.merged.gtf.gz > ${prefix}.catalogue.gtf
    """

    // Stub test for pipeline integration
    stub:
    prefix = "cohort.${ref_name}"
    """
    touch ${prefix}.merged.gtf.gz ${prefix}.track.tsv.gz ${prefix}.present_absent.tsv.gz ${prefix}.merged_info.json ${prefix}.merged_params.json
    touch ${prefix}.catalogue.gtf ${prefix}.classify.classification.txt.gz ${prefix}.classify.annotated.gtf.gz
    """
}

process COHORT_FILTER {
    tag "${ref_name}"
    label 'cpu_small'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(ref_name), val(ids), path(counts), path(catalogue), path(classification)

    output:
    tuple val(ref_name),
          path("${prefix}.{counts.tsv,cpm.tsv,high_confidence.gtf,high_confidence.classification.tsv,summary.tsv}")

    script:
    prefix = "cohort.${ref_name}"
    // One sample arrives as a single path rather than a list.
    def count_args = [[ids].flatten(), [counts].flatten()].transpose()
        .collect { id, f -> "--counts ${id}=${f}" }
        .join(' ')
    """
    cohort_cpm.py --gtf ${catalogue} --classification ${classification} \\
        ${count_args} --min ${params.cohort_min_median_cpm} --prefix ${prefix}
    """

    // Stub test for pipeline integration
    stub:
    prefix = "cohort.${ref_name}"
    """
    touch ${prefix}.counts.tsv ${prefix}.cpm.tsv ${prefix}.high_confidence.gtf ${prefix}.high_confidence.classification.tsv ${prefix}.summary.tsv
    """
}
