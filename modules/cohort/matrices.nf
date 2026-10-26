// Cohort matrices for comparing samples, per reference (bin/cohort_matrices.py,
// bin/mods_site_values.sh): one row per feature, one column per sample.
//   cohort.<ref>.gene_counts.tsv.gz, gene_cpm.tsv.gz     IsoQuant requantification, catalogue genes
//   cohort.<ref>.mods.{n_valid,n_mod,pct_fpcorr}.tsv.gz  every final-stage site in any sample;
//                                                        empty = not measured in that sample
//   cohort.<ref>.mods.summary.tsv                        final and measured sites per sample
//   cohort.<ref>.polya.{isoform,gene}.{median_pt,n_tails}.tsv.gz
//   cohort.<ref>.read_qc.tsv                             mapped / globin / mito reads per sample
// Samples are in sample-id order, whatever order their tasks finished in.

process COHORT_MATRICES {
    tag "${ref_name}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(ref_name), val(ids),
          path(gene_counts, stageAs: 'genes/*', arity: '1..*'), path(fpcorr, stageAs: 'mods/*', arity: '1..*'),
          path(polya, stageAs: 'polya/*', arity: '1..*'), path(read_qc, stageAs: 'qc/*', arity: '1..*'),
          path(catalogue), path(classification)

    output:
    tuple val(ref_name), path("${prefix}.*"), emit: out

    script:
    prefix = "cohort.${ref_name}"
    def counts_args = [ids, gene_counts].transpose().collect { id, f -> "--counts ${id}=${f}" }.join(' ')
    def mods_pairs  = [ids, fpcorr].transpose().collect { id, f -> "${id}=${f}" }.join(' ')
    def polya_args  = [ids, polya].transpose().collect { id, f -> "--features ${id}=${f}" }.join(' ')
    """
    cohort_matrices.py genes ${counts_args} --catalogue ${catalogue} --classification ${classification} --prefix ${prefix}

    mods_site_values.sh ${params.mods_min_percent} mods_sites.tsv mods_site_values.tsv ${mods_pairs}
    cohort_matrices.py mods --sites mods_sites.tsv --long mods_site_values.tsv --samples ${ids.join(' ')} \\
        --min-pct ${params.mods_min_percent} --prefix ${prefix}

    cohort_matrices.py polya ${polya_args} --prefix ${prefix}

    { head -n 1 ${read_qc[0]}; for f in ${read_qc.join(' ')}; do tail -n +2 "\$f"; done; } > ${prefix}.read_qc.tsv
    """

    // Stub test for pipeline integration
    stub:
    prefix = "cohort.${ref_name}"
    """
    touch ${prefix}.gene_counts.tsv.gz ${prefix}.read_qc.tsv
    """
}
