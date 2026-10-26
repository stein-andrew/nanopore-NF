// Per-sample read composition (bin/read_qc.sh): mapped, globin and mitochondrial reads.
// Whole-blood RNA can be mostly globin, and the share differs between samples, so it
// travels with every between-sample comparison.

process READ_QC {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_small'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(globin_bed)

    output:
    tuple val(meta), path("${meta.id}.read_qc.tsv")

    script:
    """
    read_qc.sh ${bam} ${globin_bed} ${meta.id} ${params.polya_mito_contig} > ${meta.id}.read_qc.tsv
    """

    // Stub test for pipeline integration
    stub:
    """
    printf 'sample\\tmapped_reads\\n${meta.id}\\t0\\n' > ${meta.id}.read_qc.tsv
    """
}
