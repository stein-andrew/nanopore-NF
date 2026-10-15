// StringTie, per sample, as SALRR runs it: guided long-read discovery (-L -G, no -e) with
// -c 2.5 -f 0.05 -t. params.stringtie_args holds those flags plus --fr for direct RNA
// (reads are sense-strand; SALRR's cDNA was unstranded, so it had no strand flag).

process STRINGTIE {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.stringtie.gtf")

    script:
    """
    stringtie ${bam} ${params.stringtie_args} -G ${gtf} -p ${task.cpus} -o ${meta.id}.stringtie.gtf
    """

    // Stub test for pipeline integration
    stub:
    """
    touch ${meta.id}.stringtie.gtf
    """
}
