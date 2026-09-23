// Merge flowcells of one sample, then checksum and flagstat the final BAMs.

process MERGE_RUNS {
    tag "${meta.id}:${ref_name}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), val(ref_name), path(bams, stageAs: 'runs/*'), path(bais, stageAs: 'runs/*')

    output:
    tuple val(meta), val(ref_name), path("${prefix}.bam"), path("${prefix}.bam.bai")

    script:
    // Read groups stay distinct: each run was aligned with its own RG ID.
    prefix = "${meta.id}.${ref_name}"
    """
    samtools merge -@ ${task.cpus} --write-index -o ${prefix}.bam##idx##${prefix}.bam.bai ${bams}
    """

    stub:
    prefix = "${meta.id}.${ref_name}"
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    """
}

process BAM_QC {
    tag "${meta.id}:${ref_name}"
    label 'cpu_small'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), val(ref_name), path(bam), path(bai)

    output:
    tuple val(meta), val(ref_name), path("${bam}.flagstat"), path("${bam}.md5")

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${bam}.flagstat
    md5sum ${bam} > ${bam}.md5
    """

    stub:
    """
    touch ${bam}.flagstat
    echo "00000000000000000000000000000000  ${bam}" > ${bam}.md5
    """
}
