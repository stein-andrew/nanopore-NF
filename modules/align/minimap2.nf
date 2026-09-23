// CPU alignment. Tags ride through FASTQ comments (samtools fastq -T | minimap2 -y),
// so the output already carries MM/ML and needs no tag transfer.

process MINIMAP2 {
    tag "${meta.run_id}:${ref_name}"
    label 'cpu_align'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), path(ubam), val(ref_name), path(mmi), val(prefix)
    val carry_tags

    output:
    tuple val(meta), val(ref_name), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.run_id}\\tSM:${meta.id}\\tPL:ONT"
    def sort_threads = Math.max(1, (task.cpus / 4) as int)
    """
    samtools fastq -@ 4 -T ${carry_tags} ${ubam} \\
        | minimap2 -y -ax ${params.mm2_preset} ${params.mm2_args} -R '${rg}' -t ${task.cpus} ${mmi} - \\
        | samtools view -u -F ${params.filter_flags} - \\
        | samtools sort -@ ${sort_threads} -m 2G -T sort.tmp --write-index \\
            -o ${prefix}.bam##idx##${prefix}.bam.bai -
    """

    stub:
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    """
}
