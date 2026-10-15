// SALRR: only primary alignments with MAPQ >= 40 go into transcript reconstruction.
// Stage 1 already dropped unmapped/secondary/supplementary (-F 2308); repeating it here
// keeps the filter true for a BAM from anywhere else. This BAM is for isoform curation
// only; modification and variant steps read the Stage 1 BAM.

process FILTER_MAPQ {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_small'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.flagstat"), emit: flagstat

    script:
    prefix = "${meta.id}.${meta.ref}.mapq${params.isoform_min_mapq}"
    """
    samtools view -@ ${task.cpus} -b -F 2308 -q ${params.isoform_min_mapq} -o ${prefix}.bam ${bam}
    samtools index -@ ${task.cpus} ${prefix}.bam
    samtools flagstat -@ ${task.cpus} ${prefix}.bam > ${prefix}.flagstat
    """

    // Stub test for pipeline integration
    stub:
    prefix = "${meta.id}.${meta.ref}.mapq${params.isoform_min_mapq}"
    """
    touch ${prefix}.bam ${prefix}.bam.bai ${prefix}.flagstat
    """
}
