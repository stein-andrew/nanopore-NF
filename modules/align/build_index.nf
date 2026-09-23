// Reference indexes, built once per (reference, preset) and cached in
// params.index_dir via storeDir, so the shared reference directory is never written to.

process SAMTOOLS_FAIDX {
    tag "${ref_name}"
    label 'cpu_small'
    storeDir "${params.index_dir}/${ref_name}"

    input:
    tuple val(ref_name), path(fasta)

    output:
    tuple val(ref_name), path("${fasta}.fai")

    script:
    """
    samtools faidx ${fasta} --fai-idx ${fasta}.fai
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

process MINIMAP2_INDEX {
    tag "${ref_name}"
    label 'cpu_align'
    storeDir "${params.index_dir}/${ref_name}"

    input:
    tuple val(ref_name), path(fasta)
    val index_args

    output:
    tuple val(ref_name), path("${ref_name}.${index_tag}.mmi")

    script:
    // e.g. '-x splice -k14 -I 16G' -> 'splice_k14_I16G'
    index_tag = index_args.replaceAll(/-x\s+/, '').replaceAll(/[-:\s]+/, '_').replaceAll(/^_|_$/, '')
    """
    minimap2 ${index_args} -t ${task.cpus} -d ${ref_name}.${index_tag}.mmi ${fasta}
    """

    stub:
    index_tag = index_args.replaceAll(/-x\s+/, '').replaceAll(/[-:\s]+/, '_').replaceAll(/^_|_$/, '')
    """
    touch ${ref_name}.${index_tag}.mmi
    """
}
