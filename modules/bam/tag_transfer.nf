// Put the uBAM's tags back onto a Parabricks alignment (pbrun --in-bam drops them).
//
// copy_tags.py copies every carried tag verbatim (MM/ML/MN, qs, mv, pt) and, with
// --unclip, puts the bases Parabricks hard-clipped off supplementaries back from the
// uBAM (H -> S), exactly as minimap2 -Y writes them. Positional tags then stay valid
// on every record, and the output matches the CPU path record for record. This
// replaced `modkit repair`, which aborts intermittently while finishing its output
// (0.5.0 and 0.6.4, "malloc(): unaligned fastbin chunk detected"; AICR 2026-09-23).

// The donor is name-sorted once per run and reused for every reference it was
// aligned to. Only the carried tags are kept (SEQ/QUAL stay, for --unclip).
process DONOR_NAMESORT {
    tag "${meta.run_id}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), path(ubam)
    val keep_tags

    output:
    tuple val(meta), path("${meta.run_id}.donor.ns.bam")

    script:
    """
    samtools view -u -@ ${task.cpus} --keep-tag ${keep_tags} ${ubam} \\
        | samtools sort -n -@ ${task.cpus} -m 2G -T donor.tmp -o ${meta.run_id}.donor.ns.bam -
    """

    stub:
    """
    touch ${meta.run_id}.donor.ns.bam
    """
}

process COPY_TAGS {
    tag "${meta.run_id}:${ref_name}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), val(ref_name), path(bam), path(bai), path(donor), val(prefix)
    val tags

    output:
    tuple val(meta), val(ref_name), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam

    script:
    """
    copy_tags.py \\
        --ubam ${donor} \\
        --ubam-name-sorted \\
        --unclip \\
        --aligned ${bam} \\
        --out ${prefix}.bam \\
        --tags ${tags} \\
        --threads ${task.cpus} \\
        --tmpdir \$PWD
    """

    stub:
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    """
}
