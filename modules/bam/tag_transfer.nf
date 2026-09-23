// Put the uBAM's tags back onto a Parabricks alignment.
//
//   DNA  MODKIT_REPAIR    MM/ML/MN, correct on hard-clipped supplementaries;
//        COPY_DONOR_TAGS  then qs and mv verbatim (mv on primaries, which is all
//                         Clair3's dwell-time models read).
//   RNA  COPY_TAGS        any tags (pt for IsoQuant, mv, qs), verbatim; valid
//                         because RNA keeps primary alignments only.

// The donor is name-sorted once per run and reused for every reference it was
// aligned to. It keeps only `keep_tags`: MM/ML/MN for modkit repair plus whatever
// COPY_DONOR_TAGS carries afterwards (mv is most of an emit-moves uBAM).
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

process MODKIT_REPAIR {
    tag "${meta.run_id}:${ref_name}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), val(ref_name), path(bam), path(bai), path(donor), val(prefix)
    val bam_suffix      // '' when this is the final BAM; '.repaired' when COPY_DONOR_TAGS follows

    output:
    tuple val(meta), val(ref_name), path("${prefix}${bam_suffix}.bam"), path("${prefix}${bam_suffix}.bam.bai"), emit: bam
    path "${prefix}.repair.log", emit: log

    script:
    def out = "${prefix}${bam_suffix}"
    """
    samtools sort -n -@ ${task.cpus} -m 2G -T acceptor.tmp -o acceptor.ns.bam ${bam}
    modkit repair \\
        --threads ${task.cpus} \\
        --donor-bam ${donor} \\
        --acceptor-bam acceptor.ns.bam \\
        --output-bam repaired.ns.bam \\
        --log-filepath ${prefix}.repair.log
    rm acceptor.ns.bam
    # Records modkit cannot place unambiguously are dropped. In testing these were
    # short supplementary segments in tandem repeats, never primary alignments.
    grep -o 'repaired [0-9]* records, [0-9]* failed' ${prefix}.repair.log | tail -n1 || true
    samtools sort -@ ${task.cpus} -m 2G -T cs.tmp --write-index \\
        -o ${out}.bam##idx##${out}.bam.bai repaired.ns.bam
    rm repaired.ns.bam

    # Fail loudly rather than hand Stage 2 a BAM with no modification calls.
    n=\$(samtools view ${out}.bam | head -n 2000 | grep -c \$'\\tMM:Z:' || true)
    echo "MM:Z: present in \$n of the first 2000 records"
    [[ \$n -gt 0 ]] || { echo "FATAL: repaired BAM has no MM tags" >&2; exit 1; }
    """

    stub:
    def out = "${prefix}${bam_suffix}"
    """
    touch ${out}.bam ${out}.bam.bai ${prefix}.repair.log
    """
}

// DNA, after MODKIT_REPAIR: the tags modkit repair does not carry (qs, mv), from
// the already name-sorted donor. copy_tags.py skips positional tags (mv) on the
// hard-clipped supplementaries, so mv lands on primaries only.
process COPY_DONOR_TAGS {
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

process COPY_TAGS {
    tag "${meta.run_id}:${ref_name}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), val(ref_name), path(bam), path(bai), path(ubam), val(prefix)
    val tags

    output:
    tuple val(meta), val(ref_name), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam

    script:
    """
    copy_tags.py \\
        --ubam ${ubam} \\
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
