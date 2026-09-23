// Put the uBAM's tags back onto a Parabricks alignment.
//
//   DNA  MODKIT_REPAIR  MM/ML/MN only, but correct on hard-clipped supplementaries.
//   RNA  COPY_TAGS      any tags (mv for Clair3-RNA, pt for IsoQuant), verbatim;
//                       valid because RNA keeps primary alignments only.

// The donor is name-sorted once per run and reused for every reference it was
// aligned to. Only MM/ML/MN are kept: on an emit-moves uBAM the move table is most
// of the file, and modkit repair reads nothing else.
process DONOR_NAMESORT {
    tag "${meta.run_id}"
    label 'cpu_sort'

    input:
    tuple val(meta), path(ubam)

    output:
    tuple val(meta), path("${meta.run_id}.donor.ns.bam")

    script:
    """
    samtools view -u -@ ${task.cpus} --keep-tag MM,ML,MN ${ubam} \\
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

    input:
    tuple val(meta), val(ref_name), path(bam), path(bai), path(donor), val(prefix)

    output:
    tuple val(meta), val(ref_name), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    path "${prefix}.repair.log", emit: log

    script:
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
        -o ${prefix}.bam##idx##${prefix}.bam.bai repaired.ns.bam
    rm repaired.ns.bam

    # Fail loudly rather than hand Stage 2 a BAM with no modification calls.
    n=\$(samtools view ${prefix}.bam | head -n 2000 | grep -c \$'\\tMM:Z:' || true)
    echo "MM:Z: present in \$n of the first 2000 records"
    [[ \$n -gt 0 ]] || { echo "FATAL: repaired BAM has no MM tags" >&2; exit 1; }
    """

    stub:
    """
    touch ${prefix}.bam ${prefix}.bam.bai ${prefix}.repair.log
    """
}

process COPY_TAGS {
    tag "${meta.run_id}:${ref_name}"
    label 'cpu_sort'

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
