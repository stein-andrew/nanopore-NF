// Random subsampling before modkit pileup (bin/subsample_loci.sh), so that its --max-depth cap
// never decides which reads count: with a cap, modkit keeps the first reads at each position
// in alignment order, which biased whole-blood globin site values (ST1918, 2026-10-27). Loci
// deeper than params.modkit_max_depth (reads overlapping 1 kb windows, spliced spans
// included) are randomly subsampled to just under it; every other read is kept. With no such
// locus, or with modkit uncapped, the original BAM is passed on as it is.
//
// Only modkit pileup reads the result. sample-probs (the call thresholds) still samples the
// original BAM, and MOD_FEATURES still counts every read; its pileup check exempts these loci.

process SUBSAMPLE_LOCI {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_sort'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.pileup.bam"), path("${prefix}.pileup.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.loci.bed"), emit: loci
    tuple val(meta), path("${prefix}.{loci.bed,subsample.log}"), emit: out

    script:
    prefix = meta.id
    """
    subsample_loci.sh run ${bam} ${bai} ${params.modkit_max_depth} ${params.subsample_window} \\
        ${params.subsample_merge} ${params.subsample_seed} ${task.cpus} ${prefix}
    """

    // Stub test for pipeline integration
    stub:
    prefix = meta.id
    """
    touch ${prefix}.pileup.bam ${prefix}.pileup.bam.bai ${prefix}.loci.bed ${prefix}.subsample.log
    """
}
