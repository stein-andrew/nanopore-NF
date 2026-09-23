// GPU alignment with Parabricks. Runs from the NGC image (params.parabricks_image),
// not the site module: pbrun from the module is itself an apptainer wrapper and
// cannot run inside a container.
//
// --in-bam drops MM/ML, so the output is untagged; MODKIT_REPAIR or COPY_TAGS grafts
// the tags back. The filter runs here because the NGC image ships samtools.

process PBRUN_MINIMAP2 {
    tag "${meta.run_id}:${ref_name}"
    label 'gpu_pbrun'
    container "${params.parabricks_image}"

    input:
    tuple val(meta), path(ubam), val(ref_name), path(fasta), path(fai), path(mmi)

    output:
    tuple val(meta), val(ref_name), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    path "${prefix}.pbrun.log", emit: log

    script:
    prefix = "${meta.run_id}.${ref_name}.untagged"
    """
    mkdir -p pbtmp
    pbrun minimap2 \\
        --ref ${fasta} \\
        --index ${mmi} \\
        --in-bam ${ubam} \\
        --preset ${params.mm2_preset} \\
        ${params.pbrun_args} \\
        --out-bam raw.bam \\
        --read-group-sm ${meta.id} \\
        --read-group-pl ONT \\
        --read-group-id-prefix ${meta.run_id} \\
        --num-gpus ${params.pbrun_gpus} \\
        --num-threads ${task.cpus} \\
        --gpusort \\
        --gpuwrite \\
        --gpuwrite-deflate-algo 1 \\
        --tmp-dir \$PWD/pbtmp \\
        --logfile ${prefix}.pbrun.log

    samtools quickcheck -v raw.bam

    # pbrun copies dorado's @RG (useful provenance: run id, model, mod models) but
    # ignores --read-group-sm and leaves records untagged. Rebuild it as our ID/SM/PL
    # plus dorado's other fields, and stamp it on every record, matching the CPU path.
    rg_extra=\$(samtools view -H raw.bam | grep -m1 '^@RG' | tr '\\t' '\\n' \\
        | grep -vE '^(@RG|ID:|SM:|PL:)' | paste -sd '\\t' || true)
    samtools addreplacerg -u -w -m overwrite_all \\
        -r "@RG\\tID:${meta.run_id}\\tSM:${meta.id}\\tPL:ONT\${rg_extra:+\\t\$rg_extra}" raw.bam \\
        | samtools view -@ ${task.cpus} -b -F ${params.filter_flags} --write-index \\
            -o ${prefix}.bam##idx##${prefix}.bam.bai -
    rm -rf raw.bam pbtmp
    """

    stub:
    prefix = "${meta.run_id}.${ref_name}.untagged"
    """
    touch ${prefix}.bam ${prefix}.bam.bai ${prefix}.pbrun.log
    """
}
