process DORADO_BASECALLER {
    tag "${meta.run_id}"
    label 'gpu_dorado'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(meta), path(pod5)         // a .pod5 file or a directory of them
    val mod_models

    output:
    tuple val(meta), path("${meta.run_id}.ubam.bam"), env('DORADO_VERSION'), emit: ubam

    script:
    // Models are baked into the image under /opt/models; full paths, not names, so
    // dorado never tries to download at run time.
    def models = '/opt/models'
    def mods   = mod_models
        ? '--modified-bases-models ' + mod_models.collect { m -> "${models}/${params.dorado_model}_${m}" }.join(',')
        : ''
    """
    DORADO_VERSION=\$(dorado --version 2>&1 | tail -n1)
    nvidia-smi -L

    dorado basecaller \\
        ${models}/${params.dorado_model} \\
        ${pod5} \\
        --recursive \\
        --device ${params.dorado_device} \\
        ${mods} \\
        ${params.emit_moves ? '--emit-moves' : ''} \\
        ${params.estimate_poly_a ? '--estimate-poly-a' : ''} \\
        ${params.dorado_args} \\
        > ${meta.run_id}.ubam.bam

    # -u: a uBAM has no @SQ lines, which quickcheck otherwise treats as a failure
    samtools quickcheck -u -v ${meta.run_id}.ubam.bam
    """

    // Stub test for pipeline integration
    stub:
    """
    DORADO_VERSION=stub
    touch ${meta.run_id}.ubam.bam
    """
}
