// .fai for Stage 2, from the same cache REFERENCE_INDEX writes in Stage 1
// (params.index_dir/<ref>/<fasta>.fai), so it is built at most once per cluster. No .mmi:
// Stage 2 never aligns.

process FAIDX {
    tag "${ref_name}"
    label 'cpu_small'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(ref_name), path(fasta)

    output:
    tuple val(ref_name), path(fasta), path("${fasta}.fai")

    script:
    def dir = "${params.index_dir}/${ref_name}"
    def fai = "${dir}/${fasta.name}.fai"
    """
    mkdir -p ${dir}
    if [[ ! -s ${fai} ]]; then
        samtools faidx ${fasta} --fai-idx ${fai}.tmp.\$\$
        mv ${fai}.tmp.\$\$ ${fai}
    fi
    ln -s ${fai} ${fasta}.fai
    """

    // Stub test for pipeline integration
    stub:
    """
    touch ${fasta}.fai
    """
}
