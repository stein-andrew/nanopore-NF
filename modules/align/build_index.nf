// Reference indexes (.fai and a minimap2 .mmi per preset), built once into
// params.index_dir and reused by every later run; the shared reference directory
// is never written to for overwrite protection.
//

process REFERENCE_INDEX {
    tag "${ref_name}"
    label 'cpu_index'
    container "${params.container_dir}/${params.stage1_image}"

    input:
    tuple val(ref_name), path(fasta)
    val index_args

    output:
    tuple val(ref_name), path("${fasta}.fai"), path("${ref_name}.${index_tag}.mmi")

    script:
    // e.g. '-x splice -k14 -I 16G' -> 'splice_k14_I_16G'
    index_tag = index_args.replaceAll(/-x\s+/, '').replaceAll(/[-:\s]+/, '_').replaceAll(/^_|_$/, '')
    def dir = "${params.index_dir}/${ref_name}"
    def fai = "${dir}/${fasta.name}.fai"
    def mmi = "${dir}/${ref_name}.${index_tag}.mmi"
    """
    mkdir -p ${dir}
    if [[ ! -s ${fai} ]]; then
        samtools faidx ${fasta} --fai-idx ${fai}.tmp.\$\$
        mv ${fai}.tmp.\$\$ ${fai}
    fi
    if [[ ! -s ${mmi} ]]; then
        minimap2 ${index_args} -t ${task.cpus} -d ${mmi}.tmp.\$\$ ${fasta}
        mv ${mmi}.tmp.\$\$ ${mmi}
    fi
    ln -s ${fai} ${fasta}.fai
    ln -s ${mmi} ${ref_name}.${index_tag}.mmi
    """

    // Stub test for pipeline integration
    stub:
    index_tag = index_args.replaceAll(/-x\s+/, '').replaceAll(/[-:\s]+/, '_').replaceAll(/^_|_$/, '')
    """
    touch ${fasta}.fai ${ref_name}.${index_tag}.mmi
    """
}
