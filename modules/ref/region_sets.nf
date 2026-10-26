// Region sets from the annotation, once per reference (bin/region_sets.sh): the globin genes
// (for the per-sample globin fraction) and the variant mask (IG_* / TR_* genes from the
// annotation's own gene_type, plus the MHC region as HLA).

process REGION_SETS {
    tag "${ref_name}"
    label 'cpu_small'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(ref_name), path(gtf)

    output:
    tuple val(ref_name), path("${ref_name}.globin.bed"), path("${ref_name}.mask.bed"), emit: regions
    tuple val(ref_name), path("${ref_name}.{globin.bed,mask.bed,region_sets.log}"), emit: out

    script:
    def genes = params.globin_genes instanceof List ? params.globin_genes.join(',') : params.globin_genes.toString()
    """
    region_sets.sh ${gtf} ${genes} '${params.variant_mask_mhc ?: ''}' ${ref_name}
    """

    // Stub test for pipeline integration
    stub:
    """
    touch ${ref_name}.globin.bed ${ref_name}.mask.bed ${ref_name}.region_sets.log
    """
}
