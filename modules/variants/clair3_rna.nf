// RNA small-variant calling with Clair3-RNA (step 6), per sample, on the Stage 1 BAM.
//
// SALRR's settings: the phasing model (--enable_phasing_model) and REDIportal tagging
// (--tag_variant_using_readiportal: known editing sites get FILTER=RNAEditing, not PASS).
// The model is params.clair3_rna_platform: ont_dorado_drna004 for SQK-RNA004 direct RNA
// (SALRR's cDNA data used ont_r10_dorado_cdna). Clair3-RNA calls chr1-22, X, Y by default;
// chrM is left out on purpose (a diploid model, while mitochondrial variants are
// heteroplasmic): add it through params.clair3_rna_args if wanted.
// Singularity needs absolute input/output paths and --conda_prefix (Clair3-RNA README).

process CLAIR3_RNA {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_align'
    container "${params.container_dir}/${params.clair3_rna_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.clair3_rna.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.clair3_rna.{vcf.gz,log}"), emit: out

    script:
    prefix = meta.id
    def phasing = params.clair3_rna_phasing ? '--enable_phasing_model' : ''
    """
    set -o pipefail
    /opt/bin/run_clair3_rna \\
        --bam_fn \$PWD/${bam} \\
        --ref_fn \$PWD/${fasta} \\
        --threads ${task.cpus} \\
        --platform ${params.clair3_rna_platform} \\
        ${phasing} \\
        --tag_variant_using_readiportal \\
        --sample_name ${prefix} \\
        --output_dir \$PWD/clair3 \\
        --output_prefix ${prefix} \\
        --remove_intermediate_dir \\
        --conda_prefix /opt/conda/envs/clair3_rna \\
        ${params.clair3_rna_args} \\
        2>&1 | tee ${prefix}.clair3_rna.log

    [[ -s clair3/${prefix}.vcf.gz ]] || { echo "Clair3-RNA wrote no clair3/${prefix}.vcf.gz; output dir has:" >&2; ls -la clair3 >&2; exit 1; }
    mv clair3/${prefix}.vcf.gz ${prefix}.clair3_rna.vcf.gz
    """

    // Stub test for pipeline integration
    stub:
    prefix = meta.id
    """
    echo stub | gzip > ${prefix}.clair3_rna.vcf.gz
    touch ${prefix}.clair3_rna.log
    """
}
