// Known vs novel for Clair3-RNA calls (step 6b), per sample, GRCh38 only.
//
// Calls are split to one ALT per record and left-aligned against the reference, as gnomAD and
// ClinVar are; only the call positions are then looked up in each reference, through its
// index. bin/known_novel.py matches on exact CHROM, POS, REF and ALT and classes each call:
// editing (REDIportal tag, split off first as SALRR does) / known (in gnomAD and passing QC:
// nuclear joint FILTER not BOTH_FILTERED; chrM PASS) / novel / lowqual / refcall, attaching
// ClinVar significance wherever it matches. The classes and frequencies are written back into
// the VCF as INFO fields (bcftools annotate with REF/ALT columns: exact allele matching).
// References: scripts/fetch_variant_refs.sh (cluster profile: gnomad_nuclear, gnomad_mito, clinvar).

process KNOWN_NOVEL {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_small'
    container "${params.container_dir}/${params.variants_image}"

    input:
    tuple val(meta), path(vcf), path(fasta), path(fai)
    tuple path(gnomad), path(gnomad_csi), path(mito), path(mito_csi), path(clinvar), path(clinvar_csi)

    output:
    tuple val(meta), path("${prefix}.{variants.vcf.gz,variants.vcf.gz.tbi,variants.tsv.gz,variants.summary.tsv}")

    script:
    prefix = meta.id
    """
    set -o pipefail

    bcftools norm -f ${fasta} -m -any --threads ${task.cpus} -Ob -o calls.norm.bcf ${vcf}
    bcftools index -f calls.norm.bcf
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' calls.norm.bcf | sort -k1,1 -k2,2n -u > calls.bed

    # each reference restricted to the call positions (index lookups, not a scan)
    if [[ -s calls.bed ]]; then
        bcftools view -R calls.bed ${gnomad} \\
          | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/AF_joint\\t%INFO/AF_exomes\\t%INFO/AF_genomes\\t%INFO/nhomalt_joint\\n' > gnomad.tsv
        bcftools view -R calls.bed ${mito} \\
          | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/AF_hom\\t%INFO/AF_het\\n' > mito.tsv
        bcftools view -R calls.bed ${clinvar} -Ou | bcftools norm -f ${fasta} -m -any -Ou \\
          | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\t%INFO/CLNDN\\t%INFO/GENEINFO\\n' > clinvar.tsv
    else
        : > gnomad.tsv; : > mito.tsv; : > clinvar.tsv
    fi
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%QUAL\\t[%GT\\t%GQ\\t%DP\\t%AF]\\n' calls.norm.bcf > calls.tsv

    known_novel.py --calls calls.tsv --gnomad gnomad.tsv --mito mito.tsv --clinvar clinvar.tsv \\
        --mito-contig ${params.variants_mito_contig} --dp-split ${params.variants_dp_split} --prefix ${prefix}

    bgzip -c ${prefix}.annot.tsv > annot.tsv.gz
    tabix -f -s1 -b2 -e2 annot.tsv.gz
    bcftools annotate -a annot.tsv.gz -h ${prefix}.annot.hdr \\
        -c CHROM,POS,REF,ALT,INFO/VARIANT_CLASS,INFO/gnomAD_AF,INFO/gnomAD_AF_exomes,INFO/gnomAD_AF_genomes,INFO/gnomAD_FILTER,INFO/CLNSIG,INFO/CLNREVSTAT \\
        -Oz -o ${prefix}.variants.vcf.gz calls.norm.bcf
    tabix -f -p vcf ${prefix}.variants.vcf.gz
    """

    // Stub test for pipeline integration
    stub:
    prefix = meta.id
    """
    echo stub | gzip > ${prefix}.variants.vcf.gz
    touch ${prefix}.variants.vcf.gz.tbi ${prefix}.variants.summary.tsv
    echo stub | gzip > ${prefix}.variants.tsv.gz
    """
}
