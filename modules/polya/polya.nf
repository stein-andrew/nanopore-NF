// Poly(A) tail lengths per read, per isoform and per gene (step 5b), by SALRR's rules.
//
// Tails are Dorado's pt:i on each primary alignment of the Stage 1 BAM. As in SALRR, chrM
// reads are kept apart from nuclear ones, and a tail is dropped when missing, not found
// (pt <= 0) or censored at the estimator's cap. Dorado's direct-RNA estimates have no cap
// (checked on these data: the top values are single reads, no pile at a maximum), so the
// caps default to off; set polya_cap_nuclear / polya_cap_mito (400 / 250 in SALRR) for cDNA.
// Isoform and gene come from the IsoQuant requantification exactly as for the modification
// levels: isoform = unique assignment, gene = a single gene (bin/read_features.py).

process POLYA_FEATURES {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_small'
    container "${params.container_dir}/${params.mods_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(assignments), path(catalogue), path(classification)

    output:
    tuple val(meta),
          path("${prefix}.{isoform_polya.bed,gene_polya.bed,polya_features.tsv.gz,polya_reads.tsv.gz,polya_hist.tsv,polya_summary.tsv}")

    script:
    prefix = meta.id
    """
    set -o pipefail

    read_features.py ${assignments} > ${prefix}.read_features.tsv

    # read_id, contig, pt ('NA' if the read has no pt:i tag), primary alignments only
    samtools view -F 2308 ${bam} \\
      | awk -F'\\t' 'BEGIN { OFS = "\\t" }
                     { pt = "NA"
                       for (i = 12; i <= NF; i++) if (substr(\$i, 1, 5) == "pt:i:") { pt = substr(\$i, 6); break }
                       print \$1, \$3, pt }' > ${prefix}.reads_pt.tsv

    polya_features.py \\
        --reads ${prefix}.reads_pt.tsv --features ${prefix}.read_features.tsv \\
        --catalogue ${catalogue} --classification ${classification} \\
        --mito-contig ${params.polya_mito_contig} \\
        --cap-nuclear ${params.polya_cap_nuclear} --cap-mito ${params.polya_cap_mito} \\
        --min-reads ${params.polya_min_reads} \\
        --prefix ${prefix}
    """

    // Stub test for pipeline integration
    stub:
    prefix = meta.id
    """
    touch ${prefix}.isoform_polya.bed ${prefix}.gene_polya.bed ${prefix}.polya_hist.tsv ${prefix}.polya_summary.tsv
    echo stub | gzip > ${prefix}.polya_features.tsv.gz
    echo stub | gzip > ${prefix}.polya_reads.tsv.gz
    """
}
