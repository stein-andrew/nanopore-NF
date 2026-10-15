// Per-sample consensus, SALRR: drop low-abundance models from each tool (value >
// params.isoform_min_cpm; IsoQuant by CPM, StringTie by TPM, see bin/filter_models.py),
// merge the two with isomatch, and classify the result against the annotation.
// One process because the three steps take seconds next to IsoQuant and StringTie.
//
// isomatch names its sources after the input file names (S.isoquant.filtered.gtf,
// S.stringtie.filtered.gtf), so the track/present_absent tables say which tool found what.

process ISOFORM_CONSENSUS {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_small'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(meta), path(iq_gtf), path(iq_counts), path(st_gtf), path(fasta), path(fai), path(ref_gtf)

    output:
    // merged GTF goes on to the cohort merge (step 3b); everything goes to outdir
    tuple val(meta), path("${meta.id}.merged.gtf.gz"), emit: merged
    tuple val(meta),
          path("${meta.id}.{isoquant,stringtie}.filtered.{gtf,tsv}"),
          path("${meta.id}.{merged.gtf.gz,track.tsv.gz,present_absent.tsv.gz,merged_info.json,merged_params.json}"),
          path("${meta.id}.classify.*"), emit: out

    script:
    def id = meta.id
    """
    filter_models.py --gtf ${iq_gtf} --counts ${iq_counts} --min ${params.isoform_min_cpm} \\
        --out ${id}.isoquant.filtered.gtf --report ${id}.isoquant.filtered.tsv
    filter_models.py --gtf ${st_gtf} --attr TPM --min ${params.isoform_min_cpm} \\
        --out ${id}.stringtie.filtered.gtf --report ${id}.stringtie.filtered.tsv

    isomatch merge --ref-fa ${fasta} -o ${id} ${params.isomatch_merge_args} \\
        ${id}.isoquant.filtered.gtf ${id}.stringtie.filtered.gtf
    isomatch classify --ref-fa ${fasta} --ref-gtf ${ref_gtf} -o ${id}.classify ${params.isomatch_classify_args} \\
        ${id}.merged.gtf.gz
    """

    // Stub test for pipeline integration
    stub:
    def id = meta.id
    """
    touch ${id}.isoquant.filtered.gtf ${id}.isoquant.filtered.tsv ${id}.stringtie.filtered.gtf ${id}.stringtie.filtered.tsv
    touch ${id}.merged.gtf.gz ${id}.track.tsv.gz ${id}.present_absent.tsv.gz ${id}.merged_info.json ${id}.merged_params.json
    touch ${id}.classify.classification.txt.gz ${id}.classify.annotated.gtf.gz
    """
}
