// IsoQuant discovery, per sample, as SALRR runs it: -d nanopore with the annotation as
// --genedb and --complete_genedb (GENCODE has gene and transcript lines, so gffutils
// skips its slow inference). params.isoquant_args adds anything else.
//
// The converted annotation database goes to node-local tmp and is deleted afterwards:
// it is GB-scale for GENCODE, and SQLite is unreliable on network filesystems.
// IsoQuant writes into iq/<sample>/; its files are moved up to the task directory so
// they publish as isoquant/<sample>/<file>, not isoquant/<sample>/iq/<sample>/<file>.
//
// HOME is the task directory: IsoQuant 3.10 always creates ~/.config/IsoQuant (a cache of
// past GTF conversions) and the real home is read-only inside the container. A per-task
// home also keeps concurrent tasks off one cache file and keeps that cache, which would
// point at our deleted tmp database, out of later runs.

process ISOQUANT_DISCOVER {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_align'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(gtf)

    output:
    // Models + their counts feed the per-sample CPM filter. The counts file is
    // discovered_transcript_counts.tsv in the IsoQuant 3.x we have run (3.12.2) and
    // transcript_model_counts.tsv in older 3.x; exactly one exists.
    tuple val(meta), path("${meta.id}.transcript_models.gtf"),
          path("${meta.id}.{discovered_transcript_counts,transcript_model_counts}.tsv"), emit: models
    // Everything IsoQuant wrote, plus its log, for outdir. (Staged inputs never match.)
    tuple val(meta), path("${meta.id}.*"), emit: out

    script:
    def id = meta.id
    """
    export HOME=\$PWD/home
    mkdir -p "\$HOME"
    genedb=\$(mktemp -d "\${TMPDIR:-/tmp}/iq_genedb.XXXXXX")
    isoquant \\
        --data_type nanopore \\
        --reference ${fasta} \\
        --genedb ${gtf} --complete_genedb --genedb_output "\$genedb" \\
        --bam ${bam} \\
        --labels ${id} \\
        --prefix ${id} \\
        --threads ${task.cpus} \\
        ${params.isoquant_args} \\
        -o iq
    rm -rf "\$genedb" iq/${id}/aux

    mv iq/${id}/${id}.* .
    mv iq/isoquant.log ${id}.isoquant.log
    rm -rf iq

    if [[ ! -e ${id}.discovered_transcript_counts.tsv && ! -e ${id}.transcript_model_counts.tsv ]]; then
        echo "no IsoQuant model counts file; IsoQuant wrote:" >&2
        ls ${id}.* >&2
        exit 1
    fi
    """

    // Stub test for pipeline integration
    stub:
    def id = meta.id
    """
    touch ${id}.transcript_models.gtf ${id}.discovered_transcript_counts.tsv ${id}.transcript_counts.tsv ${id}.isoquant.log
    """
}

// IsoQuant requantification, per sample, against the cohort catalogue (SALRR: every sample's
// reads re-assigned to the unified catalogue). --no_model_construction: assign and count
// only, no new models. The catalogue has gene lines (gtf_add_genes.py), so
// --complete_genedb holds. read_assignments.tsv.gz in its output maps each read to a
// catalogue transcript: that is what isoform-level mod and poly(A) summaries use later.
// params.isoquant_requant_args (default --transcript_quantification unique_only, as SALRR)
// comes after params.isoquant_args, so it wins if both set the same flag.
// HOME and the tmp database: as in ISOQUANT_DISCOVER. bin/check_requant_ids.sh then fails
// the task unless every counted transcript is a catalogue ID (see that script for why).

process ISOQUANT_REQUANT {
    tag "${meta.id}:${meta.ref}"
    label 'cpu_align'
    container "${params.container_dir}/${params.isoform_image}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai), path(catalogue)

    output:
    tuple val(meta), path("${meta.id}.transcript_counts.tsv"), emit: counts
    tuple val(meta), path("${meta.id}.*"), emit: out

    script:
    def id = meta.id
    """
    export HOME=\$PWD/home
    mkdir -p "\$HOME"
    genedb=\$(mktemp -d "\${TMPDIR:-/tmp}/iq_genedb.XXXXXX")
    isoquant \\
        --data_type nanopore \\
        --reference ${fasta} \\
        --genedb ${catalogue} --complete_genedb --genedb_output "\$genedb" \\
        --no_model_construction \\
        --bam ${bam} \\
        --labels ${id} \\
        --prefix ${id} \\
        --threads ${task.cpus} \\
        ${params.isoquant_args} \\
        ${params.isoquant_requant_args} \\
        -o iq
    rm -rf "\$genedb" iq/${id}/aux

    mv iq/${id}/${id}.* .
    mv iq/isoquant.log ${id}.isoquant.log
    rm -rf iq

    check_requant_ids.sh ${catalogue} ${id}.transcript_counts.tsv
    """

    // Stub test for pipeline integration
    stub:
    def id = meta.id
    """
    touch ${id}.transcript_counts.tsv ${id}.read_assignments.tsv.gz ${id}.isoquant.log
    """
}
