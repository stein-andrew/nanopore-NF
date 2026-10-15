#!/usr/bin/env nextflow
// Stage 2 (RNA): Stage 1 BAMs -> SALRR isoform curation.
//
//   per sample   MAPQ >= 40 -> IsoQuant + StringTie -> CPM/TPM > 1 -> isomatch merge -> classify
//   (step 3b, next: cohort isomatch merge -> IsoQuant requant -> median CPM > 1)
//
// Example Usage:
//   nextflow run RNA_analysis.nf -profile explorer,rna -params-file params/run1_rna.yaml
//   nextflow run RNA_analysis.nf -profile explorer,rna,stub -stub-run -params-file params/run1_rna.yaml \
//       --outdir /scratch/$USER/ont-pipeline/stub-out-rna
//
// Input: Stage 1's analysis_samplesheet.csv as written. bam/bai are relative to the sheet's
// own directory (Stage 1's outdir), so this works from any launch directory.
//   sample_id, reference, data_type, bam, bai, basecall_model, mod_models, dorado_version,
//   n_runs, ...extra columns (kept in meta.extra)
// Rows whose data_type is not RNA are skipped with a warning.

include { FAIDX }             from './modules/ref/faidx.nf'
include { FILTER_MAPQ }       from './modules/bam/filter_mapq.nf'
include { ISOQUANT_DISCOVER } from './modules/isoform/isoquant.nf'
include { STRINGTIE }         from './modules/isoform/stringtie.nf'
include { ISOFORM_CONSENSUS } from './modules/isoform/consensus.nf'

workflow {
    main:
    if( !params.input )
        error("RNA_analysis.nf: " + "input (Stage 1's analysis_samplesheet.csv) is required")
    def rows = parseSamplesheet(file(params.input, checkIfExists: true))
    def refs = rows.collect { row -> row[0].ref }.unique()
    validateParams(refs)
    log.info "RNA_analysis.nf: ${rows.size()} BAM(s), references=${refs.join(',')}, " +
             "MAPQ >= ${params.isoform_min_mapq}, per-sample model filter > ${params.isoform_min_cpm}"

    // references: (ref, fasta, fai) and (ref, gtf)
    ch_ref = FAIDX(channel.fromList(refs.collect { ref -> tuple(ref, file(params.genomes[ref], checkIfExists: true)) }))
    ch_gtf = channel.fromList(refs.collect { ref -> tuple(ref, file(params.annotations[ref], checkIfExists: true)) })

    FILTER_MAPQ(channel.fromList(rows))
    ch_keyed = FILTER_MAPQ.out.bam
        .map { meta, bam, bai -> tuple(meta.ref, meta, bam, bai) }
        .combine(ch_ref, by: 0)
        .combine(ch_gtf, by: 0)                                     // (ref, meta, bam, bai, fasta, fai, gtf)

    ISOQUANT_DISCOVER(ch_keyed.map { _ref, meta, bam, bai, fasta, fai, gtf -> tuple(meta, bam, bai, fasta, fai, gtf) })
    STRINGTIE(ch_keyed.map { _ref, meta, bam, bai, _fasta, _fai, gtf -> tuple(meta, bam, bai, gtf) })

    ISOFORM_CONSENSUS(
        ISOQUANT_DISCOVER.out.models
            .join(STRINGTIE.out)                                    // (meta, iq_gtf, iq_counts, st_gtf)
            .map { meta, iq_gtf, iq_counts, st_gtf -> tuple(meta.ref, meta, iq_gtf, iq_counts, st_gtf) }
            .combine(ch_ref, by: 0)
            .combine(ch_gtf, by: 0)
            .map { _ref, meta, iq_gtf, iq_counts, st_gtf, fasta, fai, gtf -> tuple(meta, iq_gtf, iq_counts, st_gtf, fasta, fai, gtf) }
    )

    publish:
    mapq_qc   = FILTER_MAPQ.out.flagstat
    isoquant  = ISOQUANT_DISCOVER.out.out
    stringtie = STRINGTIE.out
    consensus = ISOFORM_CONSENSUS.out.out
}

output {
    mapq_qc {
        path { meta, _flagstat -> "qc/${meta.ref}" }
    }
    isoquant {
        path { meta, _files -> "isoforms/${meta.ref}/isoquant/${meta.id}" }
    }
    stringtie {
        path { meta, _gtf -> "isoforms/${meta.ref}/stringtie" }
    }
    consensus {
        path { meta, _filtered, _merged, _classify -> "isoforms/${meta.ref}/consensus/${meta.id}" }
    }
}

def parseSamplesheet(sheet) {
    def records = sheet.splitCsv(header: true).collect { row -> row.collectEntries { k, v -> [(k?.trim()): v?.trim()] } }
    if( !records )
        error "Samplesheet ${sheet} has no rows"
    def missing = ['sample_id', 'reference', 'data_type', 'bam', 'bai'].findAll { col -> !(col in records[0].keySet()) }
    if( missing )
        error "Samplesheet ${sheet} is missing column(s): ${missing.join(', ')}"

    def rna = records.findAll { r -> r.data_type == 'RNA' }
    if( rna.size() < records.size() )
        log.warn "RNA_analysis.nf: skipping ${records.size() - rna.size()} non-RNA row(s) in ${sheet}"
    if( !rna )
        error "Samplesheet ${sheet} has no RNA rows"
    def dups = rna.countBy { r -> "${r.sample_id} on ${r.reference}".toString() }.findAll { _k, n -> n > 1 }.keySet()
    if( dups )
        error "Samplesheet ${sheet} repeats sample/reference pair(s): ${dups.join(', ')}"

    def known = ['sample_id', 'reference', 'data_type', 'bam', 'bai', 'basecall_model', 'mod_models', 'dorado_version', 'n_runs']
    rna.collect { r ->
        def meta = [
            id            : r.sample_id,
            ref           : r.reference,
            basecall_model: r.basecall_model,
            mod_models    : r.mod_models,
            dorado_version: r.dorado_version,
            n_runs        : r.n_runs,
            extra         : r.findAll { k, _v -> !(k in known) },
        ]
        tuple(meta, sheetPath(sheet, r.bam), sheetPath(sheet, r.bai))
    }
}

// Stage 1 writes paths relative to its outdir, where the sheet lives. Absolute paths pass.
def sheetPath(sheet, p) {
    if( !p )
        error "Samplesheet ${sheet} has an empty bam/bai cell"
    file(p.startsWith('/') ? p : "${sheet.parent}/${p}", checkIfExists: true)
}

def validateParams(refs) {
    if( params.data_type != 'RNA' )
        error("RNA_analysis.nf: " + "use -profile rna (data_type is '${params.data_type}')")
    if( !params.container_dir )
        error("RNA_analysis.nf: " + "pick a cluster with -profile explorer")
    if( !params.isoform_image )
        error("RNA_analysis.nf: " + "isoform_image is not set (conf/rna.config)")
    refs.each { ref ->
        if( !params.genomes[ref] )
            error("RNA_analysis.nf: " + "no FASTA for reference '${ref}' in params.genomes for this cluster profile")
        if( !params.annotations[ref] )
            error("RNA_analysis.nf: " + "no GTF for reference '${ref}' in params.annotations for this cluster profile")
    }
}
