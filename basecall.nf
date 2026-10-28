#!/usr/bin/env nextflow
// Stage 1: pod5 -> dorado uBAM -> alignment (Parabricks or CPU minimap2) -> mod tags
//
// Example Usage:
//   nextflow run basecall.nf -profile aicr,dna     --input runs.csv --outdir results/run1
//   nextflow run basecall.nf -profile explorer,rna --input runs.csv --outdir results/run1
//
// Samplesheet (CSV, header required):
//   sample_id  biological sample; becomes SM. Rows sharing it are merged.
//   run_id     one flowcell / pod5 set; becomes RG ID (default: sample_id). Unique.
//   pod5       a .pod5 file or a directory searched recursively
//   ubam       instead of pod5: an already basecalled, unaligned dorado BAM (e.g. a basecall
//              that finished after its run was cut off). It skips DORADO_BASECALLER and joins
//              the basecaller's output, so alignment, tag transfer, merge and QC treat it the
//              same. Each row gives exactly one of pod5 or ubam.
//   dorado_version  required with ubam: the VN of its @PG basecaller line
//              (samtools view -H <ubam> | grep -m1 '^@PG.*ID:basecaller'). Its basecall model
//              must be params.dorado_model, which the Stage 2 samplesheet records.
//   ...        any other columns (condition, replicate, ...) ride along to Stage 2

include { DORADO_BASECALLER }             from './modules/dorado/basecaller.nf'
include { REFERENCE_INDEX }               from './modules/align/build_index.nf'
include { MINIMAP2 }                      from './modules/align/minimap2.nf'
include { PBRUN_MINIMAP2 }                from './modules/align/pbrun_minimap2.nf'
include { DONOR_NAMESORT; COPY_TAGS }    from './modules/bam/tag_transfer.nf'
include { MERGE_RUNS; BAM_QC }            from './modules/bam/merge_qc.nf'

workflow {
    main:
    def align_to   = asList(params.align_to)
    def mod_models = asList(params.mod_models)
    validateParams(align_to, mod_models)
    def rows = parseSamplesheet(file(params.input, checkIfExists: true))
    def n_ubam = rows.count { r -> r[0].source == 'ubam' }
    log.info "basecall.nf: ${params.data_type}, ${rows.size()} run(s) (${rows.size() - n_ubam} to basecall, ${n_ubam} from a uBAM), aligner=${params.aligner}, " +
             "references=${align_to.join(',')}, mods=${mod_models.join(',') ?: 'none'}"

    // mv only exists when dorado emitted moves; asking for it otherwise fails the check.
    def carry = params.carry_tags.tokenize(',').findAll { t -> t != 'mv' || params.emit_moves }.join(',')

    // references: FASTA + cached .fai and .mmi
    ch_fasta = channel.fromList(align_to.collect { name -> tuple(name, file(params.genomes[name], checkIfExists: true)) })
    ch_ref = ch_fasta.join(REFERENCE_INDEX(ch_fasta, params.mm2_index_args))   // (ref_name, fasta, fai, mmi)

    // basecall the pod5 rows; uBAM rows (basecalled already) join the basecaller's output.
    // source only routes a row: it is dropped from meta first, so pod5 tasks hash as before.
    ch_rows = channel.fromList(rows).branch { meta, _input ->
        ubam: meta.source == 'ubam'
        pod5: true
    }
    DORADO_BASECALLER(ch_rows.pod5.map { meta, pod5 -> tuple(dropSource(meta), pod5) }, mod_models)
    ch_basecalled = DORADO_BASECALLER.out.ubam.map { meta, ubam, version -> tuple(meta + [dorado_version: version], ubam) }
    ch_ubam = ch_basecalled.mix(ch_rows.ubam.map { meta, ubam -> tuple(dropSource(meta), ubam) })

    // align + restore tags
    ch_tagged = channel.empty()
    ch_logs   = channel.empty()
    if( params.aligner == 'minimap2' ) {
        ch_in = ch_ubam.combine(ch_ref)
            .map { meta, ubam, ref_name, _fasta, _fai, mmi -> tuple(meta, ubam, ref_name, mmi, bamPrefix(meta, ref_name)) }
        ch_tagged = MINIMAP2(ch_in, carry).bam
    }
    else {
        // Tags come back from the uBAM, name-sorted once per run and paired with
        // each of that run's alignments by run_id.
        PBRUN_MINIMAP2(ch_ubam.combine(ch_ref))
        ch_logs = PBRUN_MINIMAP2.out.log
        ch_donor = DONOR_NAMESORT(ch_ubam, carry).map { meta, donor -> tuple(meta.run_id, donor) }
        ch_in = PBRUN_MINIMAP2.out.bam
            .map { meta, ref_name, bam, bai -> tuple(meta.run_id, meta, ref_name, bam, bai) }
            .combine(ch_donor, by: 0)
            .map { _run_id, meta, ref_name, bam, bai, donor -> tuple(meta, ref_name, bam, bai, donor, bamPrefix(meta, ref_name)) }
        ch_tagged = COPY_TAGS(ch_in, carry).bam
    }

    // merge flowcells per sample
    ch_split = ch_tagged.branch { meta, _ref_name, _bam, _bai ->
        single: meta.n_runs == 1
        multi:  true
    }
    ch_merged = MERGE_RUNS(
        ch_split.multi
            .map { meta, ref_name, bam, bai -> tuple(groupKey("${meta.id}\t${ref_name}", meta.n_runs), sampleMeta(meta), ref_name, bam, bai) }
            .groupTuple()
            .map { _key, metas, ref_names, bams, bais -> tuple(metas[0], ref_names[0], bams, bais) }
    )
    ch_final = ch_split.single
        .map { meta, ref_name, bam, bai -> tuple(sampleMeta(meta), ref_name, bam, bai) }
        .mix(ch_merged)

    ch_qc = BAM_QC(ch_final)

    // hand-off to Stage 2
    // Paths are relative to outdir, so the sheet still resolves after rsync.
    def extra = rows.collectMany { r -> r[0].extra.keySet() as List }.unique()
    def header = (['sample_id', 'reference', 'data_type', 'bam', 'bai', 'basecall_model',
                   'mod_models', 'dorado_version', 'n_runs'] + extra).join(',')
    ch_samplesheet = ch_final
        .map { meta, ref_name, bam, bai ->
            ([meta.id, ref_name, params.data_type, "bam/${ref_name}/${bam.name}", "bam/${ref_name}/${bai.name}",
              params.dorado_model, mod_models.join(' '), meta.dorado_version, meta.n_runs] +
             extra.collect { k -> meta.extra[k] ?: '' }).join(',')
        }
        .collectFile(name: 'analysis_samplesheet.csv', seed: header, newLine: true, sort: true)

    // `md5sum -c checksums.md5` from the destination outdir verifies the transfer.
    ch_checksums = ch_qc
        .map { _meta, ref_name, _flagstat, md5 -> md5.text.trim().replaceFirst(/\s+\S+$/, "  bam/${ref_name}/${md5.name - ~/\.md5$/}") }
        .collectFile(name: 'checksums.md5', newLine: true, sort: true)

    publish:
    ubam        = params.publish_ubam ? ch_ubam : channel.empty()
    bam         = ch_final
    qc          = ch_qc
    samplesheet = ch_samplesheet
    checksums   = ch_checksums
    logs        = ch_logs
}

output {
    ubam {
        path 'ubam'
    }
    bam {
        path { _meta, ref_name, _bam, _bai -> "bam/${ref_name}" }
    }
    qc {
        path { _meta, ref_name, _flagstat, _md5 -> "qc/${ref_name}" }
    }
    samplesheet {
        path '.'
    }
    checksums {
        path '.'
    }
    logs {
        path 'logs'
    }
}

// Final BAM name: <sample>.<ref> when the sample has one run, else <run>.<ref>
// until MERGE_RUNS writes <sample>.<ref>.
def bamPrefix(meta, ref_name) {
    meta.n_runs == 1 ? "${meta.id}.${ref_name}" : "${meta.run_id}.${ref_name}"
}

// Drop the run-level keys once a BAM represents the whole sample.
def sampleMeta(meta) {
    meta.subMap(['id', 'n_runs', 'extra', 'dorado_version'])
}

// A row's input kind (pod5 / ubam) only routes it; it is not part of the run's metadata.
def dropSource(meta) {
    meta.findAll { k, _v -> k != 'source' }
}

def parseSamplesheet(sheet) {
    def records = sheet.splitCsv(header: true).collect { row -> row.collectEntries { k, v -> [(k?.trim()): v?.trim()] } }
    if( !records )
        error "Samplesheet ${sheet} has no rows"
    def keys = records[0].keySet()
    def missing = ['sample_id'].findAll { col -> !(col in keys) }
    if( !('pod5' in keys) && !('ubam' in keys) )
        missing += ['pod5 (or ubam)']
    if( missing )
        error "Samplesheet ${sheet} is missing column(s): ${missing.join(', ')}"
    records.eachWithIndex { r, i ->
        if( (r.pod5 as boolean) == (r.ubam as boolean) )
            error "Samplesheet ${sheet} row ${i + 2} (${r.sample_id}): give exactly one of pod5 or ubam"
        if( r.ubam && !r.dorado_version )
            error "Samplesheet ${sheet} row ${i + 2} (${r.sample_id}): a ubam row needs dorado_version, " +
                  "the VN of its @PG basecaller line (samtools view -H ${r.ubam} | grep -m1 '^@PG.*ID:basecaller')"
    }

    def runIds = records.collect { r -> r.run_id ?: r.sample_id }
    def dups = runIds.countBy { id -> id }.findAll { _id, n -> n > 1 }.keySet()
    if( dups )
        error "run_id must be unique; repeated: ${dups.join(', ')} (set run_id per flowcell when a sample has several)"
    def nRuns = records.countBy { r -> r.sample_id }

    records.collect { r ->
        def meta = [
            id    : r.sample_id,
            run_id: r.run_id ?: r.sample_id,
            n_runs: nRuns[r.sample_id],
            extra : r.findAll { k, _v -> !(k in ['sample_id', 'run_id', 'pod5', 'ubam', 'dorado_version']) },
        ]
        if( r.ubam )
            return tuple(meta + [source: 'ubam', dorado_version: r.dorado_version], file(r.ubam, checkIfExists: true))
        return tuple(meta + [source: 'pod5'], file(r.pod5, checkIfExists: true))
    }
}

// List params arrive as strings when given on the command line (--align_to GRCh38).
def asList(value) {
    (value instanceof List) ? value : value.toString().tokenize(',').collect { v -> v.trim() }
}

def validateParams(align_to, mod_models) {
    if( !params.input )
        error("basecall.nf: " + "--input <samplesheet.csv> is required")
    if( !(params.data_type in ['DNA', 'RNA']) )
        error("basecall.nf: " + "pick a data type with -profile dna or -profile rna (data_type is '${params.data_type}')")
    if( !params.container_dir )
        error("basecall.nf: " + "pick a cluster with -profile explorer or -profile aicr")
    if( !(params.aligner in ['minimap2', 'parabricks']) )
        error("basecall.nf: " + "--aligner must be minimap2 or parabricks, not '${params.aligner}'")
    def cmods = mod_models.findAll { m -> m =~ /^(5mCG|5mC_|4mC)/ }
    if( cmods.size() > 1 )
        error("basecall.nf: " + "dorado allows one C modification model per run; got ${cmods.join(', ')}")
    align_to.each { name ->
        if( !params.genomes[name] )
            error("basecall.nf: " + "no FASTA for reference '${name}' in params.genomes for this cluster profile (set it, or narrow --align_to)")
    }
}
