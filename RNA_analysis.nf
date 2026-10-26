#!/usr/bin/env nextflow
// Stage 2 (RNA): Stage 1 BAMs -> SALRR isoform curation.
//
//   per sample   MAPQ >= 40 -> IsoQuant + StringTie -> CPM/TPM > 1 -> isomatch merge -> classify
//   cohort       isomatch merge of all samples -> classify -> catalogue
//   per sample   IsoQuant requantification against the catalogue
//   cohort       median CPM across samples > 1 -> high-confidence transcripts
//
//   per sample   modkit sample-probs -> per-base thresholds -> modkit pileup (raw)
//                -> Nvalid_cov >= 20 -> percent >= 20; gIVT FP correction -> corrected >= 20
//   per sample   at the final-stage sites: modification levels per isoform and per gene
//   per sample   poly(A) tail lengths per read, isoform and gene (chrM apart, SALRR's exclusions)
//   per sample   Clair3-RNA small variants (phasing model, REDIportal editing tags)
//                -> known / novel / editing against gnomAD v4.1 (chrM: v3.1), with ClinVar
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
include { ISOQUANT_DISCOVER; ISOQUANT_REQUANT } from './modules/isoform/isoquant.nf'
include { STRINGTIE }                           from './modules/isoform/stringtie.nf'
include { ISOFORM_CONSENSUS }                   from './modules/isoform/consensus.nf'
include { COHORT_CATALOGUE; COHORT_FILTER }     from './modules/isoform/cohort.nf'
include { MODKIT_PILEUP; MOD_STAGES }           from './modules/mods/modkit.nf'
include { MOD_FEATURES }                        from './modules/mods/features.nf'
include { POLYA_FEATURES }                      from './modules/polya/polya.nf'
include { CLAIR3_RNA }                          from './modules/variants/clair3_rna.nf'
include { KNOWN_NOVEL }                         from './modules/variants/known_novel.nf'
include { REGION_SETS }                         from './modules/ref/region_sets.nf'
include { READ_QC }                             from './modules/qc/read_qc.nf'
include { COHORT_MATRICES }                     from './modules/cohort/matrices.nf'

workflow {
    main:
    if( !params.input )
        error("RNA_analysis.nf: " + "input (Stage 1's analysis_samplesheet.csv) is required")
    def rows = parseSamplesheet(file(params.input, checkIfExists: true))
    def refs = rows.collect { row -> row[0].ref }.unique()
    def nPerRef = rows.countBy { row -> row[0].ref }
    validateParams(refs)
    log.info "RNA_analysis.nf: ${rows.size()} BAM(s), references=${refs.join(',')}, " +
             "MAPQ >= ${params.isoform_min_mapq}, per-sample model filter > ${params.isoform_min_cpm}, " +
             "cohort median CPM > ${params.cohort_min_median_cpm}; " +
             "mod thresholds ${params.mod_filter_thresholds ?: 'from sample-probs'}, " +
             "stages cov >= ${params.mods_min_coverage}, pct >= ${params.mods_min_percent}"

    // references: (ref, fasta, fai) and (ref, gtf)
    ch_ref = FAIDX(channel.fromList(refs.collect { ref -> tuple(ref, file(params.genomes[ref], checkIfExists: true)) }))
    ch_gtf = channel.fromList(refs.collect { ref -> tuple(ref, file(params.annotations[ref], checkIfExists: true)) })

    // globin genes and the variant mask (IG / TR / HLA), from the annotation
    REGION_SETS(ch_gtf)
    ch_globin = REGION_SETS.out.regions.map { ref, globin, _mask -> tuple(ref, globin) }
    ch_mask   = REGION_SETS.out.regions.map { ref, _globin, mask -> tuple(ref, mask) }

    // per-sample read composition (globin / mito share), on the Stage 1 BAM
    READ_QC(
        channel.fromList(rows)
            .map { meta, bam, bai -> tuple(meta.ref, meta, bam, bai) }
            .combine(ch_globin, by: 0)
            .map { _ref, meta, bam, bai, globin -> tuple(meta, bam, bai, globin) }
    )

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

    // ---- cohort: one catalogue per reference from every sample's consensus
    // groupKey with the sample count lets each reference's group close as soon as its
    // last sample arrives, instead of waiting for every reference.
    COHORT_CATALOGUE(
        ISOFORM_CONSENSUS.out.merged
            .map { meta, gtf -> tuple(groupKey(meta.ref, nPerRef[meta.ref]), gtf) }
            .groupTuple()
            // sorted by file name (= sample ID): isomatch numbers its inputs S1..Sn in
            // command-line order, which must not depend on which sample finished first
            .map { key, gtfs -> tuple(key.toString(), gtfs.sort { f -> f.name }) }
            .combine(ch_ref, by: 0)
            .combine(ch_gtf, by: 0)                                 // (ref, gtfs, fasta, fai, gtf)
    )

    // ---- every sample's MAPQ-filtered reads, requantified against its reference's catalogue
    ISOQUANT_REQUANT(
        ch_keyed
            .combine(COHORT_CATALOGUE.out.catalogue, by: 0)         // (ref, meta, bam, bai, fasta, fai, gtf, catalogue)
            .map { _ref, meta, bam, bai, fasta, fai, _gtf, catalogue -> tuple(meta, bam, bai, fasta, fai, catalogue) }
    )

    // ---- median CPM across samples -> high-confidence transcripts
    COHORT_FILTER(
        ISOQUANT_REQUANT.out.counts
            .map { meta, counts -> tuple(groupKey(meta.ref, nPerRef[meta.ref]), meta.id, counts) }
            .groupTuple()
            .map { key, ids, counts -> tuple(key.toString(), ids, counts) }
            .combine(COHORT_CATALOGUE.out.catalogue, by: 0)
            .combine(COHORT_CATALOGUE.out.classification, by: 0)    // (ref, ids, counts, catalogue, classification)
    )

    // ---- modifications, on the Stage 1 BAM: every primary alignment, not the MAPQ >= 40
    // subset that isoform curation uses
    ch_errtab = channel.value(file(params.givt_error_table, checkIfExists: true))
    MODKIT_PILEUP(channel.fromList(rows))
    MOD_STAGES(
        MODKIT_PILEUP.out.pileup
            .map { meta, bed, flags -> tuple(meta.ref, meta, bed, flags) }
            .combine(ch_ref, by: 0)
            .map { _ref, meta, bed, flags, fasta, fai -> tuple(meta, bed, flags, fasta, fai) },
        ch_errtab
    )

    // ---- modification levels per isoform and per gene, at the final-stage sites
    MOD_FEATURES(
        channel.fromList(rows)                                      // (meta, Stage 1 bam, bai)
            .join(MOD_STAGES.out.fpcorr)
            .join(MODKIT_PILEUP.out.flags)
            .join(ISOQUANT_REQUANT.out.assignments)                 // (meta, bam, bai, fpcorr, modkit_flags, givt_flags, assignments)
            .map { meta, bam, bai, fpcorr, mflags, gflags, assign -> tuple(meta.ref, meta, bam, bai, fpcorr, mflags, gflags, assign) }
            .combine(COHORT_CATALOGUE.out.catalogue, by: 0)
            .combine(COHORT_CATALOGUE.out.classification, by: 0)
            .combine(ch_ref, by: 0)
            .combine(ch_globin, by: 0)
            .map { _ref, meta, bam, bai, fpcorr, mflags, gflags, assign, catalogue, classification, fasta, fai, globin ->
                tuple(meta, bam, bai, fpcorr, mflags, gflags, assign, catalogue, classification, fasta, fai, globin) },
        ch_errtab
    )

    // ---- poly(A) tails per read, isoform and gene
    POLYA_FEATURES(
        channel.fromList(rows)                                      // (meta, Stage 1 bam, bai)
            .join(ISOQUANT_REQUANT.out.assignments)
            .map { meta, bam, bai, assign -> tuple(meta.ref, meta, bam, bai, assign) }
            .combine(COHORT_CATALOGUE.out.catalogue, by: 0)
            .combine(COHORT_CATALOGUE.out.classification, by: 0)
            .map { _ref, meta, bam, bai, assign, catalogue, classification ->
                tuple(meta, bam, bai, assign, catalogue, classification) }
    )

    // ---- small variants with Clair3-RNA, on the Stage 1 BAM
    CLAIR3_RNA(
        channel.fromList(rows)
            .map { meta, bam, bai -> tuple(meta.ref, meta, bam, bai) }
            .combine(ch_ref, by: 0)
            .map { _ref, meta, bam, bai, fasta, fai -> tuple(meta, bam, bai, fasta, fai) }
    )

    // ---- known vs novel; the variant references are GRCh38, so other references are skipped
    ch_varrefs = channel.value(tuple(
        file(params.gnomad_nuclear), file("${params.gnomad_nuclear}.csi"),
        file(params.gnomad_mito),    file("${params.gnomad_mito}.csi"),
        file(params.clinvar),        file("${params.clinvar}.csi")))
    KNOWN_NOVEL(
        CLAIR3_RNA.out.vcf
            .filter { meta, _vcf -> meta.ref == params.variant_refs_genome }
            .map { meta, vcf -> tuple(meta.ref, meta, vcf) }
            .combine(ch_ref, by: 0)
            .combine(ch_mask, by: 0)
            .map { _ref, meta, vcf, fasta, fai, mask -> tuple(meta, vcf, fasta, fai, mask) },
        ch_varrefs
    )

    // ---- cohort matrices: one row per gene / site / feature, one column per sample
    COHORT_MATRICES(
        ISOQUANT_REQUANT.out.gene_counts
            .join(MOD_STAGES.out.fpcorr)
            .join(POLYA_FEATURES.out.features)
            .join(READ_QC.out)
            .map { meta, genes, fpcorr, polya, qc -> tuple(groupKey(meta.ref, nPerRef[meta.ref]), meta.id, genes, fpcorr, polya, qc) }
            .groupTuple()
            .map { key, ids, genes, fpcorrs, polyas, qcs ->
                def s = [ids, genes, fpcorrs, polyas, qcs].transpose().sort { r -> r[0] }   // sample-id order
                tuple(key.toString(), s.collect { r -> r[0] }, s.collect { r -> r[1] }, s.collect { r -> r[2] },
                      s.collect { r -> r[3] }, s.collect { r -> r[4] }) }
            .combine(COHORT_CATALOGUE.out.catalogue, by: 0)
            .combine(COHORT_CATALOGUE.out.classification, by: 0)
    )

    publish:
    mapq_qc   = FILTER_MAPQ.out.flagstat
    isoquant  = ISOQUANT_DISCOVER.out.out
    stringtie = STRINGTIE.out
    consensus = ISOFORM_CONSENSUS.out.out
    cohort    = COHORT_CATALOGUE.out.out
    requant   = ISOQUANT_REQUANT.out.out
    high_conf = COHORT_FILTER.out
    mods_pileup = MODKIT_PILEUP.out.out
    mods_stages = MOD_STAGES.out.out
    mods_features = MOD_FEATURES.out
    polya = POLYA_FEATURES.out.out
    variants = CLAIR3_RNA.out.out
    variant_classes = KNOWN_NOVEL.out
    region_sets = REGION_SETS.out.out
    read_qc = READ_QC.out
    cohort_matrices = COHORT_MATRICES.out.out
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
    cohort {
        path { ref_name, _merged, _classify -> "isoforms/${ref_name}/cohort" }
    }
    requant {
        path { meta, _files -> "isoforms/${meta.ref}/requant/${meta.id}" }
    }
    high_conf {
        path { ref_name, _files -> "isoforms/${ref_name}/high_confidence" }
    }
    mods_pileup {
        path { meta, _files, _probs -> "mods/${meta.ref}/${meta.id}" }
    }
    mods_stages {
        path { meta, _files -> "mods/${meta.ref}/${meta.id}" }
    }
    mods_features {
        path { meta, _files -> "mods/${meta.ref}/${meta.id}" }
    }
    polya {
        path { meta, _files -> "polya/${meta.ref}/${meta.id}" }
    }
    variants {
        path { meta, _files -> "variants/${meta.ref}/${meta.id}" }
    }
    variant_classes {
        path { meta, _files -> "variants/${meta.ref}/${meta.id}" }
    }
    region_sets {
        path { ref_name, _files -> "qc/${ref_name}/regions" }
    }
    read_qc {
        path { meta, _tsv -> "qc/${meta.ref}" }
    }
    cohort_matrices {
        path { ref_name, _files -> "cohort/${ref_name}" }
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
    if( !params.mods_image )
        error("RNA_analysis.nf: " + "mods_image is not set (conf/rna.config)")
    if( !params.clair3_rna_image )
        error("RNA_analysis.nf: " + "clair3_rna_image is not set (conf/rna.config)")
    // the variant references must be the files themselves, each with its .csi index
    ['gnomad_nuclear', 'gnomad_mito', 'clinvar'].each { name ->
        def path = params[name]
        if( !path || !file(path).isFile() || !file("${path}.csi").isFile() )
            error("RNA_analysis.nf: " + "${name} must be an indexed BCF (with ${name}.csi beside it); got '${path}'")
    }
    if( !params.globin_genes )
        error("RNA_analysis.nf: " + "globin_genes is empty (conf/rna.config): the read QC needs at least one gene name")
    if( params.variant_mask_mhc && !(params.variant_mask_mhc ==~ /[^:]+:[0-9]+-[0-9]+/) )
        error("RNA_analysis.nf: " + "variant_mask_mhc must be chrom:start-end or '', not '${params.variant_mask_mhc}'")
    refs.findAll { r -> r != params.variant_refs_genome }.each { r ->
        log.warn("RNA_analysis.nf: known/novel annotation needs ${params.variant_refs_genome} references; skipping ${r} samples")
    }
    if( !params.givt_error_table || params.givt_error_table.toString().contains('CHANGE_ME') )
        error("RNA_analysis.nf: " + "givt_error_table is not set: point it at the gIVT 9mer table in the cluster profile")
    // exists isn't enough: a directory exists too, and gIVT only fails on it after the pileup
    if( !file(params.givt_error_table).isFile() )
        error("RNA_analysis.nf: " + "givt_error_table must be the table file itself, not a directory or a missing path: ${params.givt_error_table}")
    refs.each { ref ->
        if( !params.genomes[ref] )
            error("RNA_analysis.nf: " + "no FASTA for reference '${ref}' in params.genomes for this cluster profile")
        if( !params.annotations[ref] )
            error("RNA_analysis.nf: " + "no GTF for reference '${ref}' in params.annotations for this cluster profile")
    }
}
