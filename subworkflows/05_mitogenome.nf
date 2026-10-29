// ============================================================
// 05_mitogenome.nf
// Mitochondrial genome extraction by mapping assembly reads
// to a reference mitogenome
// ============================================================

nextflow.enable.dsl = 2

include { MITO_READ_FILTER ; MITO_MAPPING ; MITO_COVERAGE_STATS ; MITOHIFI } from '../modules/05_mitogenome.nf'
include { GETORGANELLE_DB ; MITO_GETORGANELLE ; MITO_MAPPED_READS ; MITO_READS_HQ ; MITO_MEDAKA } from '../modules/05_mitogenome.nf'
include { MITO_ROTATE ; MITO_LIFTOFF ; MITO_MITOS2 ; MITO_GENE_CHECK } from '../modules/05_mitogenome.nf'
include { FCS_ADAPTOR as MITO_FCS_ADAPTOR ; VECSCREEN as MITO_VECSCREEN ; VECSCREEN_CLASSIFY as MITO_VECSCREEN_CLASSIFY ;
          SCREEN_REVIEW as MITO_SCREEN_REVIEW ; APPLY_SCREEN_ACTIONS as MITO_APPLY_SCREEN } from '../modules/03c_screens.nf'

// ============================================================
// Workflow: MITOGENOME_PIPELINE
// Reads filter → map to reference → coverage stats;
// pre-purge assembly → MitoHiFi (extract, circularise).
// Final mitogenome: GetOrganelle on Illumina DNA when supplied, otherwise the
// MitoHiFi contig polished with high-quality mito ONT reads (Medaka). Either
// way: rotate to the annotation reference, liftover + MITOS2, protein check.
// ============================================================
workflow MITOGENOME_PIPELINE {

    take:
        assembly_reads
        mito_ref
        assembly        // pre-purge assembly (purge_dups drops the mito contig)
        model_reads     // original ONT reads (headers carry the basecaller model for Medaka)
        illumina_r1     // list of Illumina DNA R1 files, or [] -> ONT route
        illumina_r2     // matching R2 files, same order
        illumina_names  // their original file names, for the log
        annot_ref       // annotated related mitogenome (FASTA), or [] -> no liftover
        annot_ref_gff   // its GFF3

    main:
        MITO_READ_FILTER(assembly_reads)
        MITO_MAPPING(MITO_READ_FILTER.out.mito_reads, mito_ref)
        MITO_COVERAGE_STATS(MITO_MAPPING.out.mito_bam, MITO_MAPPING.out.mito_bai)
        MITOHIFI(assembly, mito_ref)

        def versions_ch = MITO_READ_FILTER.out.versions.mix(MITO_MAPPING.out.versions, MITO_COVERAGE_STATS.out.versions, MITOHIFI.out.versions)
        def report_ch   = MITOHIFI.out.stats
        def mito_seq_ch

        if (params.illumina_dna_r1) {
            GETORGANELLE_DB()
            MITO_GETORGANELLE(illumina_r1, illumina_r2, illumina_names, mito_ref, GETORGANELLE_DB.out.db)
            mito_seq_ch = MITO_GETORGANELLE.out.fasta.map { f -> ['Illumina DNA, GetOrganelle', f] }
            report_ch   = report_ch.mix(MITO_GETORGANELLE.out.log)
            versions_ch = versions_ch.mix(MITO_GETORGANELLE.out.versions)
        } else {
            MITO_MAPPED_READS(MITO_MAPPING.out.mito_bam, MITO_MAPPING.out.mito_bai)
            MITO_READS_HQ(MITO_MAPPED_READS.out.reads)
            MITO_MEDAKA(MITOHIFI.out.fasta, MITO_READS_HQ.out.reads, model_reads)
            mito_seq_ch = MITO_MEDAKA.out.fasta.map { f -> ['ONT contig (MitoHiFi), Medaka-polished', f] }
            report_ch   = report_ch.mix(MITO_MEDAKA.out.stats)
            versions_ch = versions_ch.mix(MITO_MAPPED_READS.out.versions, MITO_READS_HQ.out.versions, MITO_MEDAKA.out.versions)
        }

        // Adaptor/vector screens before annotation, as for the nuclear contigs (03c). No FCS-GX:
        // it is a genome-scale taxonomy screen and is not documented for organelle sequences.
        if (params.run_screens.toString().toBoolean()) {
            def m_in = mito_seq_ch.map { _src, f -> ['mitogenome', f] }
            MITO_FCS_ADAPTOR(m_in, file(params.fcs_adaptor_image, checkIfExists: true))
            MITO_VECSCREEN(m_in, file(params.univec, checkIfExists: true))
            MITO_VECSCREEN_CLASSIFY(MITO_VECSCREEN.out.raw)
            MITO_SCREEN_REVIEW(
                MITO_VECSCREEN_CLASSIFY.out.regions.map { s, v -> [s, file("${projectDir}/assets/NO_FCS_GX_REPORT", checkIfExists: true), v] },
                params.screen_decisions ? file(params.screen_decisions, checkIfExists: true)
                                        : file("${projectDir}/assets/NO_SCREEN_DECISIONS", checkIfExists: true))
            MITO_APPLY_SCREEN(m_in.join(MITO_FCS_ADAPTOR.out.report).join(MITO_SCREEN_REVIEW.out.actions))
            mito_seq_ch = mito_seq_ch.combine(MITO_APPLY_SCREEN.out.fasta).map { src, _raw, screened -> [src, screened] }
            versions_ch = versions_ch.mix(MITO_FCS_ADAPTOR.out.versions, MITO_VECSCREEN.out.versions)
        }

        def mitogenome_ch = mito_seq_ch.map { src, f -> f }
        def mito_gff_ch   = channel.empty()
        def mito_check_ch = channel.empty()
        if (params.mito_annot_ref_fasta && !(params.mito_annot_ref_fasta.toString() in ['false', 'null', ''])) {
            MITO_ROTATE(mito_seq_ch, annot_ref)
            mitogenome_ch = MITO_ROTATE.out.fasta
            MITO_LIFTOFF(mitogenome_ch, annot_ref, annot_ref_gff)
            MITO_GENE_CHECK(mitogenome_ch, MITO_LIFTOFF.out.gff, annot_ref, annot_ref_gff)
            report_ch   = report_ch.mix(MITO_ROTATE.out.info, MITO_LIFTOFF.out.gff, MITO_LIFTOFF.out.unmapped, MITO_GENE_CHECK.out.tsv)
            mito_gff_ch   = MITO_LIFTOFF.out.gff
            mito_check_ch = MITO_GENE_CHECK.out.tsv
            versions_ch = versions_ch.mix(MITO_ROTATE.out.versions, MITO_LIFTOFF.out.versions, MITO_GENE_CHECK.out.versions)
        }
        MITO_MITOS2(mitogenome_ch)
        report_ch   = report_ch.mix(MITO_MITOS2.out.dir)
        versions_ch = versions_ch.mix(MITO_MITOS2.out.versions)

    emit:
        mito_bam    = MITO_MAPPING.out.mito_bam
        mito_bai    = MITO_MAPPING.out.mito_bai
        coverage    = MITO_COVERAGE_STATS.out.coverage
        flagstat    = MITO_COVERAGE_STATS.out.flagstat
        mitogenome  = mitogenome_ch
        mito_gff    = mito_gff_ch      // Liftoff annotation (empty without mito_annot_ref_fasta)
        mito_check  = mito_check_ch    // per-gene check against the reference
        // outputs with no MultiQC module -> custom content
        report_files = report_ch
        versions    = versions_ch
}
