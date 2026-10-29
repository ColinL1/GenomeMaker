// ============================================================
// 01_kmer_genome_size.nf
// K-mer counting, histogram, and GenomeScope2 genome-size estimation.
//
// Split into TWO workflows on purpose:
//
//   MERYL_READ_DB     builds the read k-mer database. Needed by GenomeScope2
//                     for genome sizing AND by Merqury for polishing QV, and
//                     Merqury is wanted even on the --input_assembly path where
//                     genome sizing is skipped entirely. Keeping the database
//                     build separate means it can be built once and consumed by
//                     either, instead of being trapped inside the sizing branch.
//
//   KMER_GENOME_SIZE  the GenomeScope2 model and the size parse, which only the
//                     de novo assembly path needs.
// ============================================================

nextflow.enable.dsl = 2

include { MERYL_COUNT; MERYL_HIST; GENOMESCOPE; PARSE_GENOME_SIZE } from '../modules/01_kmer_genome_size.nf'

// ============================================================
// Workflow: MERYL_READ_DB
// Read k-mer database + histogram. Consumed by GenomeScope2 (sizing) and
// Merqury (polishing QV).
// ============================================================
workflow MERYL_READ_DB {

    take:
        reads

    main:
        MERYL_COUNT(reads)
        MERYL_HIST(MERYL_COUNT.out.meryl_db)

    emit:
        meryl_db   = MERYL_COUNT.out.meryl_db
        meryl_hist = MERYL_HIST.out.meryl_hist
        versions   = MERYL_COUNT.out.versions.mix(MERYL_HIST.out.versions)
}

// ============================================================
// Workflow: KMER_GENOME_SIZE
// GenomeScope2 model → parsed haploid genome size for the assemblers.
// ============================================================
workflow KMER_GENOME_SIZE {

    take:
        meryl_hist

    main:
        GENOMESCOPE(meryl_hist)
        // The empty-summary check lives inside GENOMESCOPE itself, so a failure
        // aborts that task directly rather than via a channel subscriber.
        PARSE_GENOME_SIZE(GENOMESCOPE.out.genomescope_out)

    emit:
        genomescope_out = GENOMESCOPE.out.genomescope_out
        genome_size     = PARSE_GENOME_SIZE.out.genome_size
        versions        = GENOMESCOPE.out.versions.mix(PARSE_GENOME_SIZE.out.versions)
}
