// ============================================================
// 01_kmer_genome_size.nf
// K-mer counting, histogram, and GenomeScope 2.0 estimation
// ============================================================

nextflow.enable.dsl = 2

include { MERYL_COUNT; MERYL_HIST; GENOMESCOPE; PARSE_GENOME_SIZE } from '../modules/01_kmer_genome_size.nf'

// ============================================================
// Workflow: KMER_GENOME_SIZE
// Orchestrates meryl count → histogram → genomescope
// Emits estimated genome size for use by assemblers
// ============================================================
workflow KMER_GENOME_SIZE {

    take:
        reads_pc
        genomescope_script

    main:
        MERYL_COUNT(reads_pc)
        MERYL_HIST(MERYL_COUNT.out.meryl_db)
        GENOMESCOPE(MERYL_HIST.out.meryl_hist, genomescope_script) //TODO: fix genoscope silent fails? 
        PARSE_GENOME_SIZE(GENOMESCOPE.out.genomescope_out)
    emit:
        meryl_db        = MERYL_COUNT.out.meryl_db
        meryl_hist      = MERYL_HIST.out.meryl_hist
        genomescope_out = GENOMESCOPE.out.genomescope_out
        genome_size     = PARSE_GENOME_SIZE.out.genome_size   // <-- new, this is the val you need
}
