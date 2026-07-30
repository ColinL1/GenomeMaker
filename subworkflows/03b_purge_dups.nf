// ============================================================
// 03b_purge_dups.nf (subworkflow)
// ============================================================

nextflow.enable.dsl = 2

include { PURGE_DUPS_COVERAGE; PURGE_DUPS_RUN } from '../modules/03b_purge_dups.nf'

workflow PURGE_DUPS {

    take:
        assembly    // value channel, single path -- the SELECTED assembly
        reads       // ONT reads (assembly_reads)

    main:
        PURGE_DUPS_COVERAGE(assembly, reads)
        PURGE_DUPS_RUN(
            assembly,
            PURGE_DUPS_COVERAGE.out.base_cov,
            PURGE_DUPS_COVERAGE.out.cutoffs
        )

    emit:
        purged_assembly = PURGE_DUPS_RUN.out.purged_assembly
        haplotigs       = PURGE_DUPS_RUN.out.haplotigs
        dups_bed        = PURGE_DUPS_RUN.out.dups_bed
}
