// ============================================================
// 02d_assembly_checkpoint.nf
// Measure an assembly at a named point in the pipeline.
// gfastats + BUSCO; the caller aliases this workflow per checkpoint.
// ============================================================

nextflow.enable.dsl = 2

include { BUSCO_ASSEMBLY; GFASTATS_QC } from '../modules/02_multi_assembler.nf'

workflow ASSEMBLY_CHECKPOINT {

    take:
        assembly
        checkpoint      // name used in every output filename, e.g. 'purged'
        busco_lineage
        busco_db_dir

    main:
        def named = assembly.map { p -> [checkpoint, p] }

        BUSCO_ASSEMBLY(named, channel.value(busco_lineage), busco_db_dir)
        GFASTATS_QC(named)

    emit:
        busco_table   = BUSCO_ASSEMBLY.out.busco_table.map   { n, p -> p }
        busco_summary = BUSCO_ASSEMBLY.out.busco_summary.map { n, p -> p }
        busco_json    = BUSCO_ASSEMBLY.out.busco_summary_json.map { n, p -> p }
        gfastats      = GFASTATS_QC.out.gfastats_report.map  { n, p -> p }
        versions      = BUSCO_ASSEMBLY.out.versions.mix(GFASTATS_QC.out.versions)
}
