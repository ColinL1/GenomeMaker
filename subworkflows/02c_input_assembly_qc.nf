// ============================================================
// 02c_input_assembly_qc.nf (subworkflow)
// QC of a supplied assembly (params.input_assembly) with the same channel
// shapes MULTI_ASSEMBLER produces.
// ============================================================

nextflow.enable.dsl = 2

include { BUSCO_ASSEMBLY; GFASTATS_QC; COMBINE_COMPARISON } from '../modules/02_multi_assembler.nf'

workflow INPUT_ASSEMBLY_QC {

    take:
        assembly
        busco_lineage
        busco_db_dir

    main:
        def name = 'input_assembly'
        def assembly_ch = assembly.map { p -> [name, p] }

        BUSCO_ASSEMBLY(assembly_ch, channel.value(busco_lineage), busco_db_dir)
        GFASTATS_QC(assembly_ch)

        COMBINE_COMPARISON(
            GFASTATS_QC.out.gfastats_report.map { n, p -> p }.collect(),
            BUSCO_ASSEMBLY.out.busco_summary.map { n, p -> p }.collect()
        )

    emit:
        assemblies       = assembly_ch                          
        selected_assembly = assembly                            
        selected_name     = channel.value(name)
        gfastats          = GFASTATS_QC.out.gfastats_report       
        busco_tables      = BUSCO_ASSEMBLY.out.busco_table        
        busco_summary     = BUSCO_ASSEMBLY.out.busco_summary      
        comparison        = COMBINE_COMPARISON.out.comparison_tsv
        comparison_json   = COMBINE_COMPARISON.out.comparison_json
        versions        = BUSCO_ASSEMBLY.out.versions.mix(GFASTATS_QC.out.versions, COMBINE_COMPARISON.out.versions)
}
