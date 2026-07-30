// ============================================================
// 02_multi_assembler.nf (subworkflow)
// ============================================================

nextflow.enable.dsl = 2

include { ASSEMBLE_NECAT; ASSEMBLE_FLYE; ASSEMBLE_HIFIASM; ASSEMBLE_RAVEN } from '../modules/02_multi_assembler.nf'
include { BUSCO_ASSEMBLY; GFASTATS_QC; QUAST_QC; COMBINE_COMPARISON }        from '../modules/02_multi_assembler.nf'

workflow MULTI_ASSEMBLER {

    take:
        assembly_reads
        necat_config
        necat_read_list
        genome_size
        busco_lineage         // single lineage, e.g. anthozoa_odb12

    main:
        def assemblers_list = params.assemblers instanceof List
            ? params.assemblers
            : params.assemblers.split(',').collect { a -> a.trim() }

        def assembler_channels = []

        if ('necat' in assemblers_list) {
            if (!params.necat_config || !file(params.necat_config).exists())
                error "Assembler 'necat' requested but params.necat_config not found: ${params.necat_config}"
            if (!params.necat_read_list || !file(params.necat_read_list).exists())
                error "Assembler 'necat' requested but params.necat_read_list not found: ${params.necat_read_list}"
            ASSEMBLE_NECAT(assembly_reads, necat_config, necat_read_list)
            assembler_channels << ASSEMBLE_NECAT.out.assembly.map { p -> ['necat', p] }
        }
        if ('flye' in assemblers_list) {
            ASSEMBLE_FLYE(assembly_reads, genome_size)
            assembler_channels << ASSEMBLE_FLYE.out.assembly.map { p -> ['flye', p] }
        }
        if ('hifiasm' in assemblers_list) {
            ASSEMBLE_HIFIASM(assembly_reads, genome_size)
            assembler_channels << ASSEMBLE_HIFIASM.out.assembly.map { p -> ['hifiasm', p] }
        }
        if ('raven' in assemblers_list) {
            ASSEMBLE_RAVEN(assembly_reads)
            assembler_channels << ASSEMBLE_RAVEN.out.assembly.map { p -> ['raven', p] }
        }

        if (!assembler_channels)
            error "No valid assemblers selected. params.assemblers=${params.assemblers} (expected one or more of: necat, flye, hifiasm, raven)"

        def assemblies_ch = assembler_channels.inject { acc, ch -> acc.mix(ch) }

        BUSCO_ASSEMBLY(assemblies_ch, channel.value(busco_lineage))

        GFASTATS_QC(assemblies_ch)
        QUAST_QC(assemblies_ch)

        // Aggregation: name-embedded filenames make plain path lists safe here.
        COMBINE_COMPARISON(
            GFASTATS_QC.out.gfastats_report.map { name, p -> p }.collect(),
            BUSCO_ASSEMBLY.out.busco_summary.map { name, p -> p }.collect()
        )

    emit:
        assemblies      = assemblies_ch                          // (name, path)
        gfastats         = GFASTATS_QC.out.gfastats_report       // (name, path) -- name preserved
        busco_tables      = BUSCO_ASSEMBLY.out.busco_table         // (name, path) -- name preserved
        busco_summary      = BUSCO_ASSEMBLY.out.busco_summary        // (name, path) -- used by SELECT_BEST_ASSEMBLY
        quast_reports      = QUAST_QC.out.quast_report              // (name, path)
        comparison         = COMBINE_COMPARISON.out.comparison_tsv
        comparison_json    = COMBINE_COMPARISON.out.comparison_json
}
