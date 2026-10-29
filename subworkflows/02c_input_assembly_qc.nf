// ============================================================
// 02c_input_assembly_qc.nf (subworkflow)
// Used when params.input_assembly is set: the pipeline is given a
// finished fasta assembly instead of building one from reads, so
// KMER_GENOME_SIZE / MULTI_ASSEMBLER / SELECT_BEST_ASSEMBLY are all
// skipped. This subworkflow runs the same gfastats/BUSCO QC that
// MULTI_ASSEMBLER would have produced, tagged with a fixed name, so
// downstream steps (filtering, polishing, summary, MultiQC) see the
// same channel shapes regardless of which path was taken.
// ============================================================

nextflow.enable.dsl = 2

include { BUSCO_ASSEMBLY; GFASTATS_QC; COMBINE_COMPARISON } from '../modules/02_multi_assembler.nf'

workflow INPUT_ASSEMBLY_QC {

    take:
        assembly       // path -- user-supplied fasta
        busco_lineage  // single lineage, e.g. anthozoa_odb12
        busco_db_dir
        combine_comparison_script

    main:
        def name = 'input_assembly'
        def assembly_ch = assembly.map { p -> [name, p] }

        BUSCO_ASSEMBLY(assembly_ch, channel.value(busco_lineage), busco_db_dir)
        GFASTATS_QC(assembly_ch)

        COMBINE_COMPARISON(
            GFASTATS_QC.out.gfastats_report.map { n, p -> p }.collect(),
            BUSCO_ASSEMBLY.out.busco_summary.map { n, p -> p }.collect(),
            combine_comparison_script
        )

    emit:
        assemblies       = assembly_ch                          // (name, path)
        selected_assembly = assembly                             // value channel, single path
        selected_name     = channel.value(name)
        gfastats          = GFASTATS_QC.out.gfastats_report       // (name, path)
        busco_tables      = BUSCO_ASSEMBLY.out.busco_table        // (name, path)
        busco_summary     = BUSCO_ASSEMBLY.out.busco_summary      // (name, path)
        comparison        = COMBINE_COMPARISON.out.comparison_tsv
        comparison_json   = COMBINE_COMPARISON.out.comparison_json
        versions        = BUSCO_ASSEMBLY.out.versions.mix(GFASTATS_QC.out.versions, COMBINE_COMPARISON.out.versions)
}
