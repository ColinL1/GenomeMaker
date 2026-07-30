// ============================================================
// 02b_select_best_assembly.nf (subworkflow)
// Picks the winning assembler's assembly path based on
// data-driven ranking rather than alphabetical order.
// ============================================================

nextflow.enable.dsl = 2

include { RANK_ASSEMBLIES } from '../modules/02b_select_best_assembly.nf'

workflow SELECT_BEST_ASSEMBLY {

    take:
        assemblies_ch          // (name, path) -- from MULTI_ASSEMBLER.out.assemblies
        gfastats_ch            // (name, path) -- from MULTI_ASSEMBLER.out.gfastats
        busco_host_summary_ch  // (name, path) -- from MULTI_ASSEMBLER.out.busco_host_summary

    main:
        RANK_ASSEMBLIES(
            gfastats_ch.map { name, p -> p }.collect(),
            busco_host_summary_ch.map { name, p -> p }.collect()
        )

        def best_name_ch = RANK_ASSEMBLIES.out.best_name
            .map { f -> f.text.trim() }
            .first()  // value channel: single winning name string

        def selected_assembly_ch = assemblies_ch
            .combine(best_name_ch)
            .filter { name, path, best -> name == best }
            .map    { name, path, best -> path }
            .first()

    emit:
        selected_assembly   = selected_assembly_ch           // value channel, single path
        selected_name       = best_name_ch
        ranking_table       = RANK_ASSEMBLIES.out.ranking_table
}
