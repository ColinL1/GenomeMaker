// ============================================================
// 02b_select_best_assembly.nf (subworkflow)
// Picks the winning assembly by ranking (bin/rank.py).
// ============================================================

nextflow.enable.dsl = 2

include { RANK_ASSEMBLIES } from '../modules/02b_select_best_assembly.nf'

workflow SELECT_BEST_ASSEMBLY {

    take:
        assemblies_ch          
        gfastats_ch            
        busco_host_summary_ch  
        quast_ch               
        expected_genome_size   // GenomeScope2 or params.genome_size

    main:
        RANK_ASSEMBLIES(
            gfastats_ch.map { name, p -> p }.collect(),
            busco_host_summary_ch.map { name, p -> p }.collect(),
            quast_ch.map { name, p -> p }.collect(),
            expected_genome_size
        )

        def best_name_ch = RANK_ASSEMBLIES.out.best_name
            .map { f -> f.text.trim() }
            .first()

        def selected_assembly_ch = assemblies_ch
            .combine(best_name_ch)
            .filter { name, path, best -> name == best }
            .map    { name, path, best -> path }
            .first()

    emit:
        selected_assembly   = selected_assembly_ch
        selected_name       = best_name_ch
        ranking_table       = RANK_ASSEMBLIES.out.ranking_table
        versions            = RANK_ASSEMBLIES.out.versions
}
