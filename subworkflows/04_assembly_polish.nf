// ============================================================
// 04_assembly_polish.nf
// Assembly polishing: Racon (ONT) → Medaka (ONT) → funannotate clean/sort
// Final BUSCO on polished assembly
// ============================================================

nextflow.enable.dsl = 2

include { RACON_POLISH ; MEDAKA_POLISH ; FUNANNOTATE_CLEAN_SORT ; BUSCO_FINAL ; BUSCO_FINAL_META } from '../modules/04_assembly_polish.nf'

// ============================================================
// Workflow: assembly_polishing
// Racon → Medaka → clean/sort → final BUSCO
// ============================================================
workflow ASSEMBLY_POLISHING {

    take:
        assembly
        reads
        busco_lineage_euk
        // busco_lineage_meta

    main:
        RACON_POLISH(assembly, reads)
        MEDAKA_POLISH(RACON_POLISH.out.polished_assembly, reads, channel.value([]))
        FUNANNOTATE_CLEAN_SORT(MEDAKA_POLISH.out.medaka_assembly, params.species_name)
        BUSCO_FINAL(FUNANNOTATE_CLEAN_SORT.out.final_assembly, busco_lineage_euk)
        // BUSCO_FINAL_META(FUNANNOTATE_CLEAN_SORT.out.final_assembly, busco_lineage_meta)

    emit:
        final_assembly  = FUNANNOTATE_CLEAN_SORT.out.final_assembly
        busco_tables    = BUSCO_FINAL.out.busco_table
        busco_summaries = BUSCO_FINAL.out.busco_summary
}
