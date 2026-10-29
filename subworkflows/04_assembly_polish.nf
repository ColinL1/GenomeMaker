// ============================================================
// 04_assembly_polish.nf
// Assembly polishing: Medaka → funannotate clean/sort → QC + gate
//
// No Racon pre-step (see modules/04_assembly_polish.nf for the reasoning).
//
// Every metric is measured before AND after polishing; on a regression
// POLISH_QC_GATE selects the unpolished assembly.
//
// QC sits after funannotate clean/sort, not after medaka, because `clean` drops
// contigs <200 bp and `sort` renames scaffolds -- so the post-polish numbers must
// describe the assembly that actually leaves this stage.
// ============================================================

nextflow.enable.dsl = 2

include { MEDAKA_POLISH ; FUNANNOTATE_CLEAN ; FUNANNOTATE_SORT ; POLISH_QC_GATE } from '../modules/04_assembly_polish.nf'
// Same clean/sort for the unpolished assembly, run only when the QC gate falls back to it.
include { FUNANNOTATE_CLEAN as FUNANNOTATE_CLEAN_UNPOLISHED } from '../modules/04_assembly_polish.nf'
include { FUNANNOTATE_SORT  as FUNANNOTATE_SORT_UNPOLISHED  } from '../modules/04_assembly_polish.nf'

// Aliased so the same process can measure both ends of the stage without its
// output filenames colliding -- each invocation gets its own `stage` label.
include { BUSCO_FINAL as BUSCO_PREPOLISH  } from '../modules/04_assembly_polish.nf'
include { BUSCO_FINAL as BUSCO_POSTPOLISH } from '../modules/04_assembly_polish.nf'
// BUSCO of the final product when the QC gate falls back to the unpolished assembly.
include { BUSCO_FINAL as BUSCO_FINAL_UNPOLISHED } from '../modules/04_assembly_polish.nf'
include { MERQURY_QV  as MERQURY_PREPOLISH  } from '../modules/04_assembly_polish.nf'
include { MERQURY_QV  as MERQURY_POSTPOLISH } from '../modules/04_assembly_polish.nf'

// Reuse the existing gfastats process rather than defining a second one.
include { GFASTATS_QC as GFASTATS_PREPOLISH  } from '../modules/02_multi_assembler.nf'
include { GFASTATS_QC as GFASTATS_POSTPOLISH } from '../modules/02_multi_assembler.nf'

workflow ASSEMBLY_POLISHING {

    take:
        assembly        // filtered, pre-polish
        reads           // ONT polishing reads
        busco_lineage
        busco_db_dir
        meryl_db        // read k-mer DB from MERYL_READ_DB; [] disables Merqury QV

    main:
        def versions_ch = channel.empty()

        MEDAKA_POLISH(assembly, reads)
        FUNANNOTATE_CLEAN(MEDAKA_POLISH.out.medaka_assembly)
        FUNANNOTATE_SORT(FUNANNOTATE_CLEAN.out.clean_assembly, params.scaffold_basename, 'polished')

        def polished = FUNANNOTATE_SORT.out.final_assembly
        // Assembly that continues downstream: the polished one, unless the QC gate
        // below finds a regression and falls back to the unpolished one.
        def selected = polished
        // the same assembly before FUNANNOTATE_SORT renamed it (for the scaffold name map)
        def selected_presort = FUNANNOTATE_CLEAN.out.clean_assembly

        // BUSCO at both ends. The post-polish table is what downstream summary and
        // MultiQC consume, so it is emitted regardless of the QC toggle.
        BUSCO_POSTPOLISH(polished, busco_lineage, busco_db_dir, 'postpolish')
        versions_ch = versions_ch.mix(
            MEDAKA_POLISH.out.versions, FUNANNOTATE_CLEAN.out.versions, FUNANNOTATE_SORT.out.versions,
            BUSCO_POSTPOLISH.out.versions)

        def qc_report = channel.empty()
        def pre_busco = channel.empty()
        def merqury_pre_qv = channel.empty()
        def merqury_post_qv = channel.empty()
        // Final-product BUSCO: postpolish unless the gate selects the unpolished assembly.
        def final_busco_table   = BUSCO_POSTPOLISH.out.busco_table
        def final_busco_summary = BUSCO_POSTPOLISH.out.busco_summary
        def final_busco_json    = BUSCO_POSTPOLISH.out.busco_summary_json
        // Every polishing-stage BUSCO summary (for the report); each file appears once.
        def stage_busco_summaries = BUSCO_POSTPOLISH.out.busco_summary
        def decision_file = channel.empty()
        def final_qv_ch   = channel.empty()   // Merqury QV of the selected assembly (README / release summary)
        def gate_mqc = channel.empty()

        if (params.run_polish_qc.toString().toBoolean()) {
            BUSCO_PREPOLISH(assembly, busco_lineage, busco_db_dir, 'prepolish')
            pre_busco = BUSCO_PREPOLISH.out.busco_summary

            GFASTATS_PREPOLISH(assembly.map  { p -> ['prepolish',  p] })
            GFASTATS_POSTPOLISH(polished.map { p -> ['postpolish', p] })

            MERQURY_PREPOLISH(assembly,  meryl_db, 'prepolish')
            MERQURY_POSTPOLISH(polished, meryl_db, 'postpolish')
            merqury_pre_qv  = MERQURY_PREPOLISH.out.qv
            merqury_post_qv = MERQURY_POSTPOLISH.out.qv

            POLISH_QC_GATE(
                MERQURY_PREPOLISH.out.qv,
                MERQURY_POSTPOLISH.out.qv,
                BUSCO_PREPOLISH.out.busco_summary,
                BUSCO_POSTPOLISH.out.busco_summary
            )
            qc_report = POLISH_QC_GATE.out.report

            // Gate decision -> pick the assembly. The unpolished branch gets the same
            // funannotate clean/sort, and those tasks only run if the gate asks for it.
            def decision = POLISH_QC_GATE.out.decision.map { f -> f.text.trim() }
            def unpolished_in = assembly.combine(decision)
                .filter { _asm, d -> d == 'unpolished' }
                .map { asm, _d -> asm }
            FUNANNOTATE_CLEAN_UNPOLISHED(unpolished_in)
            FUNANNOTATE_SORT_UNPOLISHED(FUNANNOTATE_CLEAN_UNPOLISHED.out.clean_assembly,
                                        params.scaffold_basename, 'unpolished')
            selected = polished.combine(decision)
                .filter { _asm, d -> d == 'polished' }
                .map { asm, _d -> asm }
                .mix(FUNANNOTATE_SORT_UNPOLISHED.out.final_assembly)
            selected_presort = FUNANNOTATE_CLEAN.out.clean_assembly.combine(decision)
                .filter { _asm, d -> d == 'polished' }
                .map { asm, _d -> asm }
                .mix(FUNANNOTATE_CLEAN_UNPOLISHED.out.clean_assembly)

            BUSCO_FINAL_UNPOLISHED(FUNANNOTATE_SORT_UNPOLISHED.out.final_assembly,
                                   busco_lineage, busco_db_dir, 'final_unpolished')
            final_busco_table = BUSCO_POSTPOLISH.out.busco_table.combine(decision)
                .filter { _f, d -> d == 'polished' }.map { f, _d -> f }
                .mix(BUSCO_FINAL_UNPOLISHED.out.busco_table)
            final_busco_summary = BUSCO_POSTPOLISH.out.busco_summary.combine(decision)
                .filter { _f, d -> d == 'polished' }.map { f, _d -> f }
                .mix(BUSCO_FINAL_UNPOLISHED.out.busco_summary)
            final_busco_json = BUSCO_POSTPOLISH.out.busco_summary_json.combine(decision)
                .filter { _f, d -> d == 'polished' }.map { f, _d -> f }
                .mix(BUSCO_FINAL_UNPOLISHED.out.busco_summary_json)
            stage_busco_summaries = BUSCO_PREPOLISH.out.busco_summary
                .mix(BUSCO_POSTPOLISH.out.busco_summary)
                .mix(BUSCO_FINAL_UNPOLISHED.out.busco_summary)
            decision_file = POLISH_QC_GATE.out.decision
            final_qv_ch = MERQURY_POSTPOLISH.out.qv.combine(decision)
                .filter { _stage, _qv, d -> d == 'polished' }.map { _stage, qv, _d -> qv }
                .mix(MERQURY_PREPOLISH.out.qv.combine(decision)
                    .filter { _stage, _qv, d -> d == 'unpolished' }.map { _stage, qv, _d -> qv })
            gate_mqc      = POLISH_QC_GATE.out.mqc
            versions_ch = versions_ch.mix(
                BUSCO_PREPOLISH.out.versions, GFASTATS_PREPOLISH.out.versions,
                GFASTATS_POSTPOLISH.out.versions, MERQURY_PREPOLISH.out.versions,
                MERQURY_POSTPOLISH.out.versions, POLISH_QC_GATE.out.versions,
                FUNANNOTATE_CLEAN_UNPOLISHED.out.versions, FUNANNOTATE_SORT_UNPOLISHED.out.versions,
                BUSCO_FINAL_UNPOLISHED.out.versions)
        }
        else {
            log.warn "run_polish_qc=false: polishing will NOT be checked for regression"
        }

    emit:
        final_assembly     = selected
        presort_assembly   = selected_presort   // final_assembly before renaming (same sequences)
        // busco_* describe the FINAL product (the assembly in final_assembly)
        busco_tables       = final_busco_table
        busco_summaries    = final_busco_summary
        busco_summary_json = final_busco_json
        busco_prepolish    = pre_busco
        busco_stage_summaries = stage_busco_summaries  // prepolish, postpolish, final_unpolished
        polish_decision    = decision_file              // "polished" | "unpolished"
        final_qv           = final_qv_ch
        polish_gate_mqc    = gate_mqc
        merqury_pre        = merqury_pre_qv
        merqury_post       = merqury_post_qv
        polish_qc_report   = qc_report
        versions           = versions_ch
}
