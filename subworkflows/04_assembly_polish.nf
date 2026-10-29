// ============================================================
// 04_assembly_polish.nf
// Medaka polishing of the contigs, before the BlobToolKit split, and the gate
// that decides whether the polished contigs go on.
//
// The gate compares the raw Medaka output with its own input (same contigs, same
// names), on Merqury QV and BUSCO; on a regression the unpolished contigs go on.
// No Racon pre-step (see modules/04_assembly_polish.nf for the reasoning).
// ============================================================

nextflow.enable.dsl = 2

include { MEDAKA_POLISH ; POLISH_QC_GATE } from '../modules/04_assembly_polish.nf'

// Aliased so the same process can measure both ends of the stage without its
// output filenames colliding -- each invocation gets its own `stage` label.
include { BUSCO_FINAL as BUSCO_PREPOLISH  } from '../modules/04_assembly_polish.nf'
include { BUSCO_FINAL as BUSCO_POSTPOLISH } from '../modules/04_assembly_polish.nf'
include { MERQURY_QV  as MERQURY_PREPOLISH  } from '../modules/04_assembly_polish.nf'
include { MERQURY_QV  as MERQURY_POSTPOLISH } from '../modules/04_assembly_polish.nf'

workflow ASSEMBLY_POLISHING {

    take:
        assembly        // contigs (mito contigs removed), value channel
        reads           // ONT polishing reads
        busco_lineage
        busco_db_dir
        meryl_db        // read k-mer DB from MERYL_READ_DB; [] when run_polish_qc=false

    main:
        def versions_ch = channel.empty()

        MEDAKA_POLISH(assembly, reads)
        def polished = MEDAKA_POLISH.out.medaka_assembly

        // BUSCO of the polished contigs: the gate's after-measurement, and the BlobToolKit
        // BUSCO track when the polished contigs go on (same sequence names as the BlobDir).
        BUSCO_POSTPOLISH(polished, busco_lineage, busco_db_dir, 'postpolish')
        versions_ch = versions_ch.mix(MEDAKA_POLISH.out.versions, BUSCO_POSTPOLISH.out.versions)

        def selected       = polished
        def selected_table = BUSCO_POSTPOLISH.out.busco_table
        def stage_busco    = BUSCO_POSTPOLISH.out.busco_summary
        def decision_file  = channel.empty()
        def decision_ch    = channel.value('polished')
        def qc_report      = channel.empty()
        def gate_mqc       = channel.empty()
        def merqury_pre_qv  = channel.empty()
        def merqury_post_qv = channel.empty()

        if (params.run_polish_qc.toString().toBoolean()) {
            BUSCO_PREPOLISH(assembly, busco_lineage, busco_db_dir, 'prepolish')
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
            qc_report     = POLISH_QC_GATE.out.report
            gate_mqc      = POLISH_QC_GATE.out.mqc
            decision_file = POLISH_QC_GATE.out.decision
            decision_ch   = POLISH_QC_GATE.out.decision.map { f -> f.text.trim() }
            decision_ch.subscribe { d ->
                if (d == 'unpolished') {
                    log.warn "Polish QC gate: Medaka degraded the contigs (04_polishing/qc/polish_qc_report.tsv); continuing with the UNPOLISHED contigs"
                }
            }

            // Gate decision -> which contigs, and which BUSCO table, go on
            selected = polished.combine(decision_ch).filter { _a, d -> d == 'polished' }.map { a, _d -> a }
                .mix(assembly.combine(decision_ch).filter { _a, d -> d == 'unpolished' }.map { a, _d -> a })
            selected_table = BUSCO_POSTPOLISH.out.busco_table.combine(decision_ch).filter { _t, d -> d == 'polished' }.map { t, _d -> t }
                .mix(BUSCO_PREPOLISH.out.busco_table.combine(decision_ch).filter { _t, d -> d == 'unpolished' }.map { t, _d -> t })
            stage_busco = BUSCO_PREPOLISH.out.busco_summary.mix(BUSCO_POSTPOLISH.out.busco_summary)
            versions_ch = versions_ch.mix(
                BUSCO_PREPOLISH.out.versions, MERQURY_PREPOLISH.out.versions,
                MERQURY_POSTPOLISH.out.versions, POLISH_QC_GATE.out.versions)
        }
        else {
            log.warn "run_polish_qc=false: Medaka polishing is NOT checked for regression"
        }

    emit:
        contigs            = selected.first()          // polished or unpolished contigs
        busco_table        = selected_table.first()    // BUSCO full_table.tsv of those contigs
        busco_stage_summaries = stage_busco            // prepolish, postpolish
        decision           = decision_ch               // 'polished' | 'unpolished'
        polish_decision    = decision_file             // polish_decision.txt (with run_polish_qc)
        polish_gate_mqc    = gate_mqc
        merqury_pre        = merqury_pre_qv
        merqury_post       = merqury_post_qv
        polish_qc_report   = qc_report
        versions           = versions_ch
}
