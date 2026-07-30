// ============================================================
// 07_structural_annotation.nf
// Structural gene annotation pipeline:
// tRNAscan-SE → RNA-seq trimming → STAR mapping → BRAKER3 → AGAT merge
// ============================================================

nextflow.enable.dsl = 2

include { TRNASCAN_SE ; RNASEQ_TRIM ; STAR_GENOME_INDEX ; STAR_MAPPING ; MERGE_BAMS ; SPLIT_STRANDED_BAM ; BRAKER3_RUN ; MERGE_ANNOTATIONS } from '../modules/07_structural_annotation.nf'

// ============================================================
// Workflow: STRUCTURAL_ANNOTATION
// tRNAscan-SE → RNA-seq → STAR → BRAKER3 → AGAT merge
//
// RNA-seq handling: `rnaseq_dir` (a plain directory path, not a
// channel) is used to build a channel of paired-end fastq tuples
// via `channel.fromFilePairs`, one item per sample. RNASEQ_TRIM and
// STAR_MAPPING are then forked/parallelized per-sample by Nextflow
// itself (no in-process for-loop over libraries). The per-sample
// STAR BAMs are collected and merged by a single MERGE_BAMS task.
// ============================================================
workflow STRUCTURAL_ANNOTATION {

    take:
        assembly
        rnaseq_dir       // plain directory path containing *_R1.fastq.gz / *_R2.fastq.gz
        protein_refs

    main:
        def paired_reads_ch = channel.fromFilePairs("${rnaseq_dir}/*_R{1,2}.fastq.gz", checkIfExists: true)

        TRNASCAN_SE(assembly)
        RNASEQ_TRIM(paired_reads_ch)
        STAR_GENOME_INDEX(assembly)
        STAR_MAPPING(STAR_GENOME_INDEX.out.star_index_dir.first(), RNASEQ_TRIM.out.trimmed_reads)
        MERGE_BAMS(STAR_MAPPING.out.sample_bam.collect())
        SPLIT_STRANDED_BAM(MERGE_BAMS.out.rnaseq_bam)
        BRAKER3_RUN(assembly, params.rnaseq_strandedness, SPLIT_STRANDED_BAM.out.plus_bam, SPLIT_STRANDED_BAM.out.minus_bam, protein_refs)
        MERGE_ANNOTATIONS(BRAKER3_RUN.out.braker_gtf, TRNASCAN_SE.out.trnascan_gff, assembly)

    emit:
        gff3     = MERGE_ANNOTATIONS.out.gff3
        proteins = MERGE_ANNOTATIONS.out.proteins_faa
}
