// ============================================================
// 06_repeats.nf
// Repeat identification, masking, and annotation
// barrnap rRNA → RepeatModeler → EDTA → combine → RepeatMasker
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: BARRNAP_RRNA
// Predict rRNA genes using barrnap
// ============================================================
process BARRNAP_RRNA {
    label 'repeats_env'
    tag   "barrnap_rrna"
    cpus  { 8 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::barrnap'
    container 'quay.io/biocontainers/barrnap:1.10.6--pl5321hdfd78af_0'
    publishDir "${params.outdir}/06_repeats/barrnap", mode: params.publish_mode

    input:
    path assembly

    output:
    path "rrna.fasta", emit: rrna_fasta
    path "rrna.gff",  emit: rrna_gff

    script:
    """
    barrnap -q -k euk ${assembly} --threads ${task.cpus} \\
      --outseq rrna.fasta --gff rrna.gff
    """

    stub:
    """
    touch rrna.fasta rrna.gff
    echo "BARRNAP_RRNA stub completed"
    """
}

// ============================================================
// Process: REPEATMODELER
// Build de novo repeat library using RepeatModeler
// ============================================================
process REPEATMODELER {
    label 'repeats_env'
    tag   "repeatmodeler"
    cpus  { 20 * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 72.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::repeatmodeler'
    container 'quay.io/biocontainers/repeatmodeler:2.0.8--pl5321hdfd78af_1'
    publishDir "${params.outdir}/06_repeats/repeatmodeler", mode: params.publish_mode

    input:
    path assembly

    output:
    path "repeat_modeler_output", emit: repeat_lib

    script:
    """
    mkdir -p repeat_modeler_output

    # Build repeat database (name matches assembly basename)
    _dbname=\$(basename \${assembly} .fasta)
    BuildDatabase -name \${_dbname} \${assembly}

    # Run RepeatModeler
    RepeatModeler -database \${_dbname} -LTRStruct -pa \${task.cpus}

    # Collect outputs
    cp \${_dbname}.families repeat_modeler_output/repeat_lib.fa
    cp \${_dbname}.consensi.fa.classified repeat_modeler_output/consensi.fa.classified
    cp \${_dbname}.log repeat_modeler_output/repeatmodeler.log
    """

    stub:
    """
    mkdir -p repeat_modeler_output
    touch repeat_modeler_output/repeat_lib.fa
    echo "REPEATMODELER stub completed"
    """
}

// ============================================================
// Process: EDTA_RUN
// Run EDTA for transposable element annotation
// ============================================================
process EDTA_RUN {
    label 'repeats_env'
    tag   "edta"
    cpus  { 20 * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::edta'
    container 'quay.io/biocontainers/edta:2.3.0--hdfd78af_0'
    publishDir "${params.outdir}/06_repeats/edta", mode: params.publish_mode

    input:
    path assembly

    output:
    path "edta_output", emit: edta_lib

    script:
    """
    mkdir -p edta_output

    _base=\$(basename \${assembly} .fasta)

    # Run EDTA with sensitive mode
    EDTA.pl --genome \${assembly} \\
      --sensitive 1 \\
      --anno 1 \\
      -t \${task.cpus} \\
      --overwrite 1 \\
      --force 1

    # Collect TE library
    if [ -f \${_base}.EDTA.TElib.fa ]; then
      cp \${_base}.EDTA.TElib.fa edta_output/edta_lib.fa
    else
      # Fallback: search for any EDTA output
      find . -name "*.EDTA.TElib.fa" -exec cp {} edta_output/edta_lib.fa \\;
    fi
    """

    stub:
    """
    mkdir -p edta_output
    touch edta_output/edta_lib.fa
    echo "EDTA_RUN stub completed"
    """
}

// ============================================================
// Process: COMBINE_REPEAT_DB
// Combine RepeatModeler + EDTA libraries, deduplicate
// ============================================================
process COMBINE_REPEAT_DB {
    label 'repeats_env'
    tag   "combine_repeat_db"
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::vsearch'
    container 'quay.io/biocontainers/vsearch:2.31.0--hd2be7a0_0'
    publishDir "${params.outdir}/06_repeats/combined_library", mode: params.publish_mode

    input:
    path repeatmodeler_lib
    path edta_lib

    output:
    path "combined_rep.lib", emit: combined_lib
    path "combined_rep_dedup.lib", emit: dedup_lib

    script:
    """
    # Combine all repeat libraries
    cat ${repeatmodeler_lib} ${edta_lib} > combined_rep.lib

    # Deduplicate using vsearch
    vsearch --derep_fulllength combined_rep.lib \\
      --output combined_rep_dedup.lib \\
      --strand both \\
      --id 0.95

    # Generate library distribution
    grep '^>' combined_rep_dedup.lib | \\
      sed -E 's/>.*//' | \\
      sort | uniq -c | sort -rn > library_distribution.txt
    """

    stub:
    """
    touch combined_rep.lib combined_rep_dedup.lib
    echo "COMBINE_REPEAT_DB stub completed"
    """
}

// ============================================================
// Process: REPEATMASKER_RUN
// Mask repeats using RepeatMasker
// ============================================================
process REPEATMASKER_RUN {
    label 'repeats_env'
    tag   "repeatmasker"
    cpus  { 8 * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::repeatmasker'
    container 'quay.io/biocontainers/repeatmasker:4.2.4--pl5321hdfd78af_0'
    publishDir "${params.outdir}/06_repeats/repeatmasker", mode: params.publish_mode

    input:
    path assembly
    path repeat_lib

    output:
    path "masked_assembly.fasta", emit: masked_assembly
    path "repeatmasker.out",      emit: repeat_tbl

    script:
    """
    # Run RepeatMasker with the combined deduplicated library
    RepeatMasker ${assembly} \\
      -lib ${repeat_lib} \\
      -pa ${task.cpus} \\
      -norna \\
      -xsmall

    # The output file follows the pattern: <basename>.out.masked
    _base=\$(basename \${assembly} .fasta)
    mv \${_base}.out.masked masked_assembly.fasta

    # Move the table output
    mv \${_base}.out repeatmasker.out
    """

    stub:
    """
    touch masked_assembly.fasta repeatmasker.out
    echo "REPEATMASKER_RUN stub completed"
    """
}

// // ============================================================
// // Workflow: repeats_pipeline
// // rRNA prediction → RepeatModeler + EDTA → combine → RepeatMasker
// // ============================================================
// workflow repeats_pipeline {

//     take:
//         assembly

//     main:
//         BARRNAP_RRNA(assembly)
//         REPEATMODELER(assembly)
//         EDTA_RUN(assembly)
//         COMBINE_REPEAT_DB(REPEATMODELER.out.repeat_lib, EDTA_RUN.out.edta_lib)
//         REPEATMASKER_RUN(assembly, COMBINE_REPEAT_DB.out.dedup_lib)

//     emit:
//         masked_assembly = REPEATMASKER_RUN.out.masked_assembly
//         repeat_tbl      = REPEATMASKER_RUN.out.repeat_tbl
// }
