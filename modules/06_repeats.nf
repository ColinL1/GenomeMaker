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

    // barrnap 0.9: 1.10.x has no eukaryotic database
    conda 'bioconda::barrnap=0.9'
    container 'https://depot.galaxyproject.org/singularity/barrnap:0.9--hdfd78af_4'
    publishDir "${params.outdir}/06_repeats/barrnap", mode: params.publish_mode

    input:
    path assembly

    output:
    path "rrna.fasta", emit: rrna_fasta
    path "rrna.gff",  emit: rrna_gff
    path "versions.yml", emit: versions

    script:
    """
    # barrnap has no --gff option: the GFF is written to stdout.
    barrnap -q -k euk ${assembly} --threads ${task.cpus} \\
      --outseq rrna.fasta > rrna.gff
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barrnap: "\$(barrnap --version 2>&1 | sed 's/barrnap //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch rrna.fasta rrna.gff
    echo "BARRNAP_RRNA stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barrnap: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: REPEATMODELER
// Build de novo repeat library using RepeatModeler
// ============================================================
process REPEATMODELER {
    label 'repeats_env'
    tag   "repeatmodeler"

    conda 'bioconda::repeatmodeler'
    container 'quay.io/biocontainers/repeatmodeler:2.0.8--pl5321hdfd78af_1'
    publishDir "${params.outdir}/06_repeats/repeatmodeler", mode: params.publish_mode

    input:
    path assembly

    output:
    path "repeat_modeler_output", emit: repeat_lib
    path "versions.yml", emit: versions

    script:
    def db = assembly.baseName
    """
    mkdir -p repeat_modeler_output

    BuildDatabase -name ${db} ${assembly}

    RepeatModeler -database ${db} -LTRStruct -threads ${task.cpus}

    cp ${db}-families.fa  repeat_modeler_output/repeat_lib.fa
    cp ${db}-families.stk repeat_modeler_output/repeat_families.stk
    cp ${db}-rmod.log     repeat_modeler_output/repeatmodeler.log 2>/dev/null || true

    if [ ! -s repeat_modeler_output/repeat_lib.fa ]; then
        echo "ERROR: RepeatModeler produced no repeat families for ${assembly}" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: "\$(RepeatModeler --version 2>&1 | sed 's/.*version //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p repeat_modeler_output
    printf '>stub_repeat\\nACGTACGTACGTACGT\\n' > repeat_modeler_output/repeat_lib.fa
    touch repeat_modeler_output/repeat_families.stk
    echo "REPEATMODELER stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: EDTA_RUN
// Run EDTA for transposable element annotation
// ============================================================
process EDTA_RUN {
    label 'repeats_env'
    tag   "edta"

    // TIR-Learner 4: the biocontainer's TIR-Learner 3.0.2 crashes under pandas 3 (KeyError: 0), so no TIRs
    conda 'bioconda::edta=2.3.0 bioconda::tir-learner=4.08'
    container 'oras://community.wave.seqera.io/library/edta_tir-learner:f137297b45b5e446'
    publishDir "${params.outdir}/06_repeats/edta", mode: params.publish_mode

    input:
    path assembly

    output:
    path "edta_output", emit: edta_lib
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p edta_output

    EDTA.pl --genome ${assembly} \\
      --sensitive 1 \\
      --anno 1 \\
      -t ${task.cpus} \\
      --overwrite 1 \\
      --force 1

    # EDTA names the library after the input file, extension included
    if [ -f "${assembly}.EDTA.TElib.fa" ]; then
      cp "${assembly}.EDTA.TElib.fa" edta_output/edta_lib.fa
    else
      # EDTA moves outputs between release versions; fall back to a search.
      find . -name "*.EDTA.TElib.fa" -exec cp {} edta_output/edta_lib.fa \\;
    fi

    if [ ! -s edta_output/edta_lib.fa ]; then
        echo "ERROR: EDTA produced no TE library for ${assembly}" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edta: "\$(EDTA.pl -h 2>&1 | grep -im1 'v[0-9]' | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p edta_output
    touch edta_output/edta_lib.fa
    echo "EDTA_RUN stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edta: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: COMBINE_REPEAT_DB
// Combine RepeatModeler + EDTA libraries, deduplicate
// ============================================================
process COMBINE_REPEAT_DB {
    label 'repeats_env'
    tag   "combine_repeat_db"

    conda 'bioconda::vsearch'
    container 'quay.io/biocontainers/vsearch:2.31.0--hd2be7a0_0'
    publishDir "${params.outdir}/06_repeats/combined_library", mode: params.publish_mode

    input:
    path repeatmodeler_lib
    path edta_lib

    output:
    path "combined_rep.lib", emit: combined_lib
    path "combined_rep_dedup.lib", emit: dedup_lib
    path "versions.yml", emit: versions

    script:
    """
    # both inputs are directories; cat the FASTA inside each
    cat ${repeatmodeler_lib}/repeat_lib.fa ${edta_lib}/edta_lib.fa > combined_rep.lib

    # exact duplicates on either strand; vsearch keeps <id>#<class/family> labels
    vsearch --derep_fulllength combined_rep.lib \\
      --output combined_rep_dedup.lib \\
      --strand both

    # Families per repeat class (the part after '#' in each header)
    grep '^>' combined_rep_dedup.lib | \\
      sed -E 's/^[^#]*#//; s/[[:space:]].*//' | \\
      sort | uniq -c | sort -rn > library_distribution.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: "\$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//;s/,.*//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch combined_rep.lib combined_rep_dedup.lib
    echo "COMBINE_REPEAT_DB stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: REPEATMASKER_RUN
// Mask repeats using RepeatMasker
// ============================================================
process REPEATMASKER_RUN {
    label 'repeats_env'
    tag   "repeatmasker"

    // RepeatModeler 2.0.8 image: the repeatmasker:4.2.4 image has no FamDB and aborts even with -lib
    conda 'bioconda::repeatmodeler=2.0.8'
    container 'quay.io/biocontainers/repeatmodeler:2.0.8--pl5321hdfd78af_1'
    publishDir "${params.outdir}/06_repeats/repeatmasker", mode: params.publish_mode

    input:
    path assembly
    path repeat_lib

    output:
    path "masked_assembly.fasta", emit: masked_assembly
    path "repeatmasker.out",      emit: repeat_tbl
    path "repeatmasker.tbl",      emit: repeat_summary
    path "versions.yml", emit: versions

    script:
    def pa_jobs = Math.max(1, (task.cpus / 4) as int)
    """
    # -pa counts jobs of 4 cores each (RMBlast). -xsmall softmasks for BRAKER.
    RepeatMasker ${assembly} \\
      -lib ${repeat_lib} \\
      -pa ${pa_jobs} \\
      -norna \\
      -xsmall

    # -xsmall also lower-cases gap runs; NCBI and gap detection expect 'N'
    awk '/^>/{print; next}{gsub(/n/, "N"); print}' ${assembly}.masked > masked_assembly.fasta
    mv ${assembly}.out    repeatmasker.out

    mv ${assembly}.tbl    repeatmasker.tbl
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: "\$(RepeatMasker -v 2>&1 | sed 's/.*version //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch masked_assembly.fasta repeatmasker.out repeatmasker.tbl
    echo "REPEATMASKER_RUN stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: "stub"
    END_VERSIONS
    """
}
