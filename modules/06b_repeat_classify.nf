// ============================================================
// 06b_repeat_classify.nf
// REPEAT_CLASSIFY: RepeatClassifier on the 'Unknown' families of the combined
// library, BEFORE RepeatMasker, so RepeatMasker's .out/.tbl carry the final
// labels. Labels are normalised (unknown -> Unknown, DNAnona/DNAauto -> DNA, ...).
// REPEAT_GFF: repeat GFF3 from RepeatMasker's .out, with the repeat class.
// Scripts: bin/repeat_relabel.py (library-only mode), bin/repeat_gff_classes.py.
// ============================================================

nextflow.enable.dsl = 2

process REPEAT_CLASSIFY {
    tag "repeat_classify"

    conda 'bioconda::repeatmodeler=2.0.8'
    container 'quay.io/biocontainers/repeatmodeler:2.0.8--pl5321hdfd78af_1'
    // Optional Dfam partitions (RepeatMasker Libraries/ directory) bound over the image's own.
    containerOptions "${params.rm_libraries ? "-B ${params.rm_libraries}:/usr/local/share/RepeatMasker/Libraries" : ''}"
    publishDir "${params.outdir}/06_repeats/classify", mode: params.publish_mode

    input:
    path lib            // combined, deduplicated library
    val  rm_libraries   // part of the task hash: a library change reruns the classifier

    output:
    path "repeats_classified.lib", emit: lib
    path "classify_report.tsv",    emit: report
    path "classify.log"
    path "versions.yml", emit: versions

    script:
    """
    awk '/^>/{keep=(\$0 ~ /#[Uu]nknown( |\$)/); if(keep){split(\$1,a,"#"); print a[1]; next}} keep' ${lib} > unknown_families.fa
    echo "unknown families: \$(grep -c '>' unknown_families.fa)" > classify.log
    if [ -s unknown_families.fa ]; then
        RepeatClassifier -consensi unknown_families.fa -threads ${task.cpus} >> classify.log 2>&1
    else
        touch unknown_families.fa.classified
    fi
    repeat_relabel.py ${lib} unknown_families.fa.classified - repeats_classified.lib - classify_report.tsv >> classify.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: "\$(RepeatModeler --version 2>&1 | grep -o '[0-9.]*' | head -1 || echo unknown)"
        rm_libraries: "${rm_libraries ?: 'image default'}"
    END_VERSIONS
    """

    stub:
    """
    cp ${lib} repeats_classified.lib
    touch classify_report.tsv classify.log
    printf '"${task.process}":\\n    repeatmodeler: stub\\n' > versions.yml
    """
}

process REPEAT_GFF {
    label 'res_light'
    tag   "repeat_gff"

    conda 'bioconda::repeatmodeler=2.0.8'
    container 'quay.io/biocontainers/repeatmodeler:2.0.8--pl5321hdfd78af_1'
    publishDir "${params.outdir}/06_repeats/repeatmasker", mode: params.publish_mode

    input:
    path rm_out

    output:
    path "repeats.gff3", emit: gff

    script:
    """
    # the /usr/local/bin copies cannot find RepeatMasker's Perl modules; use the util/ ones
    perl /usr/local/share/RepeatMasker/util/rmOutToGFF3.pl ${rm_out} > repeats_noclass.gff3
    repeat_gff_classes.py ${rm_out} repeats_noclass.gff3 repeats.gff3
    """

    stub:
    """
    touch repeats.gff3
    """
}
