// ============================================================
// 05_mitogenome.nf
// Mitochondrial genome extraction by mapping assembly reads
// to a reference mitogenome
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: MITO_READ_FILTER
// Filter assembly reads for high quality (Q>=10) for mitogenome
// ============================================================
process MITO_READ_FILTER {
    label 'mito_env'
    label 'res_read_tool'
    tag   "mito_read_filter"

    conda 'bioconda::chopper'
    container 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'
    publishDir "${params.outdir}/05_mitogenome", mode: params.publish_mode

    input:
    path reads

    output:
    path "mito_reads_q10.fastq.gz", emit: mito_reads
    path "versions.yml", emit: versions

    script:
    """
    # Filter for quality >= 10 for mitogenome mapping
    ( set -o pipefail
      gunzip -c ${reads} \\
        | chopper -q 10 \\
        | gzip > mito_reads_q10.fastq.gz )
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "\$(chopper --version 2>&1 | sed 's/chopper //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mito_reads_q10.fastq.gz
    echo "MITO_READ_FILTER stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_MAPPING
// Map filtered reads to mitogenome reference, produce BAM
// ============================================================
process MITO_MAPPING {
    label 'mito_env'
    tag   "mito_mapping"

    conda 'bioconda::minimap2 samtools'
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
    publishDir "${params.outdir}/05_mitogenome", mode: params.publish_mode

    input:
    path reads
    path mito_ref

    output:
    path "mito_mapped.sorted.bam",  emit: mito_bam
    path "mito_mapped.sorted.bam.bai", emit: mito_bai
    path "versions.yml", emit: versions

    script:
    """
    # Create minimap2 index for reference
    minimap2 -d mito_ref.mmi ${mito_ref}

    # Map reads to mitogenome reference
    ( set -o pipefail
      minimap2 -ax map-ont mito_ref.mmi ${reads} \\
        | samtools sort -@ ${task.cpus} -O BAM -o mito_mapped.sorted.bam - )

    samtools index mito_mapped.sorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version 2>&1 | head -1 || echo unknown)"
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mito_mapped.sorted.bam mito_mapped.sorted.bam.bai
    echo "MITO_MAPPING stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "stub"
        samtools: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_COVERAGE_STATS
// Generate coverage statistics for mitogenome
// ============================================================
process MITO_COVERAGE_STATS {
    label 'mito_env'
    label 'res_tiny'
    tag   "mito_coverage"

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/05_mitogenome", mode: params.publish_mode

    input:
    path bam
    path bai

    output:
    path "mito_coverage.txt", emit: coverage
    path "mito_flagstat.txt", emit: flagstat
    path "versions.yml", emit: versions

    script:
    """
    samtools coverage ${bam} > mito_coverage.txt
    samtools flagstat ${bam} > mito_flagstat.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mito_coverage.txt mito_flagstat.txt
    echo "MITO_COVERAGE_STATS stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITOHIFI
// Recover the mitogenome from the assembly contigs, circularise,
// rotate and annotate it (MitoHiFi contig mode, MITOS2 annotation).
//
// Takes the PRE-purge assembly: purge_dups removes the mito contig
// as HIGHCOV, so it is absent from every later assembly.
// MitoHiFi requires the related mitogenome as FASTA and GenBank. MITOS2
// annotates de novo, so a feature-less GenBank written from the FASTA is
// enough and no annotated reference is needed (MitoFinder would need one).
// ============================================================
process MITOHIFI {
    label 'mito_env'
    label 'res_read_tool'
    tag   "mitohifi"

    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'
    publishDir "${params.outdir}/05_mitogenome/mitohifi", mode: params.publish_mode

    input:
    path assembly
    path mito_ref

    output:
    path "final_mitogenome.fasta",          emit: fasta,       optional: true
    path "final_mitogenome.gff",            emit: gff,         optional: true
    path "final_mitogenome.annotation.png", emit: plot,        optional: true
    path "contigs_stats.tsv",               emit: stats,       optional: true
    path "shared_genes.tsv",                emit: shared,      optional: true
    path "potential_contigs",               emit: contigs_dir, optional: true
    path "contigs_ids.txt",                 emit: contig_ids,  optional: true   // contigs matching the reference
    path "mitohifi.log",                    emit: log
    path "versions.yml",                    emit: versions

    script:
    """
    python3 -c "
from Bio import SeqIO
r = SeqIO.read('${mito_ref}', 'fasta')
r.annotations['molecule_type'] = 'DNA'
r.annotations['topology'] = 'circular'
SeqIO.write(r, 'mito_ref.gb', 'genbank')
"

    status=0
    mitohifi.py -c ${assembly} -f ${mito_ref} -g mito_ref.gb \\
        -t ${task.cpus} -o ${params.mito_genetic_code} --mitos > mitohifi.log 2>&1 || status=\$?

    # No contig matching the reference is a result, not a failure: the stage
    # then publishes only the log. Any other error stops the run.
    if [ \$status -ne 0 ]; then
        cat mitohifi.log >&2
        grep -q "'parsed_blast.txt' and 'parsed_blast_all.txt' files are empty" mitohifi.log || exit \$status
        echo "MITOHIFI: no mitochondrial contig found in ${assembly}" >&2
    fi
    if [ -s contigs_filtering/contigs_ids.txt ]; then cp contigs_filtering/contigs_ids.txt contigs_ids.txt; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: "\$(mitohifi.py -v 2>&1 | tail -1 | sed 's/MitoHiFi //')"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch final_mitogenome.fasta final_mitogenome.gff mitohifi.log
    printf 'contig_id\tframeshifts_found\tannotation_file\tlength(bp)\tnumber_of_genes\twas_circular\n' > contigs_stats.tsv
    echo "MITOHIFI stub completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Mitogenome, final product
//
// Two routes to the sequence, chosen by whether Illumina shotgun DNA of the
// same individual is supplied (params.illumina_dna_r1/_r2):
//   supplied -> MITO_GETORGANELLE (de novo from the Illumina reads)
//   absent   -> MITOHIFI contig from the ONT assembly, polished by MITO_MEDAKA
//               with the mito-mapped ONT reads of mean quality >= mito_polish_min_q
// Both then go through the same steps: MITO_ROTATE, MITO_LIFTOFF (curated
// reference annotation), MITO_MITOS2 (de novo annotation), MITO_GENE_CHECK.
// ============================================================

// ============================================================
// Process: GETORGANELLE_DB
// One-time download of GetOrganelle's animal_mt seed/label database.
// ============================================================
process GETORGANELLE_DB {
    label 'mito_env'
    label 'res_small'
    tag   "getorganelle_db"

    container 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'
    storeDir "${params.getorganelle_db_dir}"

    output:
    path "animal_mt_db", emit: db

    script:
    """
    # a preloaded host library makes ld.so print an ERROR line that GetOrganelle
    # reads as "Blast is not available"
    unset LD_PRELOAD
    get_organelle_config.py -a animal_mt --config-dir animal_mt_db
    """

    stub:
    """
    mkdir -p animal_mt_db
    """
}

// ============================================================
// Process: MITO_GETORGANELLE
// Assemble the mitogenome from Illumina shotgun DNA. Several runs of the
// same library set are concatenated in the order given (R1 and R2 lists
// must be in the same order); the order is written to the log.
// A result that is not a complete circular genome publishes only the log.
// ============================================================
process MITO_GETORGANELLE {
    label 'mito_env'
    tag   "getorganelle"

    container 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'
    publishDir "${params.outdir}/05_mitogenome/getorganelle", mode: params.publish_mode, overwrite: true

    input:
    path r1, stageAs: 'R1_??.fastq.gz'
    path r2, stageAs: 'R2_??.fastq.gz'
    val  read_names
    path seed
    path db

    output:
    path "getorganelle_mitogenome.fasta", emit: fasta, optional: true
    path "getorganelle_mitogenome.gfa",   emit: gfa,   optional: true
    path "get_org.log.txt",               emit: log
    path "versions.yml",                  emit: versions

    script:
    """
    unset LD_PRELOAD
    cat R1_*.fastq.gz > illumina_R1.fastq.gz
    cat R2_*.fastq.gz > illumina_R2.fastq.gz

    get_organelle_from_reads.py -1 illumina_R1.fastq.gz -2 illumina_R2.fastq.gz \\
        -F animal_mt -s ${seed} -o getorganelle -R 15 -k 21,45,65,85,105 \\
        -t ${task.cpus} --config-dir ${db}

    { echo "Input read files, concatenated in this order (R1 then R2 lists): ${read_names}"
      cat getorganelle/get_org.log.txt; } > get_org.log.txt

    # keep only a complete circular genome; path 1 if GetOrganelle reports isomers
    if grep -q "Result status of animal_mt: circular genome" getorganelle/get_org.log.txt; then
        cp \$(ls getorganelle/*complete*.1.path_sequence.fasta | head -1) getorganelle_mitogenome.fasta
        cp \$(ls getorganelle/*complete*.selected_graph.gfa | head -1) getorganelle_mitogenome.gfa
    else
        echo "MITO_GETORGANELLE: no complete circular mitogenome; see get_org.log.txt" >&2
    fi
    rm -f illumina_R1.fastq.gz illumina_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: "\$(get_organelle_from_reads.py --version 2>&1 | tail -1 | sed 's/GetOrganelle v//')"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch getorganelle_mitogenome.fasta get_org.log.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_MAPPED_READS
// ONT reads that MITO_MAPPING aligned to the mito reference (primary only).
// ============================================================
process MITO_MAPPED_READS {
    label 'mito_env'
    label 'res_tiny'
    tag   "mito_mapped_reads"

    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

    input:
    path bam
    path bai

    output:
    path "mito_mapped_reads.fastq.gz", emit: reads
    path "versions.yml", emit: versions

    script:
    """
    ( set -o pipefail
      samtools fastq -F 0x904 ${bam} | gzip > mito_mapped_reads.fastq.gz )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mito_mapped_reads.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_READS_HQ
// Keep mito reads of mean quality >= params.mito_polish_min_q. Tested on this
// dataset: Q19-Q21 reproduced the Illumina mitogenome exactly; lower cut-offs
// kept both homopolymer errors, higher ones had too little depth.
// ============================================================
process MITO_READS_HQ {
    label 'mito_env'
    label 'res_tiny'
    tag   "mito_reads_q${params.mito_polish_min_q}"

    container 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'

    input:
    path reads

    output:
    path "mito_reads_hq.fastq.gz", emit: reads
    path "versions.yml", emit: versions

    script:
    """
    ( set -o pipefail
      gunzip -c ${reads} | chopper -q ${params.mito_polish_min_q} | gzip > mito_reads_hq.fastq.gz )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "\$(chopper --version 2>&1 | sed 's/chopper //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mito_reads_hq.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_MEDAKA
// Polish the ONT mitogenome with the high-quality mito reads. The model is
// params.medaka_model, or read from the headers of the original ONT reads
// exactly as MEDAKA_POLISH does (the mito reads come from a BAM and have
// lost their header tags).
// ============================================================
process MITO_MEDAKA {
    label 'mito_env'
    label 'res_read_tool'
    tag   "mito_medaka"

    container 'quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0'
    publishDir "${params.outdir}/05_mitogenome/ont_polish", mode: params.publish_mode, overwrite: true

    input:
    path draft
    path hq_reads
    path model_reads

    output:
    path "medaka_mitogenome.fasta", emit: fasta
    path "mito_polish_stats.tsv",   emit: stats
    path "versions.yml",            emit: versions

    script:
    def explicit_model = params.medaka_model ?: ''
    """
    if [ -n "${explicit_model}" ]; then
        MODEL="${explicit_model}"
    else
        BCS=\$(gunzip -c ${model_reads} 2>/dev/null | awk 'NR%4==1' | head -n 20000 \\
              | tr ' \\t' '\\n\\n' \\
              | { grep -E '^(basecall_model_version_id=|RG:Z:)' || true; } \\
              | sed -E 's/^basecall_model_version_id=//; s/^RG:Z:[^_]*_//' \\
              | sort -u)
        if [ "\$(printf '%s\\n' "\$BCS" | grep -c . || true)" -ne 1 ]; then
            echo "ERROR: expected exactly one basecaller model in the headers of ${model_reads}, found:" >&2
            echo "\$BCS" >&2
            echo "       Set --medaka_model explicitly." >&2
            exit 1
        fi
        MODEL="\${BCS}:consensus"
    fi
    medaka tools resolve_model --model "\$MODEL"

    medaka_consensus -i ${hq_reads} -d ${draft} -o medaka_out -t ${task.cpus} -m "\$MODEL"
    cp medaka_out/consensus.fasta medaka_mitogenome.fasta

    gunzip -c ${hq_reads} | awk -v q=${params.mito_polish_min_q} -v m="\$MODEL" -v OFS='\\t' '
        NR % 4 == 2 { n++; b += length(\$0) }
        END { print "min_read_quality", q; print "reads_used", n + 0; print "bases_used", b + 0; print "medaka_model", m }' > mito_polish_stats.tsv
    awk '!/^>/ { l += length(\$0) } END { print "polished_length\\t" l }' medaka_mitogenome.fasta >> mito_polish_stats.tsv
    awk '!/^>/ { l += length(\$0) } END { print "draft_length\\t" l }' ${draft} >> mito_polish_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: "\$(medaka --version 2>&1 | sed 's/medaka //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch medaka_mitogenome.fasta mito_polish_stats.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_ROTATE
// Put the circular mitogenome on the strand and origin of the annotation
// reference. Needed for the liftover: a gene cut by the linearisation point
// cannot be transferred (tRNA-Phe was lost that way in testing).
// ============================================================
process MITO_ROTATE {
    label 'mito_env'
    label 'res_tiny'
    tag   "mito_rotate"

    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
    publishDir "${params.outdir}/05_mitogenome/final", mode: params.publish_mode, overwrite: true

    input:
    tuple val(source), path(mito, stageAs: 'input_mitogenome.fa')
    path ref, stageAs: 'annotation_ref.fa'

    output:
    path "mitogenome.fasta",       emit: fasta
    path "mitogenome_source.tsv",  emit: info
    path "versions.yml",           emit: versions

    script:
    """
    awk 'NR==1{print ">mitogenome"; next} /^>/{exit} {print}' ${mito} > in.fa
    samtools faidx ${ref}
    refname=\$(cut -f1 ${ref}.fai | head -1)
    samtools faidx ${ref} "\$refname:1-1000" > ref_start.fa

    # strand of the annotation reference
    strand=\$(minimap2 -c -x asm20 in.fa ref_start.fa 2>/dev/null | sort -k10,10nr | head -1 | cut -f5)
    if [ -z "\$strand" ]; then
        echo "ERROR: the start of ${ref} does not align to the mitogenome" >&2
        exit 1
    fi
    if [ "\$strand" = "-" ]; then
        samtools faidx in.fa
        samtools faidx -i in.fa mitogenome | sed '1s/.*/>mitogenome/' > oriented.fa
    else
        cp in.fa oriented.fa
    fi
    samtools faidx oriented.fa
    len=\$(cut -f2 oriented.fa.fai)

    # origin: the base aligned to reference position 1
    read qstart tstart < <(minimap2 -c -x asm20 oriented.fa ref_start.fa 2>/dev/null | sort -k3,3n -k10,10nr | head -1 | cut -f3,8)
    start=\$(( (tstart - qstart + len) % len + 1 ))

    if [ "\$start" -eq 1 ]; then
        cp oriented.fa mitogenome.fasta
    else
        { echo ">mitogenome"
          { samtools faidx oriented.fa "mitogenome:\$start-\$len"; samtools faidx oriented.fa "mitogenome:1-\$((start - 1))"; } \\
              | grep -v ">" | tr -d '\\n' | fold -w 60
          echo; } > mitogenome.fasta
    fi

    printf 'source\\t%s\\nlength\\t%s\\nannotation_reference\\t%s\\ninput_strand\\t%s\\nrotation_start\\t%s\\n' \\
        "${source}" "\$len" "\$refname" "\$strand" "\$start" > mitogenome_source.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version 2>&1 | head -1 || echo unknown)"
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mitogenome.fasta mitogenome_source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_LIFTOFF
// Transfer the curated annotation of a closely related mitogenome.
// Liftoff's own valid_ORFs flag assumes the standard genetic code and is
// meaningless here; MITO_GENE_CHECK does the translation check.
// ============================================================
process MITO_LIFTOFF {
    label 'mito_env'
    label 'res_tiny'
    tag   "mito_liftoff"

    container 'quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0'
    publishDir "${params.outdir}/05_mitogenome/final", mode: params.publish_mode, overwrite: true

    input:
    path mito
    path ref,     stageAs: 'annotation_ref.fa'
    path ref_gff, stageAs: 'annotation_ref.gff3'

    output:
    path "mitogenome.liftoff.gff3",        emit: gff
    path "liftoff_unmapped_features.txt",  emit: unmapped
    path "versions.yml",                   emit: versions

    script:
    """
    liftoff -g ${ref_gff} -o mitogenome.liftoff.gff3 -u liftoff_unmapped_features.txt ${mito} ${ref}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: "\$(liftoff --version 2>&1 | tail -1 | sed 's/^v//')"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mitogenome.liftoff.gff3 liftoff_unmapped_features.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_MITOS2
// Independent de novo annotation (MITOS2 and its RefSeq 89 metazoan data,
// both shipped in the MitoHiFi image).
// ============================================================
process MITO_MITOS2 {
    label 'mito_env'
    label 'res_light'
    tag   "mito_mitos2"

    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'
    publishDir "${params.outdir}/05_mitogenome/final", mode: params.publish_mode, overwrite: true

    input:
    path mito

    output:
    path "mitos2",       emit: dir
    path "versions.yml", emit: versions

    script:
    """
    unset LD_PRELOAD
    mkdir -p mitos2
    runmitos.py --noplots -i ${mito} -c ${params.mito_genetic_code} --refdir /opt/databases -r refseq89m -o mitos2
    test -s mitos2/result.gff
    mv mitos2/result.gff mitos2/mitogenome.mitos2.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitos: "2 (MitoHiFi image, refseq89m)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p mitos2 && touch mitos2/mitogenome.mitos2.gff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitos: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MITO_GENE_CHECK
// Translate every lifted CDS with the mitochondrial genetic code and compare
// it to the reference protein: same length and no internal stop = OK.
// This is the only self-check available without Illumina data; it cannot
// see errors outside protein-coding genes.
// ============================================================
process MITO_GENE_CHECK {
    label 'mito_env'
    label 'res_tiny'
    tag   "mito_gene_check"

    container 'quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'
    publishDir "${params.outdir}/05_mitogenome/final", mode: params.publish_mode, overwrite: true

    input:
    path mito
    path gff
    path ref,     stageAs: 'annotation_ref.fa'
    path ref_gff, stageAs: 'annotation_ref.gff3'

    output:
    path "mito_gene_check.tsv",     emit: tsv
    path "mitogenome.proteins.faa", emit: proteins
    path "versions.yml",            emit: versions

    script:
    """
    agat_sp_extract_sequences.pl -g ${gff}     -f ${mito} -p --table ${params.mito_genetic_code} --keep_attributes -o mitogenome.proteins.faa > agat_lifted.log 2>&1
    agat_sp_extract_sequences.pl -g ${ref_gff} -f ${ref}  -p --table ${params.mito_genetic_code} --keep_attributes -o ref.faa > agat_ref.log 2>&1

    # one line per protein: feature id, gene name, length without the final stop, internal stops
    tab() {
        awk '/^>/{ if (id != "") print id "\\t" name "\\t" seq
                   id = substr(\$1, 2); name = id; seq = ""
                   for (i = 2; i <= NF; i++) if (\$i ~ /^gene=/ && \$i !~ /^gene=agat-/) name = substr(\$i, 6)
                   next }
             { seq = seq \$0 }
             END { if (id != "") print id "\\t" name "\\t" seq }' "\$1" \\
        | awk -F'\\t' -v OFS='\\t' '{ s = \$3; sub(/\\*\$/, "", s); n = gsub(/\\*/, "*", s); print \$1, \$2, length(s), n }'
    }
    tab mitogenome.proteins.faa > lifted.tsv
    tab ref.faa > ref.tsv

    awk -F'\\t' -v OFS='\\t' '
        NR == FNR { len[\$1] = \$3; stops[\$1] = \$4; next }
        FNR == 1  { print "gene", "aa_length", "ref_aa_length", "internal_stops", "status" }
        { if (!(\$1 in len))                          print \$2, "NA", \$3, "NA", "MISSING"
          else if (stops[\$1] > 0 || len[\$1] != \$3) print \$2, len[\$1], \$3, stops[\$1], "BROKEN"
          else                                        print \$2, len[\$1], \$3, stops[\$1], "OK" }' lifted.tsv ref.tsv > mito_gene_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "\$(agat_sp_extract_sequences.pl --version 2>&1 | grep -oim1 'v[0-9.]*' | sed 's/v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch mito_gene_check.tsv mitogenome.proteins.faa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Remove the contigs MitoHiFi identified as mitochondrial from the NUCLEAR assembly
// (bin/remove_mito_contigs.py), so the mitogenome is not released twice. Listed contigs
// longer than params.mito_contig_max_ratio x the reference mitogenome are kept.
// ============================================================
process REMOVE_MITO_CONTIGS {
    label 'res_light'
    tag   "remove_mito_contigs"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/05_mitogenome/nuclear_filter", mode: params.publish_mode, pattern: '*.tsv'

    input:
    path assembly
    path contig_ids,   stageAs: 'mitohifi/*'   // or NO_ placeholder: nothing found
    path contig_stats, stageAs: 'mitohifi/*'
    path mito_ref

    output:
    path "nuclear_contigs.fasta",    emit: fasta
    path "mito_contigs_removed.tsv", emit: report
    path "versions.yml",             emit: versions

    script:
    """
    remove_mito_contigs.py ${assembly} ${contig_ids} ${contig_stats} ${mito_ref} ${params.mito_contig_max_ratio} \\
        nuclear_contigs.fasta mito_contigs_removed.tsv

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    cp ${assembly} nuclear_contigs.fasta
    printf 'contig\\tlength\\tcircular\\tgenes\\taction\\n' > mito_contigs_removed.tsv
    printf '"${task.process}":\\n    python: stub\\n' > versions.yml
    """
}
