// ============================================================
// 07_structural_annotation.nf
// Structural gene annotation pipeline:
// tRNAscan-SE → RNA-seq trimming → STAR mapping → BRAKER3 → AGAT merge
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: TRNASCAN_SE
// Predict tRNA genes using tRNAscan-SE
// ============================================================
process TRNASCAN_SE {
    label 'structural_env'
    tag   "trnascan"
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::trnascan-se'
    container 'quay.io/biocontainers/trnascan-se:2.0.13--pl5321hab16a5f_0'
    publishDir "${params.outdir}/07_structural_annotation/trnascan", mode: params.publish_mode

    input:
    path assembly

    output:
    path "trnascan.out", emit: trnascan_out
    path "trnascan.tbl", emit: trnascan_tbl
    path "trnascan.gff", emit: trnascan_gff

    script:
    """
    tRNAscan-SE -E -I -H --detail --thread ${task.cpus} \\
      -o trnascan.out \\
      -f trnascan.tbl \\
      -m trnascan.log \\
      ${assembly}

    # Convert to GFF3 format using a separate Perl script
    perl ${params.bin_path}trnascan_to_gff3.pl
    """

    stub:
    """
    touch trnascan.out trnascan.tbl trnascan.gff
    echo "TRNASCAN_SE stub completed"
    """
}

// ============================================================
// Process: RNASEQ_TRIM
// Trim one paired-end RNA-seq library with Trimmomatic.
// Forked per-sample: the calling workflow supplies a channel of
// [sample_id, [R1, R2]] tuples (e.g. via channel.fromFilePairs),
// so Nextflow parallelizes across libraries instead of looping
// inside a single task.
// ============================================================
process RNASEQ_TRIM {
    label 'structural_env'
    tag   "rnaseq_trim_${sample_id}"
    cpus  { params.threads_default * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::trimmomatic'
    container 'quay.io/biocontainers/trimmomatic:0.41--hdfd78af_0'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/trimmed", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1P.fastq.gz"), path("${sample_id}_2P.fastq.gz"), emit: trimmed_reads

    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    export _JAVA_OPTIONS="-Xmx${avail_mem}M"
    trimmomatic PE -threads ${task.cpus} \\
      ${reads[0]} ${reads[1]} \\
      ${sample_id}_1P.fastq.gz ${sample_id}_1U.fastq.gz \\
      ${sample_id}_2P.fastq.gz ${sample_id}_2U.fastq.gz \\
      ILLUMINACLIP:${params.trimmomatic_adapters}:2:30:10 \\
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    """

    stub:
    """
    touch ${sample_id}_1P.fastq.gz ${sample_id}_2P.fastq.gz
    echo "RNASEQ_TRIM stub completed for ${sample_id}"
    """
}

// ============================================================
// Process: STAR_GENOME_INDEX
// Generate STAR genome index from assembly
// ============================================================
process STAR_GENOME_INDEX {
    label 'structural_env'
    tag   "star_index"
    cpus  { params.threads_star * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::star'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/star_index", mode: params.publish_mode

    input:
    path assembly

    output:
    path "star_index", emit: star_index_dir

    script:
    """
    STAR --runThreadN ${task.cpus} \\
      --runMode genomeGenerate \\
      --genomeDir star_index \\
      --genomeFastaFiles ${assembly} \\
      --genomeSAindexNbases 10 \\
      --limitGenomeGenerateRAM 107374182400
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SAindex
    echo "STAR_GENOME_INDEX stub completed"
    """
}

// ============================================================
// Process: STAR_MAPPING
// Map one trimmed paired-end RNA-seq library to the genome with
// STAR. Forked per-sample (paired with the shared STAR index via
// `.first()` in the calling workflow); outputs are combined and
// merged downstream by MERGE_BAMS.
// ============================================================
process STAR_MAPPING {
    label 'structural_env'
    tag   "star_mapping_${sample_id}"
    cpus  { params.threads_star * task.attempt }
    memory { 48.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::star'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/star_mapping", mode: params.publish_mode

    input:
    path star_index
    tuple val(sample_id), path(r1p), path(r2p)

    output:
    path "${sample_id}.sorted.bam", emit: sample_bam

    script:
    """
    STAR \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${star_index} \\
      --readFilesIn ${r1p} ${r2p} \\
      --readFilesCommand zcat \\
      --outSAMtype BAM SortedByCoordinate \\
      --outSAMstrandField intronMotif \\
      --twopassMode Basic \\
      --outFileNamePrefix ${sample_id}_ \\
      --limitBAMsortRAM 50000000000

    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam
    echo "STAR_MAPPING stub completed for ${sample_id}"
    """
}

// ============================================================
// Process: MERGE_BAMS
// Combine all per-sample STAR BAMs (collected from the forked
// STAR_MAPPING tasks) into a single merged, indexed BAM.
// ============================================================
process MERGE_BAMS {
    label 'structural_env'
    tag   "merge_rnaseq_bams"
    cpus  { params.threads_star * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/merged", mode: params.publish_mode

    input:
    path bams

    output:
    path "merged_rnaseq.bam", emit: rnaseq_bam
    path "merged_rnaseq.bam.bai", emit: rnaseq_bai

    script:
    """
    samtools merge -@ ${task.cpus} merged_rnaseq.bam ${bams}
    samtools index merged_rnaseq.bam
    """

    stub:
    """
    touch merged_rnaseq.bam merged_rnaseq.bam.bai
    echo "MERGE_BAMS stub completed"
    """
}

// ============================================================
// Process: SPLIT_STRANDED_BAM
// Split merged BAM by read direction for stranded RNA-seq
// ============================================================
process SPLIT_STRANDED_BAM {
    label 'structural_env'
    tag   "split_strand"
    cpus  { 4 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/stranded", mode: params.publish_mode

    input:
    path bam

    output:
    path "plus_strand.bam",  emit: plus_bam
    path "minus_strand.bam", emit: minus_bam

    script:
    """
    # For dUTP strandedness (most common):
    #   Read 1 (R1) on plus strand: flags 99, 83, 147, 163
    #   Read 1 (R1) on minus strand: flags 147, 163
    # We split by the expected strand of the transcript

    # Plus strand transcripts (R1 mapped to minus strand of genome)
    samtools view -h ${bam} | \\
      awk 'BEGIN{OFS="\\t"} /^@/ || (\$2==147 || \$2==163)' | \\
      samtools view -b -o plus_strand.bam

    # Minus strand transcripts (R1 mapped to plus strand of genome)
    samtools view -h ${bam} | \\
      awk 'BEGIN{OFS="\\t"} /^@/ || (\$2==99 || \$2==83)' | \\
      samtools view -b -o minus_strand.bam
    """

    stub:
    """
    touch plus_strand.bam minus_strand.bam
    echo "SPLIT_STRANDED_BAM stub completed"
    """
}

// ============================================================
// Process: BRAKER3_RUN
// Structural gene prediction with BRAKER3
// ============================================================
process BRAKER3_RUN {
    label 'structural_env'
    tag   "braker3"
    cpus  { params.threads_braker * task.attempt }
    memory { 128.GB * task.attempt }
    time   { 72.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::braker'
    container 'teambraker/braker3:v3.1.1'
    publishDir "${params.outdir}/07_structural_annotation/braker", mode: params.publish_mode

    input:
    path assembly
    val  rnaseq_strandedness
    path plus_bam
    path minus_bam
    path protein_refs

    output:
    path "braker/braker.gtf", emit: braker_gtf
    path "braker/braker.log", emit: braker_log

    script:
    def strand_spec = rnaseq_strandedness.toLowerCase()
    def bam_spec = plus_bam && minus_bam
        ? "--bam=${plus_bam},${minus_bam}"
        : "--bam=${plus_bam ?: minus_bam}"
    """
    # BRAKER3: Gene prediction using RNA-seq and protein hints
    braker.pl \\
      --species=${params.species_name} \\
      --genome=${assembly} \\
      ${bam_spec} \\
      --stranded=${strand_spec} \\
      --threads ${task.cpus} \\
      --prot_seq=${protein_refs} \\
      --busco_lineage=${params.busco_lineage_meta} \\
      --softmasking

    # Copy outputs for easy access
    mkdir -p braker
    cp braker.gtf braker/braker.gtf 2>/dev/null || true
    cp braker.log braker/braker.log 2>/dev/null || true
    """

    stub:
    """
    mkdir -p braker
    touch braker/braker.gtf braker/braker.log
    echo "BRAKER3_RUN stub completed"
    """
}

// ============================================================
// Process: MERGE_ANNOTATIONS
// Merge BRAKER3 GTF, tRNA GFF, and run AGAT fix/merge
// ============================================================
process MERGE_ANNOTATIONS {
    label 'structural_env'
    tag   "merge_annotations"
    cpus  { 4 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::agat bioconda::gffread'
    container 'quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'
    publishDir "${params.outdir}/07_structural_annotation/merged_annotation", mode: params.publish_mode

    input:
    path braker_gtf
    path trnascan_gff
    path assembly

    output:
    path "PAG_${params.species_name}_1.1.gff3", emit: gff3
    path "Phar.braker.prot.fasta",            emit: proteins_faa

    script:
    """
    # Convert GTF to GFF3
    agat_sp_convert_format.pl -g ${braker_gtf} -f gtf -o braker.gff3

    # Merge annotations
    agat_sp_merge_annotations.pl \\
      --gff braker.gff3 \\
      --gff ${trnascan_gff} \\
      --out merged.gff

    # Fix overlapping genes
    agat_sp_fix_overlaping_genes.pl \\
      -f merged.gff \\
      -o PAG_${params.species_name}_1.1.gff3

    # Validate GFF3
    gt gff3validator PAG_${params.species_name}_1.1.gff3

    # Extract proteins
    gffread PAG_${params.species_name}_1.1.gff3 \\
      -g ${assembly} \\
      -y Phar.braker.prot.fasta
    """

    stub:
    """
    touch PAG_${params.species_name}_1.1.gff3 Phar.braker.prot.fasta
    echo "MERGE_ANNOTATIONS stub completed"
    """
}

// ============================================================
// Workflow composition (STRUCTURAL_ANNOTATION) now lives in
// ../subworkflows/07_structural_annotation.nf — this file only
// defines the individual processes.
// ============================================================
