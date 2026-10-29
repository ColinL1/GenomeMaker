// ============================================================
// 07_structural_annotation.nf
// Structural gene annotation pipeline:
// tRNAscan-SE → RNA-seq trimming → STAR mapping → BRAKER4 → AGAT merge
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: TRNASCAN_SE
// Predict tRNA genes using tRNAscan-SE
// ============================================================
process TRNASCAN_SE {
    label 'structural_env'
    tag   "trnascan"

    conda 'bioconda::trnascan-se'
    container 'quay.io/biocontainers/trnascan-se:2.0.13--pl5321hab16a5f_0'
    publishDir "${params.outdir}/07_structural_annotation/trnascan", mode: params.publish_mode

    input:
    path assembly

    output:
    path "trnascan.out",  emit: trnascan_out
    path "trnascan.tbl",  emit: trnascan_tbl
    path "trnascan.gff",  emit: trnascan_gff
    path "trnascan.log",  emit: trnascan_log
    path "versions.yml", emit: versions

    script:
    // -j/--gff is tRNAscan-SE's own GFF3 writer (verified present in 2.0.13),
    // used instead of converting the -o tabular output.
    """
    tRNAscan-SE -E -I -H --detail --thread ${task.cpus} \\
      -o trnascan.out \\
      -f trnascan.tbl \\
      -j trnascan.gff \\
      -m trnascan.log \\
      ${assembly}

    # A genome with zero tRNAs is not plausible for a metazoan assembly, so an
    # empty GFF means the run failed rather than legitimately finding nothing.
    if [ ! -s trnascan.gff ]; then
        echo "ERROR: tRNAscan-SE produced an empty GFF3 for ${assembly}" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trnascan-se: "\$(tRNAscan-SE -h 2>&1 | grep -im1 'tRNAscan-SE [0-9]' | sed 's/tRNAscan-SE //;s/ .*//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch trnascan.out trnascan.tbl trnascan.log
    printf '##gff-version 3\\nstub\\ttRNAscan-SE\\ttRNA\\t1\\t72\\t.\\t-\\t.\\tID=stub_tRNA\\n' > trnascan.gff
    echo "TRNASCAN_SE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trnascan-se: "stub"
    END_VERSIONS
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

    conda 'bioconda::trimmomatic'
    container 'quay.io/biocontainers/trimmomatic:0.41--hdfd78af_0'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/trimmed", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)
    path trimmomatic_adapters

    output:
    tuple val(sample_id), path("${sample_id}_1P.fastq.gz"), path("${sample_id}_2P.fastq.gz"), emit: trimmed_reads
    path "versions.yml", emit: versions

    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    export _JAVA_OPTIONS="-Xmx${avail_mem}M"
    trimmomatic PE -threads ${task.cpus} \\
      ${reads[0]} ${reads[1]} \\
      ${sample_id}_1P.fastq.gz ${sample_id}_1U.fastq.gz \\
      ${sample_id}_2P.fastq.gz ${sample_id}_2U.fastq.gz \\
      ILLUMINACLIP:${trimmomatic_adapters}:2:30:10 \\
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: "\$(trimmomatic -version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_1P.fastq.gz ${sample_id}_2P.fastq.gz
    echo "RNASEQ_TRIM stub completed for ${sample_id}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: STAR_GENOME_INDEX
// Generate STAR genome index from assembly
// ============================================================
process STAR_GENOME_INDEX {
    label 'structural_env'
    tag   "star_index"

    conda 'bioconda::star'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/star_index", mode: params.publish_mode

    input:
    path assembly

    output:
    path "star_index", emit: star_index_dir
    path "versions.yml", emit: versions

    script:
    // --limitGenomeGenerateRAM must track the ALLOCATION, not a hardcoded value.
    // It was 100 GiB against a 64 GB request, so STAR believed it had 36 GB more
    // than the scheduler would ever give it -- an OOM kill waiting to happen.
    // 85% leaves headroom for STAR's non-buffer allocations.
    def gen_ram = (task.memory.toBytes() * 0.85) as long
    """
    # --genomeSAindexNbases must scale with genome size: STAR's own formula is
    # min(14, log2(genomeLength)/2 - 1). The previous hardcoded 10 is a value for
    # ~4 Mb genomes; for a ~600 Mb coral the correct value is 14, and too small a
    # value gives a degraded suffix-array index -- slower mapping and worse memory
    # behaviour. Computed from the actual assembly so it stays right if the
    # assembly size changes.
    GLEN=\$(grep -v '^>' ${assembly} | tr -d '\\n' | wc -c)
    SA_N=\$(awk -v g="\$GLEN" 'BEGIN{n=int(log(g)/log(2)/2-1); if(n>14)n=14; if(n<4)n=4; print n}')
    echo "genome length \$GLEN bp -> --genomeSAindexNbases \$SA_N"

    STAR --runThreadN ${task.cpus} \\
      --runMode genomeGenerate \\
      --genomeDir star_index \\
      --genomeFastaFiles ${assembly} \\
      --genomeSAindexNbases \$SA_N \\
      --limitGenomeGenerateRAM ${gen_ram}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: "\$(STAR --version 2>&1 | tail -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SAindex
    echo "STAR_GENOME_INDEX stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: "stub"
    END_VERSIONS
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

    conda 'bioconda::star'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/star_mapping", mode: params.publish_mode

    input:
    path star_index
    tuple val(sample_id), path(r1p), path(r2p)

    output:
    path "${sample_id}.sorted.bam",        emit: sample_bam
    path "${sample_id}_Log.final.out",     emit: star_log
    path "versions.yml", emit: versions

    script:
    // --limitBAMsortRAM was 50 GB against a 48 GB request, and the loaded genome
    // index needs memory on top of the sort buffer -- so STAR could exceed the
    // allocation on its own accounting. Give the sort buffer half the allocation
    // and let the index have the rest.
    def sort_ram = (task.memory.toBytes() * 0.5) as long
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
      --limitBAMsortRAM ${sort_ram}

    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: "\$(STAR --version 2>&1 | tail -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}_Log.final.out
    echo "STAR_MAPPING stub completed for ${sample_id}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: "stub"
    END_VERSIONS
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

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/07_structural_annotation/rnaseq/merged", mode: params.publish_mode

    input:
    path bams

    output:
    path "merged_rnaseq.bam", emit: rnaseq_bam
    path "merged_rnaseq.bam.bai", emit: rnaseq_bai
    path "versions.yml", emit: versions

    script:
    """
    samtools merge -@ ${task.cpus} merged_rnaseq.bam ${bams}
    samtools index merged_rnaseq.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch merged_rnaseq.bam merged_rnaseq.bam.bai
    echo "MERGE_BAMS stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BRAKER4_RUN
// Structural gene prediction with BRAKER4 (Snakemake re-implementation of BRAKER3)
//
// Why not BRAKER3: braker.pl 3.0.8 (teambraker/braker3 v3.1.1 and v3.0.7.6)
// cannot use --busco_lineage: its compleasm_to_hints.py forces an
// '_odb10' folder name, while the bundled compleasm 0.2.8 only uses odb12
// (and odb10 datasets are no longer downloadable). BRAKER4 handles odb12.
//
// BRAKER4 is a Snakemake workflow that launches its own containers, so this
// process has NO container: it runs on the host with the Snakemake env in
// params.braker4_snakemake_env and apptainer on PATH. BRAKER4 itself is pinned
// by params.braker4_dir (a checkout of Gaius-Augustus/BRAKER4 at a fixed commit).
// Deliberately not a `conda` directive: that needs conda.enabled=true, which
// changes the resume hash of every containerised task in the pipeline.
//
// RNA-seq: the merged STAR BAM goes in as is. BRAKER4 has no --stranded option;
// the transcript strand comes from STAR's XS tags (--outSAMstrandField intronMotif).
// ============================================================
process BRAKER4_RUN {
    label 'structural_env'
    tag   "braker4"

    publishDir "${params.outdir}/07_structural_annotation/braker", mode: params.publish_mode

    input:
    path assembly        // soft-masked (RepeatMasker -xsmall)
    path rnaseq_bam      // merged, coordinate-sorted STAR BAM
    path protein_refs

    output:
    path "braker/braker.gtf", emit: braker_gtf
    path "braker/braker.log", emit: braker_log
    path "braker/results",    emit: braker_results   // BRAKER4 results dir (UTR gtf, gff3, QC, report)
    path "versions.yml", emit: versions

    script:
    def sample = params.species_name
    """
    set -euo pipefail
    export PATH="${params.braker4_snakemake_env}/bin:\$PATH"
    BRAKER4_DIR="${params.braker4_dir}"

    GENOME=\$(readlink -f ${assembly})
    BAM=\$(readlink -f ${rnaseq_bam})
    # GeneMark-ETP aborts on a repeated protein ID. coral_proteins_ncbi.faa has
    # 29,556 IDs present twice, all with identical sequences, so keep the first
    # record per ID (the source file is left untouched).
    awk '/^>/{ keep = !seen[\$1]++ } keep' ${protein_refs} > proteins_dedup.faa
    PROT=\$(readlink -f proteins_dedup.faa)

    # The soft-masked assembly is given as both genome and genome_masked: it is
    # the same sequence (repeats lowercase), and genome_masked is what gene
    # prediction uses.
    cat > samples.csv <<-EOF
    sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
    ${sample},\$GENOME,\$GENOME,\$PROT,\$BAM,,,,,,,,${params.busco_lineage},
    EOF

    cat > config.ini <<-EOF
    [paths]
    samples_file = samples.csv
    augustus_config_path = augustus_config
    [PARAMS]
    fungus = 0
    min_contig = ${params.braker4_min_contig}
    use_varus = 0
    skip_optimize_augustus = 0
    use_compleasm_hints = 1
    skip_busco = 0
    run_omark = 0
    run_ncrna = 0
    no_cleanup = 0
    EOF

    # Bind every directory the containers must see: the task dir, the real
    # locations of the staged inputs, and the BRAKER4 checkout (download cache).
    BINDS=\$(printf '%s\\n' "\$PWD" "\$(dirname \$GENOME)" "\$(dirname \$BAM)" "\$(dirname \$PROT)" "\$BRAKER4_DIR" | sort -u | sed 's/^/-B /' | tr '\\n' ' ')

    snakemake -s "\$BRAKER4_DIR/Snakefile" \\
      --cores ${task.cpus} \\
      --use-apptainer \\
      --apptainer-prefix "\${NXF_APPTAINER_CACHEDIR:-\$PWD/apptainer_cache}" \\
      --apptainer-args "\$BINDS" \\
      --latency-wait 120

    RES=output/${sample}/results
    mkdir -p braker
    gunzip -c \$RES/braker.gtf.gz > braker/braker.gtf
    cp -r \$RES braker/results
    cat .snakemake/log/*.snakemake.log > braker/braker.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker4: "\$(git -C "\$BRAKER4_DIR" describe --tags --always 2>/dev/null || echo unknown)"
        snakemake: "\$(snakemake --version)"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p braker/results
    touch braker/braker.gtf braker/braker.log
    echo "BRAKER4_RUN stub completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker4: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MERGE_ANNOTATIONS
// Merge BRAKER GTF, tRNA GFF, and run AGAT fix/merge
// ============================================================
process MERGE_ANNOTATIONS {
    label 'structural_env'
    label 'res_mid_6h'
    tag   "merge_annotations"

    conda 'bioconda::agat'
    container 'quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'
    publishDir "${params.outdir}/07_structural_annotation/merged_annotation", mode: params.publish_mode

    input:
    path braker_gtf
    path trnascan_gff
    path ncrna_tblout   // Infernal cmscan --fmt 2 tblout, or a NO_ placeholder (Rfam off)
    path rfam_family    // Rfam family.txt (per-family type), or a NO_ placeholder

    output:
    path "${params.species_name}.v${params.version}.gff3", emit: gff3
    path "versions.yml", emit: versions

    script:
    // Output names are derived from species_name and version.
    def gff_out = "${params.species_name}.v${params.version}.gff3"
    def use_rfam = !ncrna_tblout.name.startsWith('NO_') && !rfam_family.name.startsWith('NO_')
    """
    # GTF -> GFF3. (agat_sp_convert_format.pl does not exist in AGAT;
    # agat_convert_sp_gxf2gxf.pl is its GTF/GFF -> GFF3 converter.)
    agat_convert_sp_gxf2gxf.pl -g ${braker_gtf} -o braker.gff3

    # Non-coding RNAs (bin/ncrna_to_gff3.pl):
    #  - tRNAscan-SE: one explicit gene per tRNA. Without it AGAT groups tRNAs
    #    under shared genes and funannotate keeps one tRNA per gene.
    #  - Infernal/Rfam: gene/ncRNA|rRNA/exon for Rfam "Gene" families (not tRNA,
    #    not ribozymes), typed by Rfam family.txt; rRNAs come from here, not barrnap.
    ncrna_to_gff3.pl trnascan ${trnascan_gff} > trnascan_genes.gff3
    ${use_rfam ? "ncrna_to_gff3.pl rfam ${ncrna_tblout} ${rfam_family} > rfam_ncrna.gff3" : "# Rfam ncRNAs: skipped (no Infernal output / family.txt)"}

    # Merge gene models with the ncRNA annotations
    agat_sp_merge_annotations.pl \\
      --gff braker.gff3 \\
      --gff trnascan_genes.gff3 \\
      ${use_rfam ? '--gff rfam_ncrna.gff3' : ''} \\
      --out merged.gff

    # Resolve overlapping gene models
    agat_sp_fix_overlaping_genes.pl \\
      -f merged.gff \\
      -o fixed.gff3

    # Normalise/sanity-pass with AGAT's own converter. The previous version called
    # `gt gff3validator` here, but GenomeTools is not in the AGAT container, so that
    # line could only fail. Strict spec validation now lives in GT_GFF3_VALIDATE,
    # which has its own container.
    agat_convert_sp_gxf2gxf.pl -g fixed.gff3 -o normalised.gff3

    # BRAKER's protein-coding models come through as type 'transcript'. The SO
    # type for a transcript with a CDS is mRNA, and funannotate (util stats,
    # annotate) only counts mRNA children of genes -- with 'transcript' it saw
    # zero proteins. tRNA/ncRNA features keep their own types.
    awk -F'\\t' 'BEGIN{OFS="\\t"} !/^#/ && \$3=="transcript"{\$3="mRNA"} {print}' \\
      normalised.gff3 > ${gff_out}

    # Protein extraction (gffread) is in EXTRACT_PROTEINS: gffread is not in the
    # AGAT container.

    if [ ! -s "${gff_out}" ]; then
        echo "ERROR: ${gff_out} is empty after annotation merge" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "\$(agat_sp_merge_annotations.pl --version 2>&1 | grep -oim1 'v[0-9.]*' | sed 's/v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > ${params.species_name}.v${params.version}.gff3
    echo "MERGE_ANNOTATIONS stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: EXTRACT_PROTEINS
// Protein FASTA from the merged GFF3 (consumed by functional annotation).
// Separate from MERGE_ANNOTATIONS because gffread is not in the AGAT container.
// ============================================================
process EXTRACT_PROTEINS {
    label 'structural_env'
    label 'res_light'
    tag   "extract_proteins"

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'
    publishDir "${params.outdir}/07_structural_annotation/merged_annotation", mode: params.publish_mode

    input:
    path gff3
    path assembly

    output:
    path "${params.species_name}.proteins.faa", emit: proteins_faa
    path "versions.yml", emit: versions

    script:
    """
    gffread ${gff3} \\
      -g ${assembly} \\
      -y ${params.species_name}.proteins.faa

    if [ ! -s "${params.species_name}.proteins.faa" ]; then
        echo "ERROR: ${params.species_name}.proteins.faa is empty" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "\$(gffread --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${params.species_name}.proteins.faa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Workflow composition (STRUCTURAL_ANNOTATION) now lives in
// ../subworkflows/07_structural_annotation.nf — this file only
// defines the individual processes.
// ============================================================

// ============================================================
// Process: GT_GFF3_VALIDATE
// Strict GFF3 spec validation with GenomeTools.
//
// WHY ITS OWN PROCESS: the AGAT container does not ship GenomeTools. This keeps
// one tool per container, which is how the rest of the pipeline is built.
// ============================================================
process GT_GFF3_VALIDATE {
    label 'structural_env'
    tag   "gff3_validate"

    conda 'bioconda::genometools-genometools=1.6.6'
    container 'quay.io/biocontainers/genometools-genometools:1.6.6--py311h21ec246_1'
    publishDir "${params.outdir}/07_structural_annotation/merged_annotation", mode: params.publish_mode

    input:
    path gff3

    output:
    path "gff3_validation.txt", emit: report
    path "versions.yml", emit: versions

    script:
    """
    # gt gff3validator exits non-zero on a spec violation, which should fail the
    # run: a malformed GFF3 propagates into functional annotation and submission.
    gt gff3validator ${gff3} > gff3_validation.txt 2>&1 || {
        echo "ERROR: ${gff3} failed GFF3 spec validation:" >&2
        cat gff3_validation.txt >&2
        exit 1
    }
    echo "GFF3 validation passed: ${gff3}" >> gff3_validation.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: "\$(gt -version 2>&1 | head -1 | sed 's/gt (GenomeTools) //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    echo "GT_GFF3_VALIDATE stub completed" > gff3_validation.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: "stub"
    END_VERSIONS
    """
}
