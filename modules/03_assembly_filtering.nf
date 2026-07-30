// ============================================================
// 03_assembly_filtering.nf
// BlobToolKit-based assembly filtering
// Removes contamination and splits out symbiont contigs based on
// BLASTN/DIAMOND taxonomy hits, ONT coverage, and BUSCO evidence
//
// NOTE ON VERSIONS: this uses BlobToolKit (BTK) v3/4 CLI syntax
// (`blobtools create/add/filter` against a BlobDir folder), NOT
// legacy BlobTools v1. Pin the container/conda package accordingly.
//
// NOTE ON CLI FLAGS: `blobtools filter` category params operate as
// EXCLUDE lists by default (--param FIELD--Keys=Value drops matching
// records); add --invert to flip that into an INCLUDE ("keep only
// this taxon") filter, which is what we want here. Flag names have
// shifted across BTK releases -- if you're on a different BTK
// version than the one pinned below, run `blobtools filter --help`
// inside the container and confirm --param/--Keys/--invert/--output/
// --suffix still match before trusting this at scale.
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: BLAST_AND_COV
// Build a BlobDir with taxonomy hits, ONT coverage, and BUSCO evidence
// ============================================================
process BLAST_AND_COV {
    label 'blobtools_env'
    tag   "blobtools_build"
    cpus  { params.threads_default * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::blobtoolkit bioconda::blast+ bioconda::diamond bioconda::minimap2 bioconda::samtools'
    container 'genomehubs/blobtoolkit:4.4.4'  // verify current tag against your install
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode

    input:
    path assembly
    path reads             // ONT long reads (assembly_reads / porechopped)
    path busco_table       // single BUSCO full_table.tsv, anthozoa_odb12
    path meta_yaml         // sample metadata YAML
    val  taxid             // NCBI taxid of the host (anchors lineage in BlobDir, not a filter)

    output:
    path "blastn.out",         emit: blastn_out
    path "diamond.blastx.out", emit: blastx_out
    path "coverage.bam",       emit: cov_bam
    path "btk",                emit: blobdir

    script:
    """
    # BLASTN vs NCBI nt
    blastn -db /db/nt/nt \\
      -query ${assembly} \\
      -outfmt "6 qseqid staxids bitscore std" \\
      -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 \\
      -num_threads ${task.cpus} \\
      -out blastn.out

    # DIAMOND blastx vs UniProt reference proteomes
    diamond blastx \\
      --query ${assembly} \\
      --db /db/uniprot/reference_proteomes.dmnd \\
      --outfmt 6 qseqid staxids bitscore \\
      --sensitive --max-target-seqs 1 --evalue 1e-25 \\
      --threads ${task.cpus} \\
      > diamond.blastx.out

    # Coverage via minimap2 -- ONT long reads, NOT short-read preset
    minimap2 -a -x map-ont -t ${task.cpus} ${assembly} ${reads} | \\
      samtools sort -@ ${task.cpus} -O BAM -o coverage.bam -
    samtools index coverage.bam

    # Build BlobDir
    blobtools create \\
      --fasta ${assembly} \\
      --meta ${meta_yaml} \\
      --taxid ${taxid} \\
      btk

    # Add taxonomy hits
    blobtools add \\
      --hits diamond.blastx.out \\
      --hits blastn.out \\
      --taxrule bestsumorder \\
      btk

    # Add coverage
    blobtools add --cov coverage.bam btk

    # Add BUSCO (host-specific: anthozoa_odb12) -- flags host contigs unambiguously
    blobtools add --busco ${busco_table} btk
    """

    stub:
    """
    touch blastn.out diamond.blastx.out coverage.bam
    mkdir -p btk
    echo "BLAST_AND_COV stub completed"
    """
}

// ============================================================
// Process: GENERATE_TAXON_FILTER
// Generate a non-interactive "keep this taxon" filter JSON + summary
// Fans out over (category, taxon_name) pairs, e.g. ('host','Cnidaria')
// ============================================================
process GENERATE_TAXON_FILTER {
    label 'blobtools_env'
    tag   "filter_${category}"
    cpus  { 2 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 1.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode

    input:
    tuple val(category), val(taxon_name), path(blobdir)

    output:
    tuple val(category), path("${category}_filter.json"), emit: filter_json
    tuple val(category), path("${category}_summary.json"), emit: summary

    script:
    """
    # Keep-only-this-taxon filter: exclude everything NOT matching
    # bestsumorder_phylum == taxon_name, then --invert to make it
    # an inclusive "keep" filter instead of an exclusion.
    blobtools filter \\
      --param bestsumorder_phylum--Keys=${taxon_name} \\
      --invert \\
      --json ${category}_filter.json \\
      --summary ${category}_summary.json \\
      ${blobdir}
    """

    stub:
    """
    touch ${category}_filter.json ${category}_summary.json
    echo "GENERATE_TAXON_FILTER stub completed (${category})"
    """
}

// ============================================================
// Process: BLOBTOOLS_APPLY_FILTER
// Apply a category's filter JSON to produce a category-specific FASTA
// Each category writes to its own --output dir, so no filename
// collisions across host/symbiont/etc runs of the same process.
// ============================================================
process BLOBTOOLS_APPLY_FILTER {
    label 'blobtools_env'
    tag   "apply_filter_${category}"
    cpus  { 4 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    publishDir "${params.outdir}/03_filtering/filtered", mode: params.publish_mode

    input:
    tuple val(category), path(filter_json), path(assembly), path(blobdir)

    output:
    tuple val(category), path("${category}_assembly.fasta"), emit: filtered_assembly

    script:
    """
    blobtools filter \\
      --json ${filter_json} \\
      --fasta ${assembly} \\
      --output ${category}_filtered \\
      --suffix ${category} \\
      ${blobdir}

    # Each category has its own --output dir, so exactly one fasta is
    # expected here -- no cross-category glob/overwrite risk.
    find ${category}_filtered -name "*.fasta" -exec cp {} ${category}_assembly.fasta \\;

    if [ ! -s ${category}_assembly.fasta ]; then
      echo "ERROR: no filtered fasta produced for category '${category}'" >&2
      exit 1
    fi
    """

    stub:
    """
    touch ${category}_assembly.fasta
    echo "BLOBTOOLS_APPLY_FILTER stub completed (${category})"
    """
}

// ============================================================
// Process: PASSTHROUGH
// Pass assembly through unchanged (when blobtools is skipped)
// ============================================================
process PASSTHROUGH {
    tag   "passthrough"
    cpus  { 1 * task.attempt }
    memory { 4.GB * task.attempt }
    time   { 1.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/03_filtering/filtered", mode: params.publish_mode

    input:
    path assembly

    output:
    path "filtered_assembly.fasta", emit: filtered_assembly

    script:
    """
    cp ${assembly} filtered_assembly.fasta
    """

    stub:
    """
    touch filtered_assembly.fasta
    echo "PASSTHROUGH stub completed"
    """
}


// // ============================================================
// // Workflow: assembly_filtering
// // BlobTools filtering or passthrough
// // ============================================================
// workflow assembly_filtering {

//     take:
//         assembly
//         reads
//         busco_table
//         meta_yaml

//     main:
//         if (params.skip_blobtools) {
//             // Skip filtering, pass assembly through
//             PASSTHROUGH(assembly)
//         } else {
//             // Full BlobTools workflow
//             BLAST_AND_COV(assembly, reads, busco_table, meta_yaml, params.taxon_id)
//             BLOBTOOLS_FILTER(assembly, BLAST_AND_COV.out.blobdir, channel.value([]))
//         }

//     emit:
//         filtered_assembly = params.skip_blobtools
//             ? PASSTHROUGH.out.filtered_assembly
//             : BLOBTOOLS_FILTER.out.filtered_assembly
// }
