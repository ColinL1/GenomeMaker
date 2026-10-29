// ============================================================
// 10_comparison.nf
// Processes for the Porites comparative analysis (porites_compare.nf):
// reduce an annotated NCBI genome to bare gene models for re-annotation with
// stage 08, and run `funannotate compare`. The NCBI download itself is done by
// comparison/run_compare.sh with the site's `ncbi_datasets` conda env (the
// biocontainers ncbi-datasets-cli image is too old for the current NCBI API).
//
// Resources are declared here rather than in process_resource.config: these
// processes are not part of main.nf and the selectors there would be unused.
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PREP_REF_ANNOTATION
// NCBI FASTA/GFF3 -> the inputs stage 08 expects: sequence names of at most
// 16 characters and a GFF3 holding only gene models (see bin/ncbi_ref_prep.pl).
// Isoforms are kept as delivered, as they are for the pipeline's own annotation.
// ============================================================
process PREP_REF_ANNOTATION {
    tag   "${params.species_name}"
    cpus   4
    memory 32.GB
    time   6.h

    conda 'bioconda::agat=1.7.0'
    container 'quay.io/biocontainers/agat:1.7.0--pl5321hdfd78af_0'
    publishDir "${params.outdir}/01_prepared", mode: params.publish_mode

    input:
    path fasta
    path gff

    output:
    path "${params.species_name}.genome.fa",    emit: assembly
    path "${params.species_name}.models.gff3",  emit: gff3
    path "${params.species_name}.seqmap.tsv",   emit: seqmap
    path "${params.species_name}.prep_summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    def sp = params.species_name
    """
    ncbi_ref_prep.pl fasta ${fasta} ${sp}.genome.fa ${sp}.seqmap.tsv

    # Normalise first (complete gene/mRNA/exon/CDS hierarchy, one Parent per
    # line), then cut down to the models to re-annotate.
    agat_convert_sp_gxf2gxf.pl -g ${gff} -o normalised.gff3
    ncbi_ref_prep.pl gff normalised.gff3 ${sp}.seqmap.tsv \\
      > ${sp}.models.gff3 2> ${sp}.prep_summary.txt
    cat ${sp}.prep_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "\$(agat_convert_sp_gxf2gxf.pl --version 2>&1 | grep -oim1 'v[0-9.]*' | sed 's/v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    def sp = params.species_name
    """
    touch ${sp}.genome.fa ${sp}.models.gff3 ${sp}.seqmap.tsv ${sp}.prep_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: DROP_UNTRANSLATABLE
// Remove mRNAs whose CDS does not translate cleanly from the genome (in-frame
// stop). RefSeq has a few hundred: models whose protein NCBI corrected relative
// to the assembly, and selenoproteins. gffread writes those stops as '.', which
// Phobius rejects outright, and funannotate would translate them the same way.
// Genes left without any transcript are removed with them.
// ============================================================
process DROP_UNTRANSLATABLE {
    tag   "${params.species_name}"
    cpus   2
    memory 8.GB
    time   2.h

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'
    publishDir "${params.outdir}/01_prepared", mode: params.publish_mode

    input:
    path gff3
    path assembly

    output:
    path "${params.species_name}.gff3",                emit: gff3
    path "${params.species_name}.untranslatable.txt",  emit: dropped
    path "versions.yml", emit: versions

    script:
    def sp = params.species_name
    """
    gffread ${gff3} -g ${assembly} -y all_models.faa
    awk '/^>/{id=substr(\$1,2); next} /\\./{print id}' all_models.faa | sort -u > ${sp}.untranslatable.txt

    # ${gff3} carries only "ID=..;Parent=..[;product=..]" (bin/ncbi_ref_prep.pl).
    # pass 1: genes that keep a transcript; pass 2: print what survives.
    awk -F'\\t' '
      function attr(key,   n, a, i) {
          n = split(\$9, a, ";")
          for (i = 1; i <= n; i++) if (index(a[i], key "=") == 1) return substr(a[i], length(key) + 2)
          return ""
      }
      pass == 0 { bad[\$1] = 1; next }
      /^#/      { if (pass == 2) print; next }
      pass == 1 { if (\$3 != "gene" && \$3 != "exon" && \$3 != "CDS" && !(attr("ID") in bad)) keep[attr("Parent")] = 1; next }
      \$3 == "gene"                 { if (attr("ID") in keep) print; next }
      \$3 == "exon" || \$3 == "CDS" { if (!(attr("Parent") in bad)) print; next }
                                    { if (!(attr("ID") in bad)) print }
    ' pass=0 ${sp}.untranslatable.txt pass=1 ${gff3} pass=2 ${gff3} > ${sp}.gff3

    echo "mRNAs dropped (in-frame stop): \$(wc -l < ${sp}.untranslatable.txt)"
    echo "genes kept: \$(awk -F'\\t' '\$3=="gene"' ${sp}.gff3 | wc -l)"
    echo "mRNAs kept: \$(awk -F'\\t' '\$3=="mRNA"' ${sp}.gff3 | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "\$(gffread --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    def sp = params.species_name
    """
    touch ${sp}.gff3 ${sp}.untranslatable.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: FUNANNOTATE_GENE_IDS
// Reshape gene IDs for `funannotate annotate`. Before merging its annotation
// sources it keeps only lines that start with a "base" it derives from each
// protein header "<transcript> <gene>" (library.getGeneBasename): with a "_" in
// the gene ID the base is everything up to the first "_" of the whole header,
// which spans the space and matches nothing -- every Pfam/InterPro/eggNOG/GO
// annotation is then dropped silently ("adding 0 valid annotations"). NCBI
// locus tags (gene-ABFA07_000001) hit exactly that. With "-" and no "_" in the
// gene ID the base is the transcript's own prefix (rna), which always matches.
// Transcript IDs are left alone, so the proteins and every search run on them
// stay valid.
// ============================================================
process FUNANNOTATE_GENE_IDS {
    tag   "${params.species_name}"
    cpus   1
    memory 4.GB
    time   1.h

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'
    publishDir "${params.outdir}/01_prepared", mode: params.publish_mode

    input:
    path gff3

    output:
    path "${params.species_name}.funannotate.gff3", emit: gff3

    script:
    """
    awk -F'\\t' -v OFS='\\t' '
      function fix(id) { gsub(/_/, ".", id); if (id !~ /-/) id = "gene-" id; return id }
      /^#/ { print; next }
      \$3 == "gene" { \$9 = "ID=" fix(substr(\$9, 4)); print; next }
      \$3 == "mRNA" || \$3 == "tRNA" || \$3 == "rRNA" {
          n = split(\$9, a, ";")
          if (a[1] !~ /^ID=[^;]*-/) { print "transcript ID without a dash: " a[1] > "/dev/stderr"; exit 1 }
          a[2] = "Parent=" fix(substr(a[2], 8))
          \$9 = a[1]; for (i = 2; i <= n; i++) \$9 = \$9 ";" a[i]
          print; next
      }
      { print }
    ' ${gff3} > ${params.species_name}.funannotate.gff3
    """

    stub:
    """
    touch ${params.species_name}.funannotate.gff3
    """
}

// ============================================================
// Process: FUNANNOTATE_COMPARE
// Comparative genomics across funannotate-annotated genomes (GenBank files
// written by `funannotate annotate`).
// ============================================================
process FUNANNOTATE_COMPARE {
    tag   "funannotate_compare"
    // proteinortho's all-vs-all search is the parallel part; the run has the node
    // to itself by default (params.compare_cpus)
    // (params set in the entry script do not reach an included module, hence the defaults here)
    cpus   { (params.compare_cpus ?: 240) as int }
    memory { 256.GB * task.attempt }
    time   { 96.h * task.attempt }

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path gbks
    path funannotate_db   // configured funannotate database (funannotate setup)
    path proteinortho_tsv // precomputed ProteinOrtho table, or a NO_ placeholder

    output:
    path "funannotate_compare", emit: compare_dir
    path "versions.yml",        emit: versions

    script:
    // By default compare runs `proteinortho -synteny` itself. On fragmented
    // assemblies that synteny (FFAdj-MCS) step does not finish: after 2.5 h only the
    // 3 pairs involving the chromosome-level P. lutea were done, one job had found
    // 190 million candidate runs and the 12 jobs held 340 GB. A table made without
    // -synteny (comparison/proteinortho_nosynteny.sh) can be handed in instead.
    def po_arg = proteinortho_tsv.name.startsWith('NO_') ? '' : "--proteinortho ${proteinortho_tsv}"
    // --ml_method raxml: the image ships RAxML but not IQ-TREE (compare's
    // default). No --run_dnds: codeml is not in the image either. No --outgroup:
    // FUNANNOTATE_DB/outgroups only holds fungal BUSCO sets.
    """
    export FUNANNOTATE_DB=\$(readlink -f ${funannotate_db})

    funannotate compare \\
      -i ${gbks} \\
      -o funannotate_compare \\
      --ml_method raxml \\
      ${po_arg} \\
      --num_orthos ${params.compare_num_orthos ?: 500} \\
      --cpus ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "\$(funannotate --version 2>&1 | grep -o '[0-9][0-9.]*\$' | tail -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p funannotate_compare
    echo "FUNANNOTATE_COMPARE stub completed" > funannotate_compare/index.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "stub"
    END_VERSIONS
    """
}
