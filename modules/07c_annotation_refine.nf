// ============================================================
// 07c_annotation_refine.nf
// Refinement of the merged gene models BEFORE protein extraction, so OMArk,
// StringTie, functional annotation and the NCBI files all see one canonical set.
//   UTRs from BRAKER's stringtie2utr output; transcript span = exon extent
//   gene membership recomputed with a frame-aware rule (undoes AGAT's merges)
//   tRNAs restricted to the EukHighConfidenceFilter set; pseudogenes to a separate file
//   Rfam filters (miRNA, non-eukaryotic rRNA, opposite-strand duplicates, gaps, repeats)
//   evidence / haplotig tags; split of gene models fused across assembly gaps
//   isoform filter: alternative isoforms without evidence of their own are removed
// Scripts: bin/refine_annotation.py, bin/isoform_filter.py, bin/haplotig_scan_classify.py,
// bin/gap_junction_segments.py, bin/classify_gap_junctions.py, bin/protein_identity.py.
// ============================================================

nextflow.enable.dsl = 2

process TRNA_HC_FILTER {
    label 'res_light'
    tag   "trna_hc"

    conda 'bioconda::trnascan-se=2.0.13'
    container 'quay.io/biocontainers/trnascan-se:2.0.13--pl5321hab16a5f_0'
    publishDir "${params.outdir}/07_structural_annotation/refine/trna_hc", mode: params.publish_mode

    input:
    path trna_out
    path trna_tbl

    output:
    path "trna_hc.out",     emit: hc_out
    path "trna_hc_ids.txt", emit: ids
    path "trna_hc.log"
    path "versions.yml",    emit: versions

    script:
    """
    EukHighConfidenceFilter -i ${trna_out} -s ${trna_tbl} -o . -p trna_hc -r
    # IDs as in the tRNAscan-SE GFF3: <seq>.trna<N>
    awk -F'\\t' 'NR>3{gsub(/ /,"",\$1); gsub(/ /,"",\$2); print \$1".trna"\$2}' trna_hc.out > trna_hc_ids.txt
    test -s trna_hc_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trnascan-se: "\$(tRNAscan-SE -h 2>&1 | grep -im1 'tRNAscan-SE [0-9]' | sed 's/tRNAscan-SE //;s/ .*//' || echo unknown)"
    END_VERSIONS
    """

    stub:
    """
    touch trna_hc.out trna_hc_ids.txt trna_hc.log
    printf '"${task.process}":\\n    trnascan-se: stub\\n' > versions.yml
    """
}

process DIAMOND_BLASTP {
    tag "${sname}"

    conda 'bioconda::diamond=2.2.4'
    container 'quay.io/biocontainers/diamond:2.2.4--he361c42_0'

    input:
    tuple val(sname), path(query), path(db_fasta, stageAs: "dbfa/*"), val(args)

    output:
    tuple val(sname), path("${sname}.hits.tsv"), emit: hits
    path "versions.yml", emit: versions

    script:
    """
    diamond makedb --in ${db_fasta} -d db -p ${task.cpus} --quiet
    diamond blastp -q ${query} -d db -o ${sname}.hits.tsv ${args} -p ${task.cpus} --quiet \\
      -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: "\$(diamond version 2>&1 | sed 's/diamond version //')"
    END_VERSIONS
    """

    stub:
    """
    touch ${sname}.hits.tsv
    printf '"${task.process}":\\n    diamond: stub\\n' > versions.yml
    """
}

// Longest protein per BRAKER gene, and whether its scaffold is chromosome-scale
// (scaffold_1..n_chromosome_scaffolds).
process HAPLOTIG_LONGEST {
    label 'res_light'
    tag   "haplotig_longest"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'

    input:
    path braker_results

    output:
    path "longest.faa", emit: faa
    path "genes.tsv",   emit: genes

    script:
    """
    haplotig_scan_classify.py longest ${braker_results}/braker.gtf.gz ${braker_results}/braker.aa.gz \\
      longest.faa genes.tsv ${params.n_chromosome_scaffolds ?: 0}
    """

    stub:
    """
    touch longest.faa genes.tsv
    """
}

process HAPLOTIG_CLASSIFY {
    label 'res_light'
    tag   "haplotig_classify"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/07_structural_annotation/refine/haplotig_scan", mode: params.publish_mode

    input:
    path genes
    tuple val(sname), path(hits)

    output:
    path "haplotig_tags.tsv", emit: tags
    path "haplotig_summary.txt"

    script:
    """
    haplotig_scan_classify.py classify ${genes} ${hits} haplotig_tags.tsv > haplotig_summary.txt
    """

    stub:
    """
    touch haplotig_tags.tsv haplotig_summary.txt
    """
}

// BRAKER transcripts cut into protein segments at introns that contain an assembly
// gap; the related proteomes (one FASTA per species, the file name is the species
// label) are the references for the split test.
process GAP_SEGMENTS {
    label 'res_light'
    tag   "gap_segments"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/07_structural_annotation/refine/gap_fusion_scan", mode: params.publish_mode, pattern: 'junctions.tsv'

    input:
    path braker_results
    path assembly
    path proteomes   // directory of *.faa / *.fa / *.fasta

    output:
    path "segments.faa",  emit: faa
    path "junctions.tsv", emit: junctions
    path "refs.faa",      emit: refs

    script:
    """
    for f in ${proteomes}/*.faa ${proteomes}/*.fa ${proteomes}/*.fasta; do
      [ -e "\$f" ] || continue
      sp=\$(basename "\$f" | sed -E 's/\\.(faa|fa|fasta)\$//')
      awk -v p="\$sp" '/^>/{print ">"p"|"substr(\$1,2); next}{print}' "\$f"
    done > refs.faa
    test -s refs.faa
    gap_junction_segments.py --gtf ${braker_results}/braker.gtf.gz --aa ${braker_results}/braker.aa.gz \\
      --hints ${braker_results}/hintsfile.gff.gz --fasta ${assembly} --out-faa segments.faa --out-tsv junctions.tsv
    """

    stub:
    """
    touch segments.faa junctions.tsv refs.faa
    """
}

process GAP_CLASSIFY {
    label 'res_light'
    tag   "gap_classify"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/07_structural_annotation/refine/gap_fusion_scan", mode: params.publish_mode

    input:
    path junctions
    tuple val(sname), path(hits)

    output:
    path "gap_junction_tags.tsv", emit: tags
    path "classify_summary.txt"

    script:
    """
    classify_gap_junctions.py ${junctions} ${hits} gap_junction_tags.tsv \\
      ${params.gap_split_min_species} ${params.gap_split_min_cov} > classify_summary.txt
    """

    stub:
    """
    touch gap_junction_tags.tsv classify_summary.txt
    """
}

process REFINE_GFF {
    label 'res_mid_4h'
    tag   "refine_gff"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    // with the isoform filter on, its output is the refined.gff3; this one is kept under another name
    publishDir "${params.outdir}/07_structural_annotation/refine", mode: params.publish_mode,
        saveAs: { fn -> fn == 'refined.gff3' && params.related_proteomes && params.isoform_filter.toString().toBoolean()
                        ? 'refined_before_isoform_filter.gff3' : fn }

    input:
    path braker_results
    path validated_gff    // AGAT-merged GFF3 that passed GT_GFF3_VALIDATE
    path assembly
    path trna_ids
    path rm_out
    path haplotig_tags    // or an empty NO_ placeholder
    path gap_tags         // or an empty NO_ placeholder
    path transcript_edits // manual transcript edits, or the NO_TRANSCRIPT_EDITS placeholder

    output:
    path "refined.gff3",          emit: gff
    path "trna_pseudogenes.gff3", emit: pseudo
    path "splits.tsv",            emit: splits
    path "gene_tags.tsv",         emit: tags
    path "refine_report.tsv",     emit: report

    script:
    def opts = [params.refine_keep_mirna.toString().toBoolean() ? '--keep-mirna' : '',
                "--frame-min ${params.refine_frame_min}", "--repeat-max ${params.refine_repeat_max}",
                "--haplotig-pid ${params.refine_haplotig_pid}", "--haplotig-qcov ${params.refine_haplotig_qcov}",
                "--min-intron-reads ${params.isoform_min_reads}"].findAll().join(' ')
    // Nextflow does not hash bin/ scripts: the script's checksum in the command makes -resume rerun this step after an edit
    def script_md5 = file("${projectDir}/bin/refine_annotation.py").text.md5()
    """
    # refine_annotation.py md5 ${script_md5}
    refine_annotation.py --braker-utr-gtf ${braker_results}/braker_utr.gtf.gz --validated-gff ${validated_gff} \\
      --fasta ${assembly} --trna-hc-ids ${trna_ids} --repeatmasker-out ${rm_out} \\
      --gene-support ${braker_results}/gene_support.tsv --haplotig-tags ${haplotig_tags} --gap-tags ${gap_tags} \\
      --out-gff refined.gff3 --out-pseudo trna_pseudogenes.gff3 --out-report refine_report.tsv \\
      --out-splits splits.tsv --out-tags gene_tags.tsv --transcript-edits ${transcript_edits} \\
      --hints ${braker_results}/hintsfile.gff.gz ${opts}
    """

    stub:
    """
    printf '##gff-version 3\\n' > refined.gff3
    touch trna_pseudogenes.gff3 splits.tsv gene_tags.tsv refine_report.tsv
    """
}

// Isoform evidence filter (bin/isoform_filter.py): an alternative isoform stays only with evidence of its
// own (RNA-seq introns and depth, related proteomes); every gene keeps one isoform.
process ISOFORM_PREP {
    label 'res_light'
    tag   "isoform_prep"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'

    input:
    path gff, stageAs: 'in/refined.gff3'
    path assembly

    output:
    path "isoform_proteins.faa", emit: faa
    path "isoform_regions.bed",  emit: bed

    script:
    def script_md5 = file("${projectDir}/bin/isoform_filter.py").text.md5()
    """
    # isoform_filter.py md5 ${script_md5}
    isoform_filter.py prep ${gff} ${assembly} isoform_proteins.faa isoform_regions.bed
    """

    stub:
    """
    touch isoform_proteins.faa isoform_regions.bed
    """
}

// DIAMOND with the aligned sequences (the short-exon test reads the alignment around the exon)
process DIAMOND_ISOFORMS {
    tag "isoforms"

    conda 'bioconda::diamond=2.2.4'
    container 'quay.io/biocontainers/diamond:2.2.4--he361c42_0'

    input:
    path query
    path db_fasta, stageAs: "dbfa/*"

    output:
    path "isoform_hits.tsv", emit: hits
    path "versions.yml",     emit: versions

    script:
    """
    diamond makedb --in ${db_fasta} -d db -p ${task.cpus} --quiet
    diamond blastp -q ${query} -d db -o isoform_hits.tsv --more-sensitive -e 1e-10 -k 50 -p ${task.cpus} --quiet \\
      -f 6 qseqid sseqid pident qstart qend sstart send qlen slen evalue bitscore qseq sseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: "\$(diamond version 2>&1 | sed 's/diamond version //')"
    END_VERSIONS
    """

    stub:
    """
    touch isoform_hits.tsv
    printf '"${task.process}":\\n    diamond: stub\\n' > versions.yml
    """
}

process ISOFORM_BEDCOV {
    label 'res_light'
    tag   "isoform_bedcov"

    conda 'bioconda::samtools=1.21'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

    input:
    path bed
    path bam
    path bai

    output:
    path "isoform_bedcov.tsv", emit: cov

    script:
    """
    samtools bedcov -j ${bed} ${bam} > isoform_bedcov.tsv
    """

    stub:
    """
    touch isoform_bedcov.tsv
    """
}

process ISOFORM_FILTER {
    label 'res_light'
    tag   "isoform_filter"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/07_structural_annotation/refine", mode: params.publish_mode

    input:
    path gff, stageAs: 'in/refined.gff3'
    path assembly
    path braker_results
    path hits
    path bedcov
    path splits, stageAs: 'in/splits.tsv'

    output:
    path "refined.gff3",               emit: gff
    path "isoform_filter_removed.tsv", emit: removed
    path "splits.tsv",                 emit: splits
    path "isoform_filter_summary.txt", emit: summary

    script:
    def script_md5 = file("${projectDir}/bin/isoform_filter.py").text.md5()
    """
    # isoform_filter.py md5 ${script_md5}
    isoform_filter.py filter ${gff} ${assembly} ${braker_results}/hintsfile.gff.gz ${hits} ${bedcov} ${splits} \\
      refined.gff3 isoform_filter_removed.tsv splits.tsv isoform_filter_summary.txt \\
      --min-reads ${params.isoform_min_reads} --cov-diff ${params.isoform_cov_diff} \\
      --start-tol ${params.isoform_start_tol} --short-exon ${params.isoform_short_exon}
    """

    stub:
    """
    cp ${gff} refined.gff3
    cp ${splits} splits.tsv
    touch isoform_filter_removed.tsv isoform_filter_summary.txt
    """
}

// Strict GFF3 validation: any warning or error fails (as NCBI would).
process GT_VALIDATE_REFINED {
    label 'res_light'
    tag   "${gff.baseName}"

    conda 'bioconda::genometools-genometools=1.6.6'
    container 'quay.io/biocontainers/genometools-genometools:1.6.6--py311h21ec246_1'
    publishDir "${params.outdir}/07_structural_annotation/refine", mode: params.publish_mode

    input:
    path gff

    output:
    path "gt_${gff.baseName}.txt", emit: report

    script:
    """
    gt gff3validator ${gff} > gt_${gff.baseName}.txt 2>&1 || true
    if grep -qiE "warning|error" gt_${gff.baseName}.txt; then cat gt_${gff.baseName}.txt >&2; exit 1; fi
    """

    stub:
    """
    echo 'input is valid GFF3' > gt_${gff.baseName}.txt
    """
}

// The refined proteins must equal the unrefined ones, except the split genes.
process PROTEIN_IDENTITY {
    label 'res_light'
    tag   "protein_identity"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/07_structural_annotation/refine", mode: params.publish_mode

    input:
    path refined_proteins, stageAs: 'refined/*'
    path unrefined_proteins, stageAs: 'unrefined/*'
    path splits

    output:
    path "protein_identity.tsv", emit: report

    script:
    """
    protein_identity.py ${refined_proteins} ${unrefined_proteins} ${splits} protein_identity.tsv
    """

    stub:
    """
    printf 'problems\\t0\\n' > protein_identity.tsv
    """
}

// Proteins of the unrefined (AGAT-merged) models, for the protein identity check only; not published.
process GFFREAD_UNREFINED {
    label 'res_light'
    tag   "gffread_unrefined"

    conda 'bioconda::gffread=0.12.7'
    container 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4'

    input:
    path gff3
    path assembly

    output:
    path "unrefined.proteins.faa", emit: proteins

    script:
    """
    gffread ${gff3} -g ${assembly} -y unrefined.proteins.faa
    test -s unrefined.proteins.faa
    """

    stub:
    """
    touch unrefined.proteins.faa
    """
}
