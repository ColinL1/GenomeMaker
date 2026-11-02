// ============================================================
// modules/03_assembly_filtering.nf
// BlobToolKit (v4 CLI) contamination filtering and symbiont split.
// Wired together in subworkflows/03_assembly_filtering.nf.
// ============================================================

nextflow.enable.dsl = 2


// ============================================================
// Process: WRITE_META_YAML
// Build the BlobToolKit meta.yaml describing the assembly.
// ============================================================
def meta_yaml_lines() {
    return ['assembly:',
            "    alias: ${params.species_name}",
            '    record_type: contig'] +
           (params.bioproject ? ["    bioproject: ${params.bioproject}"] : []) +
           (params.biosample ? ["    biosample: ${params.biosample}"] : []) +
           ['taxon:',
            "    name: ${params.species_name}",
            'reads:',
            '    coverage:',
            '        platform: ONT']
}

process WRITE_META_YAML {
    tag   "meta_yaml"

    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03_filtering", mode: params.publish_mode

    output:
    path "meta.yaml", emit: meta_yaml
    path "versions.yml", emit: versions

    script:
    // Without `reads:`, `blobtools filter --summary` crashes on a BlobDir with coverage.
    def meta = meta_yaml_lines().collect { l -> "'${l}'" }.join(' ')
    """
    printf '%s\\n' ${meta} > meta.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    def meta = meta_yaml_lines().collect { l -> "'${l}'" }.join(' ')
    """
    printf '%s\\n' ${meta} > meta.yaml
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: PREPARE_BLAST_QUERY
// Reshape the assembly into a BLASTN query set.
//
// Large scaffolds are cut into overlapping windows so blastn and the
// scheduler have many queries to spread across threads and tasks.
// Windowing makes E-values more permissive (shorter query) and applies
// -max_target_seqs per window, so large-scaffold taxonomy is sampled along
// the whole scaffold.
//
// Do not add -g to `seqkit seq`: there it strips gap characters from the sequence.
// ============================================================
process PREPARE_BLAST_QUERY {
    label 'seqkit_env'
    tag   { assembly.baseName }

    conda 'bioconda::seqkit'
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode_intermediate

    input:
    path assembly

    output:
    path 'blast_query.fasta', emit: query
    path "versions.yml", emit: versions

    script:
    def minlen = params.blast_window_min_len
    """
    seqkit seq -M \$(( ${minlen} - 1 )) ${assembly} > small.fa

    # sliding -g (greedy) keeps the short final window
    seqkit seq -m ${minlen} ${assembly} > big.fa
    seqkit sliding \\
      -s ${params.blast_window_step} \\
      -W ${params.blast_window_size} \\
      -g big.fa > big_win.fa

    cat small.fa big_win.fa > blast_query.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: "\$(seqkit version 2>&1 | sed 's/seqkit v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    # non-empty so splitFasta downstream has a record to chunk
    printf '>stub_contig\\nACGTACGTACGTACGTACGT\\n' > blast_query.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLASTN_NT
// BLASTN assembly vs NCBI nt
//
// -mt_mode 0: on 5 MB windowed chunks, mode 1 was >7x slower and left most
// cores idle. Peak RSS (~224 GB) is reclaimable page cache from the mmapped
// database, so the memory directive is the real working set.
// ============================================================
process BLASTN_NT {
    label 'blast_env'
    tag   { chunk.baseName }

    conda 'bioconda::blast'
    container 'quay.io/biocontainers/blast:2.17.0--h66d330f_0'

    input:
    path chunk
    path blast_db_dir

    output:
    path "${chunk.baseName}.blastn.out", emit: blastn_out
    path "versions.yml", emit: versions

    script:
    """
    blastn -db ${blast_db_dir}/core_nt \\
      -query ${chunk} \\
      -outfmt "6 qseqid staxids bitscore std" \\
      -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 \\
      -mt_mode 0 \\
      -num_threads ${task.cpus} \\
      -out ${chunk.baseName}.raw.out

    # Rebase window coordinates (cols 10/11) onto the parent scaffold and strip
    # the _sliding:<start>-<end> suffix from both qseqid columns (1 and 4).
    awk 'BEGIN{FS=OFS="\\t"}
    {
      if (match(\$1, /_sliding:[0-9]+-[0-9]+\$/)) {
        split(substr(\$1, RSTART+9), c, "-")
        off = c[1] - 1
        \$10 += off; \$11 += off
        sub(/_sliding:[0-9]+-[0-9]+\$/, "", \$1)
        sub(/_sliding:[0-9]+-[0-9]+\$/, "", \$4)
      }
      print
    }' ${chunk.baseName}.raw.out > ${chunk.baseName}.blastn.out

    # blobtools silently drops hits whose qseqid is not in the BlobDir.
    if grep -q '_sliding:' ${chunk.baseName}.blastn.out; then
        echo "ERROR: window IDs survived rebasing in ${chunk.baseName}.blastn.out" >&2
        echo "       blobtools would silently drop these hits -- aborting." >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "\$(blastn -version 2>&1 | head -1 | sed 's/blastn: //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${chunk.baseName}.blastn.out
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: DIAMOND_BLASTX
// DIAMOND blastx vs UniProt over the same windowed chunks as BLASTN_NT.
// ============================================================
process DIAMOND_BLASTX {
    label 'diamond_env'
    tag   { chunk.baseName }

    conda 'bioconda::diamond'
    container 'quay.io/biocontainers/diamond:2.2.4--he361c42_0'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode_intermediate

    input:
    path chunk
    path diamond_db

    output:
    path "${chunk.baseName}.blastx.out", emit: blastx_out
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${params.diamond_tmpdir}
    diamond blastx \\
      --query ${chunk} \\
      --db ${diamond_db} \\
      --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
      -k 1 --evalue 1e-25 \\
      --block-size ${params.diamond_block_size} \\
      --index-chunks ${params.diamond_index_chunks} \\
      --threads ${task.cpus} \\
      --tmpdir ${params.diamond_tmpdir} \\
      --verbose \\
      > ${chunk.baseName}.raw.out

    # Same rebasing as BLASTN_NT (blastx may report qstart > qend; the offset still applies).
    awk 'BEGIN{FS=OFS="\\t"}
    {
      if (match(\$1, /_sliding:[0-9]+-[0-9]+\$/)) {
        split(substr(\$1, RSTART+9), c, "-")
        off = c[1] - 1
        \$10 += off; \$11 += off
        sub(/_sliding:[0-9]+-[0-9]+\$/, "", \$1)
        sub(/_sliding:[0-9]+-[0-9]+\$/, "", \$4)
      }
      print
    }' ${chunk.baseName}.raw.out > ${chunk.baseName}.blastx.out

    # blobtools silently drops hits whose qseqid is not in the BlobDir.
    if grep -q '_sliding:' ${chunk.baseName}.blastx.out; then
        echo "ERROR: window IDs survived rebasing in ${chunk.baseName}.blastx.out" >&2
        echo "       blobtools would silently drop these hits -- aborting." >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: "\$(diamond --version 2>&1 | sed 's/.*version //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${chunk.baseName}.blastx.out
    echo "DIAMOND_BLASTX stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MINIMAP2_SORT
// ONT alignment piped straight into a sorted, indexed BAM.
// ============================================================
process MINIMAP2_SORT {
    label 'minimap2_env'
    tag   "minimap2_sort"

    conda 'bioconda::minimap2 bioconda::samtools'
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode_intermediate

    input:
    path assembly
    path reads

    output:
    path "coverage.bam",     emit: cov_bam
    path "coverage.bam.bai", emit: cov_bai
    path "coverage.flagstat", emit: flagstat
    path "versions.yml", emit: versions

    script:
    // samtools sort -m is per thread; leave the rest for minimap2's index.
    def sort_threads = Math.max(1, (task.cpus / 4) as int)
    def sort_mem_mb  = Math.max(256, (task.memory.toMega() * 0.4 / sort_threads) as long)
    """
    ( set -o pipefail
      minimap2 -a -x map-ont -t ${task.cpus} ${assembly} ${reads} \\
        | samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -O BAM -o coverage.bam - )
    samtools index -@ ${sort_threads} coverage.bam

    # an empty coverage track would be added silently
    samtools flagstat -@ ${sort_threads} coverage.bam > coverage.flagstat
    MAPPED=\$(awk '/ mapped \\(/{print \$1; exit}' coverage.flagstat)
    if [ "\${MAPPED:-0}" -eq 0 ]; then
        echo "ERROR: no reads mapped to ${assembly} -- coverage track would be empty" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version 2>&1 | head -1 || echo unknown)"
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch coverage.bam coverage.bam.bai
    printf '100 + 0 in total (QC-passed reads + QC-failed reads)\\n95 + 0 mapped (95.00%% : N/A)\\n' > coverage.flagstat
    echo "MINIMAP2_SORT stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "stub"
        samtools: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_CREATE
// Initialise the BlobDir
// ============================================================
process BLOBTOOLS_CREATE {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_create"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    // intermediate BlobDir, not published

    input:
    path assembly
    path meta_yaml
    val  taxid   // NCBI taxid of the host (anchors lineage in BlobDir, not a filter)
    path taxdump_dir

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    blobtools create \\
      --fasta ${assembly} \\
      --meta ${meta_yaml} \\
      --taxid ${taxid} \\
      --taxdump ${taxdump_dir} \\
      btk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p btk
    echo "BLOBTOOLS_CREATE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_ADD_HITS
// Add BLASTN + DIAMOND taxonomy hits to a copy of the BlobDir
// ============================================================
process BLOBTOOLS_ADD_HITS {
    label 'blobtools_env'
    label 'res_small_2h'
    tag   "blobtools_add_hits"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    // intermediate BlobDir, not published

    input:
    path btk_in, stageAs: 'btk_in'   // copied: `blobtools add` edits the BlobDir in place
    path blastn_out
    path blastx_out
    path taxdump_dir

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    # --taxdump is required by `add`, not only `create`.
    # Hit order matters: for bestsumorder the first file with a hit wins, so
    # nt comes first and DIAMOND is the fallback.
    cp -rL btk_in btk
    blobtools add \\
      --hits ${blastn_out} \\
      --hits ${blastx_out} \\
      --taxrule bestsumorder \\
      --hit-count ${params.blobtools_hit_count} \\
      --taxdump ${taxdump_dir} \\
      btk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    cp -rL btk_in btk
    echo "BLOBTOOLS_ADD_HITS stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_ADD_COV
// Add ONT coverage to a copy of the BlobDir
// ============================================================
process BLOBTOOLS_ADD_COV {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_add_cov"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    // intermediate BlobDir, not published

    input:
    path btk_in, stageAs: 'btk_in'   // copied: see BLOBTOOLS_ADD_HITS
    path cov_bam
    // no .bai: blobtools 4.4.4 crashes deleting its own .csi when a .bai is present

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    cp -rL btk_in btk
    blobtools add --cov ${cov_bam} btk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    cp -rL btk_in btk
    echo "BLOBTOOLS_ADD_COV stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_ADD_BUSCO
// Add BUSCO evidence to a copy of the BlobDir
// ============================================================
process BLOBTOOLS_ADD_BUSCO {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_add_busco"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode

    input:
    path btk_in, stageAs: 'btk_in'   // copied: see BLOBTOOLS_ADD_HITS
    path busco_table

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    cp -rL btk_in btk
    blobtools add --busco ${busco_table} btk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    cp -rL btk_in btk
    echo "BLOBTOOLS_ADD_BUSCO stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_APPLY_FILTER
// Filter one category (host/symbiont) by a manual Viewer JSON
// (params.filter_json) or the automatic best-hit taxon filter.
// ============================================================
process BLOBTOOLS_APPLY_FILTER {
    label 'blobtools_env'
    label 'res_mid_4h'
    tag   "apply_filter_${category}"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    publishDir "${params.outdir}/03_filtering/filtered", mode: params.publish_mode

    input:
    tuple val(category), val(taxon_rank), val(taxon_name), path(filter_json), path(assembly), path(blobdir)

    output:
    tuple val(category), path("${category}_assembly.fasta"), emit: filtered_assembly
    tuple val(category), path("${category}_filtered"), emit: filtered_blobdir
    tuple val(category), path("${category}_summary.json"), emit: summary
    path "${category}_scaffolds.tsv", emit: scaffolds
    path "versions.yml", emit: versions

    script:
    // --param Keys=A,B excludes A,B; --invert keeps them instead.
    def selection = filter_json ? "--json ${filter_json}" : "--param ${taxon_rank}--Keys=\$KEYS --invert"
    def keys_cmd = filter_json ? "" : """
    # host also keeps "no-hit", only if present (blobtools errors on an absent key)
    KEYS="${taxon_name}"
    ${category == 'host' ? """blobtools filter --table ${category}_ranks.tsv --table-fields ${taxon_rank} ${blobdir}
    # blobtools tables use CRLF
    if awk -F'\\t' '{sub(/\\r\$/, "")} NR>1 && \$3=="no-hit" {f=1} END {exit !f}' ${category}_ranks.tsv; then KEYS="\$KEYS,no-hit"; fi""" : ""}
    echo "${category}: keeping ${taxon_rank} = \$KEYS"
    """
    """
    ${keys_cmd}
    # the filtered FASTA is written next to the input, not into --output
    blobtools filter \\
      ${selection} \\
      --fasta ${assembly} \\
      --output ${category}_filtered \\
      --suffix ${category} \\
      --summary ${category}_summary.json \\
      --table ${category}_scaffolds.tsv \\
      --table-fields length \\
      ${blobdir}

    N=\$(( \$(wc -l < ${category}_scaffolds.tsv) - 1 ))
    if [ "\$N" -le 0 ]; then
      echo "ERROR: selection for category '${category}' matched 0 scaffolds (${filter_json ? 'json ' + filter_json : taxon_rank + '=' + taxon_name})" >&2
      exit 1
    fi
    echo "${category}: \$N scaffolds selected"

    OUT=\$(find . -maxdepth 1 -name "*.${category}.fa*" | head -1)
    [ -n "\$OUT" ] && cp "\$OUT" ${category}_assembly.fasta

    if [ ! -s ${category}_assembly.fasta ]; then
      echo "ERROR: no filtered fasta produced for category '${category}'" >&2
      exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${category}_assembly.fasta ${category}_summary.json ${category}_scaffolds.tsv
    mkdir -p ${category}_filtered
    echo "BLOBTOOLS_APPLY_FILTER stub completed (${category})"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_STATIC_PLOTS
// Static blob/snail/cumulative PNGs of the stage's final BlobDir.
// ============================================================
process BLOBTOOLS_STATIC_PLOTS {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_plots"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    publishDir "${params.outdir}/03_filtering/plots", mode: params.publish_mode

    input:
    path blobdir

    output:
    path "*.png", emit: plots
    path "versions.yml", emit: versions

    script:
    """
    # blobtools view needs a writable HOME
    export HOME=\$PWD

    for VIEW in blob snail cumulative; do
      blobtools view \\
        --plot \\
        --view \$VIEW \\
        --format png \\
        --out . \\
        ${blobdir}
    done

    for VIEW in blob snail cumulative; do
      ls *\${VIEW}*.png >/dev/null 2>&1 || { echo "ERROR: no \$VIEW plot produced" >&2; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch stub.blob.circle.png stub.snail.png stub.cumulative.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: PASSTHROUGH
// Used when params.run_blobtools = false
// ============================================================
process PASSTHROUGH {
    tag   "passthrough"

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/03_filtering/filtered", mode: params.publish_mode

    input:
    path assembly

    output:
    path "filtered_assembly.fasta", emit: filtered_assembly
    path "versions.yml", emit: versions

    script:
    """
    cp ${assembly} filtered_assembly.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch filtered_assembly.fasta
    echo "PASSTHROUGH stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}
