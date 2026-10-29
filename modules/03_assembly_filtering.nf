// ============================================================
// modules/03_assembly_filtering.nf
// BlobToolKit-based assembly filtering -- PROCESS DEFINITIONS ONLY
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
//
// MODULARITY NOTE: each tool below is its own process:
//   - BLASTN_NT / DIAMOND_BLASTX: taxonomy hits
//   - MINIMAP2_ALIGN / SAMTOOLS_SORT_INDEX: coverage, split so minimap2
//     and samtools each get their own container/conda env
//   - BLOBTOOLS_CREATE / _ADD_HITS / _ADD_COV / _ADD_BUSCO: BlobDir build,
//     chained by passing the BlobDir folder through as both input and
//     output of each step (`blobtools add` modifies it in place)
//   - BLOBTOOLS_APPLY_FILTER: per-category filter + FASTA (host/symbiont),
//     by automatic taxon filter or an optional manually exported JSON
//   - PASSTHROUGH: no-op used when params.skip_blobtools = true
//
// This file has NO workflow{} block -- see
// subworkflows/03_assembly_filtering.nf for the ASSEMBLY_FILTERING
// workflow that wires these processes together.
// ============================================================

nextflow.enable.dsl = 2


// ============================================================
// Process: WRITE_META_YAML
// Build the BlobToolKit meta.yaml describing the assembly.
//
// WHY A PROCESS: main.nf used to write this file into ${workDir} directly
// during DAG construction:
//     def meta_yaml_file = file("${workDir}/meta_yaml/meta.yaml")
//     meta_yaml_file.parent.mkdirs()
//     meta_yaml_file.text = meta_yaml_content
// That is a side effect at parse time rather than a task, so it is outside
// Nextflow's provenance and caching, it pollutes the work dir root, and on a
// resumed run it silently rewrites a file that completed tasks already consumed.
// As a process it is cached, provenanced and published like everything else.
// ============================================================
process WRITE_META_YAML {
    tag   "meta_yaml"

    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03_filtering", mode: params.publish_mode

    output:
    path "meta.yaml", emit: meta_yaml
    path "versions.yml", emit: versions

    script:
    """
    # `reads:` declares the coverage library (named after the BAM, coverage.bam,
    # by MINIMAP2_SORT). Without it `blobtools filter --summary` crashes in
    # cov.py (AttributeError: reads) on any BlobDir that has coverage added.
    cat > meta.yaml <<'META'
assembly:
    alias: ${params.species_name}
    bioproject: ${params.bioproject}
    biosample: ${params.biosample}
    record_type: contig
taxon:
    name: ${params.species_name}
reads:
    coverage:
        platform: ONT
META
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf 'assembly:\\n    alias: ${params.species_name}\\ntaxon:\\n    name: ${params.species_name}\\nreads:\\n    coverage:\\n        platform: ONT\\n' > meta.yaml
    
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
// WHY: the purged assembly is ~14 chromosome-scale scaffolds (18-47 Mb,
// 75% of the sequence mass) plus ~16.7k small contigs. Each big scaffold
// is a SINGLE fasta record, which is bad twice over:
//   1. Nextflow's unit of work is the chunk, so one 47 Mb record pins a
//      task until it finishes -- the step's wall time is floored by the
//      slowest scaffold no matter how many cores sit idle.
//   2. blastn -mt_mode 1 parallelises by splitting QUERIES across threads.
//      With nseq=1 there is nothing to split, so it silently falls back to
//      -mt_mode 0 (database split), which flattens out past ~8-16 threads.
//
// Windowing the big scaffolds turns 14 indivisible records into ~4.5k
// independent queries, which both -mt_mode 1 and the scheduler can spread.
// Contigs below params.blast_window_min_len are passed through untouched,
// so the small-contig fraction keeps its original IDs and statistics.
//
// Windows OVERLAP by (window_size - window_step) = 10 kb by default. Any
// HSP shorter than the overlap is guaranteed to fall entirely inside at
// least one window, which is what stops boundary-straddling alignments
// from being truncated into two sub-threshold pieces.
//
// ACCURACY NOTE -- this is not a free optimisation, it changes results:
//   * BLAST E-values scale with query length (E ~ K*m*n*e^-lambda*S), so a
//     110 kb window makes the same alignment ~400x more significant than it
//     was against a 47 Mb query. Against a fixed -evalue the search becomes
//     MORE permissive. BlobTools ranks by bitscore, which is length
//     independent, so the taxonomic ordering itself is unaffected.
//   * -max_target_seqs is applied PER QUERY. A 47 Mb scaffold previously
//     got 10 hits total to represent 47 Mb; it now gets up to 10 per
//     window. Under the bestsumorder taxrule this means big-scaffold
//     taxonomy is decided on evidence sampled across the whole scaffold
//     rather than 10 arbitrary hits. Expect different (better-sampled)
//     calls for the large scaffolds -- worth stating in methods.
//
// Verified on this assembly: 16,710 records -> 16,696 passthrough + 14
// scaffolds -> 4,464 windows; bases 146,561,500 + 445,753,291 =
// 592,314,791 = original exactly; window ranges tile each scaffold 1..L
// with no gaps; stripping the _sliding: suffix reproduces all 16,710
// original IDs exactly.
//
// NOTE: do NOT add -g/--remove-gaps to `seqkit seq` here. In `seqkit seq`
// -g means remove-gaps and strips "- \t." from the SEQUENCE (it is only
// "greedy" in `seqkit sliding`). It would silently alter the contigs.
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
    # contigs below the threshold: passed through with original IDs
    seqkit seq -M \$(( ${minlen} - 1 )) ${assembly} > small.fa

    # scaffolds at/above the threshold: cut into overlapping windows.
    # -g here is sliding's GREEDY flag, so the final short window at the
    # end of each scaffold is emitted rather than dropped.
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
    # emit a minimal valid record so splitFasta downstream has something to
    # chunk -- an empty file would make every downstream stub vacuous
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
// -mt_mode 0 (database split), NOT -mt_mode 1 (query split). MEASURED.
//
// This pipeline previously used -mt_mode 1, on the reasoning that windowing
// "guarantees" many query records for blastn to split across threads. That
// reasoning does not survive contact with the actual chunks: a 5 MB chunk of the
// WINDOWED region holds only ~48 records of ~110 kb, i.e. 1.5 records per thread
// at 32 threads, so there is essentially nothing to split.
//
// Measured on real chunks against core_nt, same 5 MB windowed chunk, 32 threads:
//     -mt_mode 0 : 17:59, 224 GB peak RSS, 480 hits
//     -mt_mode 1 : still running after 2h12m when the test was abandoned,
//                  sustaining only 222-335% CPU of the 3200% requested
// i.e. mode 1 is at least 7.4x slower here and leaves ~90% of its cores idle.
//
// Passthrough chunks (~580 records, 18 per thread) would suit mode 1 far better,
// but they are only 25% of the sequence mass, so one setting has to serve both
// and mode 0 is the one that does not collapse on the expensive 75%.
//
// MEMORY: mode 0 mmaps database partitions, so peak RSS reads ~224 GB -- but
// smaps_rollup shows Anonymous = 0 and Private_Dirty = 0, i.e. ALL of it is clean
// file-backed page cache the kernel reclaims on demand (observed dropping
// 183 -> 87 GB mid-run). The memory directive therefore declares the real working
// set, not the mapped database. Do not "fix" it upward from a peak-RSS number.
//
// NOT MEASURED: whether a chunk size other than 5 MB suits mode 0 better. Mode 0
// pays a full database pass per invocation, so larger chunks would likely amortise
// better -- but that was not tested, so 5 MB stands.
//
// Also note the database changed from full nt to core_nt (264 GB / 80 volumes,
// vs 723 GB / 246), so the chunk-size tuning in nextflow.config predates the
// database currently in use.
//
// Resources are flat rather than branching on chunk.size(): after
// windowing, chunks are uniform (params.blastn_chunk_size of many small
// records), so the old "is this one giant scaffold?" branch no longer has
// anything to detect.
//
// The trailing sed restores the parent scaffold ID by stripping the
// _sliding:<start>-<end> suffix that PREPARE_BLAST_QUERY added, so
// downstream BlobTools sees the original assembly IDs.
//
// The strip is GLOBAL, not anchored to field 1, because `std` expands to
// "qseqid sseqid pident length mismatch gapopen qstart qend sstart send
// evalue bitscore" -- so qseqid appears in BOTH column 1 and column 4 and
// the two must stay consistent. `_sliding:<int>-<int>` is a suffix this
// pipeline creates itself, so a global strip cannot collide with an nt
// subject accession.
//
// CAVEAT: qstart/qend (columns 10/11) remain WINDOW-relative -- they are
// not rebased onto the parent scaffold. BlobTools only reads qseqid,
// staxids and bitscore, so this is harmless here, but do not reuse the
// published blastn.out for coordinate-based analysis without adding the
// window offset back.
// ============================================================
process BLASTN_NT {
    label 'blast_env'
    tag   { chunk.baseName }

    conda 'bioconda::blast'
    container 'quay.io/biocontainers/blast:2.17.0--h66d330f_0'  // verify current tag against your install
    // Per-chunk output is NOT published: there are ~120 of them and the
    // subworkflow collectFile()s them into one blastn.out. Publishing both
    // duplicates the whole search result as directory noise.

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

    # Rebase window-local query coordinates onto the parent scaffold, THEN
    # collapse the window ID. Both must happen together: once the suffix is
    # gone the offset is unrecoverable.
    #
    # Columns for -outfmt "6 qseqid staxids bitscore std" are
    #   1 qseqid  2 staxids  3 bitscore  4 qseqid  5 sseqid  6 pident
    #   7 length  8 mismatch 9 gapopen  10 qstart 11 qend   12 sstart
    #  13 send   14 evalue                       (14 columns, evalue last)
    # which is exactly what blobtools' default --hits-cols expects.
    # qstart/qend always ascend for blastn (strand is carried by
    # sstart/send), so a flat offset is correct for both strands.
    # "_sliding:" is 9 characters, hence RSTART+9.
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

    # GUARD -- do not remove. blobtools SILENTLY drops hits whose qseqid is
    # absent from the BlobDir: exit 0, no warning, the sequence just becomes
    # "no-hit". The windowed scaffolds are ~75% of this assembly, and
    # BLOBTOOLS_APPLY_FILTER keeps only the host taxon, so a rename that
    # failed to fire would quietly delete most of the genome with a clean
    # exit code. Verified against blobtoolkit 4.4.4. Fail loud instead.
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
// DIAMOND blastx vs UniProt reference proteomes, over CHUNKS of the same
// windowed query BLASTN_NT uses.
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

    # Rebase window-local query coordinates onto the parent scaffold, THEN
    # collapse the window ID (both qseqid columns). Same layout as BLASTN_NT's
    # outfmt: 1 qseqid, 2 staxids, 3 bitscore, 4 qseqid, 5 sseqid, ...,
    # 10 qstart, 11 qend, 14 evalue -- blobtools' default --hits-cols, and it
    # needs qstart/qend to bin hits. Unlike blastn, blastx reports qstart > qend
    # on reverse-strand frames; a flat offset is still correct because both
    # coordinates are window-local and shift together.
    # "_sliding:" is 9 characters, hence RSTART+9.
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

    # GUARD -- do not remove. blobtools SILENTLY drops hits whose qseqid is absent
    # from the BlobDir: exit 0, no warning, the sequence just becomes "no-hit".
    # The windowed scaffolds are ~75% of this assembly, so a rename that failed to
    # fire would quietly mis-call most of the genome with a clean exit code.
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
// ONT long-read alignment (map-ont preset) piped STRAIGHT into a coordinate
// sorted, indexed BAM.
//
// This replaces the former MINIMAP2_ALIGN + SAMTOOLS_SORT_INDEX pair. That split
// was documented as deliberate ("so minimap2 and samtools each get their own
// container"), but the cost was severe: minimap2 wrote an uncompressed SAM of a
// ~600 Mb assembly against 14.5 GB of ONT reads to the work dir, and a second
// task then read the whole thing back -- likely several hundred GB of
// write-then-read on the rotational array holding work/.
//
// Piping through samtools sort means the SAM never exists on disk. The container
// is the same mulled minimap2+samtools image MITO_MAPPING already uses, so no new
// dependency is introduced.
//
// BAM, never SAM -- here and everywhere else in this pipeline.
// ============================================================
process MINIMAP2_SORT {
    label 'minimap2_env'
    tag   "minimap2_sort"

    conda 'bioconda::minimap2 bioconda::samtools'
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode_intermediate

    input:
    path assembly
    path reads   // ONT long reads (porechopped)

    output:
    path "coverage.bam",     emit: cov_bam
    path "coverage.bam.bai", emit: cov_bai
    path "coverage.flagstat", emit: flagstat
    path "versions.yml", emit: versions

    script:
    // Sort threads are capped and given an explicit per-thread memory budget:
    // samtools sort -m is PER THREAD, so -@ N -m M must satisfy N*M < the
    // allocation or the task is OOM-killed. Half the allocation is left for
    // minimap2's index, which for a 600 Mb genome is the larger consumer.
    def sort_threads = Math.max(1, (task.cpus / 4) as int)
    def sort_mem_mb  = Math.max(256, (task.memory.toMega() * 0.4 / sort_threads) as long)
    """
    ( set -o pipefail
      minimap2 -a -x map-ont -t ${task.cpus} ${assembly} ${reads} \\
        | samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -O BAM -o coverage.bam - )
    samtools index -@ ${sort_threads} coverage.bam

    # Cheap sanity check plus a MultiQC-parseable alignment summary: a coverage BAM
    # with no mapped reads would make BLOBTOOLS_ADD_COV silently add an empty track.
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
// Initialize the BlobDir (fasta + metadata + host taxid anchor)
// ============================================================
process BLOBTOOLS_CREATE {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_create"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    // NOT published: this is an intermediate BlobDir state. All four
    // blobtools steps used to publish `btk` to the SAME directory, so the
    // 4 GB+ BlobDir was copied four times, each overwriting the last.
    // Only BLOBTOOLS_ADD_BUSCO -- the final state -- publishes now.

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
// Add BLASTN + DIAMOND taxonomy hits to the BlobDir (in place)
// ============================================================
process BLOBTOOLS_ADD_HITS {
    label 'blobtools_env'
    label 'res_small_2h'
    tag   "blobtools_add_hits"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    // NOT published: this is an intermediate BlobDir state. All four
    // blobtools steps used to publish `btk` to the SAME directory, so the
    // 4 GB+ BlobDir was copied four times, each overwriting the last.
    // Only BLOBTOOLS_ADD_BUSCO -- the final state -- publishes now.

    input:
    path btk
    path blastn_out
    path blastx_out
    path taxdump_dir

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    # --taxdump is REQUIRED here, not optional: without it blobtools dies in
    # taxdump.py with "TypeError: expected str, bytes or os.PathLike object,
    # not NoneType" when it tries to join None with nodes.dmp. `create`
    # taking --taxdump does not cover `add`.
    #
    # ORDER MATTERS: for the *order taxrules (bestsumorder), the first
    # --hits file that has any hit for a sequence wins, and later files are
    # only consulted as a fallback. blastn/nt is listed FIRST so nucleotide
    # hits take precedence and the DIAMOND/UniProt blastx hits act as the
    # fallback for sequences nt could not place.
    # --hit-count caps how many hit rows (top N by bitscore) feed the
    # per-taxon bitscore sum. The default of 10 was sized for ONE blast
    # call per sequence; with windowing a single 110 kb region can supply
    # all 10 slots and decide a 47 Mb scaffold. Raising it lets the whole
    # scaffold vote. Measured effect of this knob on a synthetic conflict
    # (20 hits @ bitscore 500 vs 5 @ 1500): --hit-count 10 -> Cnidaria,
    # --hit-count 25 -> Pseudomonadota. It changes filtering outcomes, so
    # it is a param rather than a hard-coded value.
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
    echo "BLOBTOOLS_ADD_HITS stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_ADD_COV
// Add ONT coverage to the BlobDir (in place)
// ============================================================
process BLOBTOOLS_ADD_COV {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_add_cov"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    // NOT published: this is an intermediate BlobDir state. All four
    // blobtools steps used to publish `btk` to the SAME directory, so the
    // 4 GB+ BlobDir was copied four times, each overwriting the last.
    // Only BLOBTOOLS_ADD_BUSCO -- the final state -- publishes now.

    input:
    path btk
    path cov_bam
    // The .bai is deliberately NOT staged: blobtools 4.4.4 builds its own
    // .csi index and deletes it afterwards, but if a .bai is present it
    // builds no .csi and then crashes on os.remove('<bam>.csi').

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    blobtools add --cov ${cov_bam} btk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    echo "BLOBTOOLS_ADD_COV stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_ADD_BUSCO
// Add BUSCO evidence to the BlobDir (in place) -- host-specific
// lineage (e.g. anthozoa_odb12), flags host contigs unambiguously
// ============================================================
process BLOBTOOLS_ADD_BUSCO {
    label 'blobtools_env'
    label 'res_small'
    tag   "blobtools_add_busco"

    conda 'bioconda::blobtoolkit'
    container 'genomehubs/blobtoolkit:4.4.4'
    publishDir "${params.outdir}/03_filtering/blobtools", mode: params.publish_mode

    input:
    path btk
    path busco_table   // single BUSCO full_table.tsv, anthozoa_odb12

    output:
    path "btk", emit: blobdir
    path "versions.yml", emit: versions

    script:
    """
    blobtools add --busco ${busco_table} btk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "\$(blobtools --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    echo "BLOBTOOLS_ADD_BUSCO stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BLOBTOOLS_APPLY_FILTER
// Filter the assembly for one category (host/symbiont) and write the
// category FASTA plus a filtered BlobDir, summary and scaffold table.
//
// Selection is EITHER:
//   - a manual BlobtoolKit Viewer selection, exported from
//     `blobtools view --remote btk` as a JSON list file and supplied via
//     params.filter_json (optional; applied to the host assembly), or
//   - the automatic best-hit filter (params.auto_filter, default true):
//     keep the target taxon (+ no-hit for the host), drop all other hits.
// With neither, the subworkflow skips this process and only the BlobDir is
// produced.
// The BlobDir published by BLOBTOOLS_ADD_BUSCO is what the manual review
// is done on.
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
    // Manual JSON selection (kept as-is, no --invert) when supplied. Otherwise
    // the automatic filter, driven by the BEST HIT only: keep the target taxon,
    // and for the host also keep "no-hit" scaffolds (no evidence they are
    // contamination); everything with a hit to anything else is dropped.
    // (--param Keys=A,B alone EXCLUDES A,B; --invert flips it to "keep".)
    def selection = filter_json ? "--json ${filter_json}" : "--param ${taxon_rank}--Keys=\$KEYS --invert"
    def keys_cmd = filter_json ? "" : """
    # Keys to keep. "no-hit" is added for the host only, and only if it exists
    # in this BlobDir: blobtools errors on a key absent from the field.
    KEYS="${taxon_name}"
    ${category == 'host' ? """blobtools filter --table ${category}_ranks.tsv --table-fields ${taxon_rank} ${blobdir}
    # (blobtools writes CRLF line endings, so strip the \\r before comparing)
    if awk -F'\\t' '{sub(/\\r\$/, "")} NR>1 && \$3=="no-hit" {f=1} END {exit !f}' ${category}_ranks.tsv; then KEYS="\$KEYS,no-hit"; fi""" : ""}
    echo "${category}: keeping ${taxon_rank} = \$KEYS"
    """
    """
    ${keys_cmd}
    # blobtools writes <assembly stem>.<suffix>.fasta next to the input FASTA
    # (the task dir), NOT into --output; --output builds the filtered BlobDir.
    blobtools filter \\
      ${selection} \\
      --fasta ${assembly} \\
      --output ${category}_filtered \\
      --suffix ${category} \\
      --summary ${category}_summary.json \\
      --table ${category}_scaffolds.tsv \\
      --table-fields length \\
      ${blobdir}

    # Fail loudly if the selection matched nothing (e.g. a misnamed taxon or
    # empty JSON) instead of silently producing an empty downstream input.
    # (Counted from the table: the summary JSON always has a >0 "n" field.)
    N=\$(( \$(wc -l < ${category}_scaffolds.tsv) - 1 ))
    if [ "\$N" -le 0 ]; then
      echo "ERROR: selection for category '${category}' matched 0 scaffolds (${filter_json ? 'json ' + filter_json : taxon_rank + '=' + taxon_name})" >&2
      exit 1
    fi
    echo "${category}: \$N scaffolds selected"

    # <stem>.<category>.<ext>; the category suffix keeps host/symbiont apart.
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
// Render the static PNGs (blob, snail, cumulative) of the stage's FINAL output
// BlobDir, once, at the end of the blobtools stage: the filtered host BlobDir
// when a filter ran, else the unfiltered final BlobDir. Interactive review is done
// separately with `blobtools view --remote btk`.
//
// `--plot` makes blobtools render with blobtk, so no browser is needed
// (the 4.4.4 container ships no chromium; the --driver option is only used
// by the non-`--plot` viewer path and is ignored here).
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
    # blobtools view mkdirs HOME, which is read-only / absent inside the
    # container under Nextflow -- give it a writable one.
    export HOME=\$PWD

    for VIEW in blob snail cumulative; do
      blobtools view \\
        --plot \\
        --view \$VIEW \\
        --format png \\
        --out . \\
        ${blobdir}
    done

    # Require one PNG per view instead of publishing a silent partial set.
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
// Pass assembly through unchanged (when blobtools is skipped)
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
