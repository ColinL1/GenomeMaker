// ============================================================
// 08_functional_annotation.nf
// Functional annotation pipeline:
// Phobius → InterProScan (funannotate iprscan) → eggNOG-mapper → funannotate annotate
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PHOBIUS_RUN
// Signal peptide and transmembrane topology prediction.
//
// PORTABILITY: Phobius is licence-restricted, so no public image may legally
// ship it -- but that is not a reason to abandon the container. phobius.pl is
// plain Perl plus a small bundled binary, so it runs in a stock perl image with
// the LICENSED DIRECTORY STAGED AS A PATH INPUT. Anyone with their own licence
// drops it in and the container is unchanged. This is the same pattern the
// pipeline already uses for blast_db, diamond_db and the BUSCO offline dir.
//
// phobius.pl resolves its model and decodeanhmm relative to its own location
// (FindBin), so the whole directory must be staged, not just the script.
// ============================================================
process PHOBIUS_RUN {
    label 'functional_env'
    tag   "phobius"

    conda 'conda-forge::perl'
    container 'quay.io/biocontainers/perl:5.32'
    publishDir "${params.outdir}/08_functional_annotation/phobius", mode: params.publish_mode

    input:
    path proteins
    path phobius_dir

    output:
    path "phobius_results.txt", emit: phobius_results
    path "versions.yml",        emit: versions

    script:
    """
    # Phobius rejects the '*' stop characters gffread emits for partial models.
    sed '/^>/!s/\\*//g' ${proteins} > proteins_clean.faa

    # phobius.pl calls ./decodeanhmm next to itself, which is the 32-bit binary and
    # cannot run in a 64-bit container (Phobius then dies with the misleading
    # "Could not read provided fasta sequence"). Run from a local copy of the
    # install with the shipped decodeanhmm.64bit in its place; the licensed
    # install itself is not modified.
    cp -rL ${phobius_dir} phobius_local
    if [ -x phobius_local/decodeanhmm.64bit ]; then
        cp phobius_local/decodeanhmm.64bit phobius_local/decodeanhmm
    fi

    # -short: funannotate annotate --phobius parses the one-line-per-protein short
    # format (ID TM SP PREDICTION); the default long format crashes its parser.
    perl phobius_local/phobius.pl -short proteins_clean.faa > phobius_results.txt

    if [ ! -s phobius_results.txt ]; then
        echo "ERROR: Phobius produced no output. Check the licensed install in" >&2
        echo "       ${phobius_dir} -- it must contain phobius.pl, phobius.model," >&2
        echo "       phobius.options and decodeanhmm." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phobius: "1.01"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch phobius_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phobius: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: INTERPROSCAN
// InterProScan run directly, not through funannotate's `-m docker` wrapper --
// that made funannotate shell out to `docker run` from inside a container, which
// cannot work under Apptainer on a compute node.
//
// PORTABILITY. InterProScan needs ~45 GB of analysis data that no image carries,
// so the DATA always comes from the host. The SOFTWARE still comes from a
// container by default. Two supported modes:
//
//   container mode (default, portable): params.interproscan_container supplies
//     the software and params.interproscan_data the analysis data. InterProScan
//     validates that the two are the SAME RELEASE, so they must match.
//
//   local mode (params.interproscan_dir set): use an existing site install and
//     bypass the container. Chosen on this system because a complete 5.63-95.0
//     install already exists and no public image matches that release
//     (biocontainers stops at 5.59-91.0; interpro/interproscan starts at
//     5.73-104.0), so container mode would mean re-downloading ~50 GB of
//     version-matched data for no scientific gain.
// ============================================================
process INTERPROSCAN {
    label 'functional_env'
    tag   "interproscan"

    conda 'bioconda::interproscan'
    container "${params.interproscan_dir ? '' : params.interproscan_container}"
    // Container mode: the official image expects its analysis data at
    // /opt/interproscan/data (docs: -v <data>:/opt/interproscan/data).
    containerOptions "${params.interproscan_dir ? '' : "-B ${params.interproscan_data}:/opt/interproscan/data"}"
    publishDir "${params.outdir}/08_functional_annotation/interproscan", mode: params.publish_mode

    input:
    path proteins
    path ipr_dir     // local install dir, or the staged data dir in container mode

    output:
    path "iprscan.xml", emit: iprscan_xml
    path "iprscan.tsv", emit: iprscan_tsv, optional: true
    path "versions.yml", emit: versions

    script:
    def ipr_exe = params.interproscan_dir ? "${params.interproscan_dir}/interproscan.sh" : "/opt/interproscan/interproscan.sh"
    """
    # InterProScan rejects '*' stop characters, which gffread emits for partial
    # models -- strip them rather than letting a long run die late.
    sed '/^>/!s/\\*//g' ${proteins} > proteins_clean.faa

    ${ipr_exe} \\
      --input proteins_clean.faa \\
      --seqtype p \\
      --formats XML,TSV \\
      --output-file-base iprscan \\
      --cpu ${task.cpus} \\
      --goterms --pathways \\
      --disable-precalc \\
      --tempdir \$PWD/ipr_temp
    # --tempdir: the default is temp/ inside the InterProScan install, which is
    # read-only in the container (and would write into a shared install locally).

    if [ ! -s iprscan.xml ]; then
        echo "ERROR: InterProScan produced no XML." >&2
        echo "       In container mode, params.interproscan_data must be the SAME" >&2
        echo "       RELEASE as params.interproscan_container -- InterProScan refuses" >&2
        echo "       to run on mismatched software/data versions." >&2
        echo "       In local mode, ${ipr_dir} must hold interproscan.sh and a" >&2
        echo "       populated data/ directory (~45 GB)." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: "\$(${ipr_exe} --version 2>&1 | grep -oim1 'version [0-9.-]*' | sed 's/version //' || echo unknown)"
        container: "${task.container ?: params.interproscan_dir ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf '<?xml version="1.0"?>\\n<protein-matches/>\\n' > iprscan.xml
    touch iprscan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: "stub"
    END_VERSIONS
    """
}

// ============================================================
// NOTE: IPRSCAN_RUN (funannotate iprscan -m docker) was REMOVED. It launched
// the InterProScan *Docker* image from inside a container, which cannot work
// under Apptainer/SLURM. Replaced by the INTERPROSCAN process below.
// ============================================================

// ============================================================
// Process: EGGNOG_MAPPER
// Run eggNOG-mapper for orthology and GO annotations
// ============================================================
process EGGNOG_MAPPER {
    label 'functional_env'
    tag   "eggnog"

    conda 'bioconda::eggnog-mapper'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0'
    publishDir "${params.outdir}/08_functional_annotation/eggnog", mode: params.publish_mode

    input:
    path proteins
    val  eggnog_db

    output:
    path "eggnog.emapper.annotations", emit: eggnog_anno
    path "eggnog.emapper.seed_orthologs", emit: eggnog_ortho
    path "versions.yml", emit: versions

    script:
    // --data_dir must point at the eggNOG database, not the funannotate one --
    // it was previously handed params.funannotate_db, a different database
    // entirely, so emapper either failed to find its data or silently searched
    // the wrong thing.
    //
    // `--no_egeny` was not a valid emapper flag (no such option in 2.1.x) and has
    // been dropped. Annotation is the point of running this, so there is nothing
    // here that should be suppressed.
    """
    emapper.py \\
      --cpu ${task.cpus} \\
      -m mmseqs \\
      --data_dir ${eggnog_db} \\
      -i ${proteins} \\
      -o eggnog \\
      --override
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: "\$(emapper.py --version 2>&1 | grep -oim1 'emapper-[0-9.]*' | sed 's/emapper-//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch eggnog.emapper.annotations eggnog.emapper.seed_orthologs
    echo "EGGNOG_MAPPER stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: FUNANNOTATE_ANNOTATE
// Final functional annotation synthesis
// ============================================================
process FUNANNOTATE_ANNOTATE {
    label 'functional_env'
    tag   "funannotate_annotate"

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/08_functional_annotation/funannotate", mode: params.publish_mode

    input:
    path gff3
    path assembly
    path eggnog_anno
    path iprscan_xml
    path phobius_results
    path funannotate_db   // configured funannotate database (funannotate setup)

    output:
    path "final_anno", emit: anno_dir
    // product names tbl2asn had to rewrite; written to the task dir, read by the release's product curation
    path "*.fixedproducts", emit: fixedproducts, optional: true
    path "versions.yml", emit: versions

    script:
    // InterProScan and Phobius are OPTIONAL evidence: InterProScan needs a ~50 GB
    // data directory that may not be installed, and Phobius is licence-restricted
    // so its image is not always pullable. Both arrive as a NO_FILE placeholder
    // when skipped -- the same pattern MULTIQC_FINAL uses for its optional config.
    def ipr_arg = (iprscan_xml     && !iprscan_xml.name.startsWith('NO_'))     ? "--iprscan ${iprscan_xml}"     : ""
    def pho_arg = (phobius_results && !phobius_results.name.startsWith('NO_')) ? "--phobius ${phobius_results}" : ""
    //
    // --busco_db must be a lineage inside FUNANNOTATE_DB. funannotate 1.x only
    // supports its own odb9-era sets (metazoa, dikarya, ...), not odb12, so it
    // has its own param; the pipeline's BUSCO/compleasm QC keeps busco_lineage.
    """
    # funannotate finds its database only through FUNANNOTATE_DB; without it,
    # annotate aborts with "Funannotate database not properly configured".
    export FUNANNOTATE_DB=\$(readlink -f ${funannotate_db})

    # funannotate's BUSCO step runs AUGUSTUS and needs AUGUSTUS_CONFIG_PATH (unset
    # in the image), and it creates a species there with new_species.pl -- the
    # image's /usr/local/config is read-only, so use a writable copy.
    cp -r "\${AUGUSTUS_CONFIG_PATH:-/usr/local/config}" augustus_config
    export AUGUSTUS_CONFIG_PATH="\$PWD/augustus_config"

    funannotate annotate \\
      --gff ${gff3} \\
      --fasta ${assembly} \\
      -s ${params.species_name} \\
      --busco_db ${params.funannotate_busco_db} \\
      --eggnog ${eggnog_anno} \\
      ${ipr_arg} \\
      ${pho_arg} \\
      --cpus ${task.cpus} \\
      -o final_anno
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "\$(funannotate --version 2>&1 | sed 's/funannotate v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p final_anno/annotate_results
    for x in gff3 proteins.fa stats.json tbl scaffolds.fa annotations.txt cds-transcripts.fa mrna-transcripts.fa gbk; do
        touch final_anno/annotate_results/stub.\$x
    done
    touch final_anno/annotate_results/Gene2Products.need-curating.txt
    echo "FUNANNOTATE_ANNOTATE stub completed" > final_anno/annotate.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "stub"
    END_VERSIONS
    """
}

// ============================================================
// NOTE: a FINAL_BUSCO process lived here. It was dead code -- never
// imported or called by subworkflows/08_functional_annotation.nf -- and
// carried three bugs: BUSCO output paths that BUSCO never creates, a GFF3
// passed to -m transcriptome (which needs a transcript FASTA), and
// --offline given an argument when it takes none. Annotation-level
// completeness is covered by BUSCO_FINAL (stage 04) plus the compleasm and
// OMArk QC added in stage 07.
// ============================================================

