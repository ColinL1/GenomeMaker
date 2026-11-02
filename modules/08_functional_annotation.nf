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
// Licence-restricted: the licensed directory is staged as a path input into a
// stock perl image. phobius.pl locates its files via FindBin, so stage the whole directory.
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

    # ./decodeanhmm is 32-bit; run a local copy with decodeanhmm.64bit in its place
    cp -rL ${phobius_dir} phobius_local
    if [ -x phobius_local/decodeanhmm.64bit ]; then
        cp phobius_local/decodeanhmm.64bit phobius_local/decodeanhmm
    fi

    # -short: funannotate's --phobius parser needs the short format
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
// InterProScan run directly (funannotate's `-m docker` cannot work under Apptainer).
// The analysis data always comes from the host. Container mode (default): the
// data must match params.interproscan_container's release. Local mode:
// params.interproscan_dir points at a site install.
// ============================================================
process INTERPROSCAN {
    label 'functional_env'
    tag   "interproscan"

    conda 'bioconda::interproscan'
    container "${params.interproscan_dir ? '' : params.interproscan_container}"
    containerOptions "${params.interproscan_dir ? '' : "-B ${params.interproscan_data}:/opt/interproscan/data"}"
    publishDir "${params.outdir}/08_functional_annotation/interproscan", mode: params.publish_mode

    input:
    path proteins
    path ipr_dir     // local install, or the data dir in container mode

    output:
    path "iprscan.xml", emit: iprscan_xml
    path "iprscan.tsv", emit: iprscan_tsv, optional: true
    path "versions.yml", emit: versions

    script:
    def ipr_exe = params.interproscan_dir ? "${params.interproscan_dir}/interproscan.sh" : "/opt/interproscan/interproscan.sh"
    """
    # InterProScan rejects '*' stop characters (gffread emits them for partial models)
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
    # --tempdir: the default is inside the (read-only) install

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
    path funannotate_db

    output:
    path "final_anno", emit: anno_dir
    // product names tbl2asn rewrote, read by the release's product curation
    path "*.fixedproducts", emit: fixedproducts, optional: true
    path "versions.yml", emit: versions

    script:
    // InterProScan and Phobius are optional (NO_FILE placeholder when skipped)
    def ipr_arg = (iprscan_xml     && !iprscan_xml.baseName.startsWith('NO_'))     ? "--iprscan ${iprscan_xml}"     : ""
    def pho_arg = (phobius_results && !phobius_results.baseName.startsWith('NO_')) ? "--phobius ${phobius_results}" : ""
    """
    export FUNANNOTATE_DB=\$(readlink -f ${funannotate_db})

    # funannotate's BUSCO step writes into AUGUSTUS_CONFIG_PATH, read-only in the image
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
        funannotate: "\$(funannotate --version 2>&1 | grep -m1 '^funannotate v' | sed 's/funannotate v//' || echo unknown)"
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
