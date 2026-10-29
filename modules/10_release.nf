// ============================================================
// 10_release.nf
// Release curation: the NCBI-ready release, built from the pipeline's own outputs:
//   product curation -> funannotate NCBI pass -> .tbl fixes -> table2asn (+ mitogenome)
//   repeat files, high-confidence subset, symbiont BUSCO, summary.tsv, release QC checks,
//   release folder <outdir>/release_curation/<release_name>/.
// Scripts: bin/curate_products.py, make_fix_files.py, fix_tbl.py, mito_tbl.py, drop_loci.py,
// sync_release_files.py, release_sets.py, genome_summary.py, annotation_qc_gate.py.
// The released sequence is the final masked assembly as it is: nothing edits it after annotation.
// ============================================================

nextflow.enable.dsl = 2

// "<flag> <file>", or nothing for a NO_* placeholder
def opt_arg(flag, f) {
    return f.baseName.startsWith('NO_') ? '' : "${flag} ${f}"
}

// NCBI linkage-evidence terms for scaffold gaps (table2asn/tbl2asn -l)
def ncbi_linkage_terms() {
    return ['paired-ends', 'align-genus', 'align-xgenus', 'align-trnscpt', 'within-clone', 'clone-contig',
            'map', 'strobe', 'unspecified', 'pcr', 'proximity-ligation']
}

// Linkage evidence of this run's scaffold gaps: params.ncbi_linkage_evidence ('_' read as '-'), otherwise
// align-genus when RagTag scaffolded against a reference (scaffold_reference) and unspecified when not.
def linkage_evidence() {
    if (params.ncbi_linkage_evidence) {
        return params.ncbi_linkage_evidence.toString().replace('_', '-')
    }
    return params.scaffold_reference ? 'align-genus' : 'unspecified'
}

// Product names that break NCBI rules (GO phrases, truncations, 'hypothetical protein' on a gene
// with a symbol, funannotate's need-curating list) are looked up in UniProt by gene symbol:
// the recommended name (action 'fix'), or 'hypothetical protein' without a symbol ('remove').
// Lookups are cached in uniprot_cache.json: a symbol in the cache is never queried again, so a rerun
// with the cache of an earlier run gives the same names (params.product_cache; by default the cache an
// earlier run published in outdir). params.uniprot_offline: cache only, no network.
// Run twice: 'pass1' before funannotate's NCBI pass (funannotate's own *.fixedproducts), 'final' with
// the names table2asn rewrote in the preflight run (FIX_TBL applies them to the .tbl).
process CURATE_PRODUCTS {
    label 'res_light'
    tag   "curate_products_${pass}"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir path: { "${params.outdir}/10_release_curation${params.release_folder_suffix}/product_curation" + (pass == 'final' ? '' : "/${pass}") }, mode: params.publish_mode,
               pattern: '{product_curation.tsv,uniprot_cache.json}'

    input:
    val  pass     // 'pass1' or 'final'
    val  lookup_all   // '--lookup-all' (params.product_lookup_all) or '': an input, so changing it reruns this step
    path anno     // funannotate annotate output (pass 1)
    path fixed, stageAs: 'fixed/fixed*.txt'   // *.fixedproducts (names tbl2asn/table2asn rewrote), or a NO_ placeholder
    path draft    // optional symbol -> product overrides, or NO_ placeholder
    path cache, stageAs: 'cache_in/*'   // UniProt cache from an earlier run or pass, or NO_ placeholder
    path funannotate_db                 // go.obo and interpro.tsv, staged so containers see them

    output:
    path "product_curation.tsv", emit: curation
    path "uniprot_cache.json",   emit: cache
    path "versions.yml",         emit: versions

    script:
    def offline = params.uniprot_offline.toString().toBoolean() ? '--offline' : ''
    """
    if [ -s ${cache} ]; then cp ${cache} uniprot_cache.json; fi
    ann=\$(ls ${anno}/annotate_results/*.annotations.txt | head -1)
    need=\$(ls ${anno}/annotate_results/Gene2Products.need-curating.txt 2>/dev/null || echo NO_NEED_CURATING)
    cat fixed/* > fixedproducts.txt
    # GO term and InterPro domain names (funannotate database): products equal to one are curated too
    curate_products.py ${offline} ${lookup_all} \$ann \$need ${draft} uniprot_cache.json product_curation.tsv fixedproducts.txt \
      ${funannotate_db}/go.obo ${funannotate_db}/interpro.tsv
    test -s product_curation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //')"
        uniprot: "\$(python3 -c "import json; print(', '.join(json.load(open('uniprot_cache.json')).get('#uniprot_release', [])) or 'cache only')")"
    END_VERSIONS
    """

    stub:
    """
    printf 'gene\\tsymbol\\told_product\\tnew_product\\taction\\tsource\\n' > product_curation.tsv
    echo '{}' > uniprot_cache.json
    printf '"${task.process}":\\n    uniprot: stub\\n' > versions.yml
    """
}

// --fix / --remove files keyed by the locus tags funannotate's NCBI pass will assign.
process MAKE_FIX_FILES {
    label 'res_light'
    tag   "make_fix_files"

    conda 'bioconda::funannotate=1.8.17'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/10_release_curation${params.release_folder_suffix}/product_curation", mode: params.publish_mode

    input:
    path refined_gff
    path fasta
    path curation

    output:
    path "funannotate_fix.tsv",    emit: fix
    path "funannotate_remove.txt", emit: remove
    path "locus_tag_map.tsv"
    path "make_fix_files.log"

    script:
    """
    make_fix_files.py ${refined_gff} ${fasta} ${params.locus_tag} ${curation} \\
      funannotate_fix.tsv funannotate_remove.txt locus_tag_map.tsv > make_fix_files.log 2>&1
    """

    stub:
    """
    touch funannotate_fix.tsv funannotate_remove.txt locus_tag_map.tsv make_fix_files.log
    """
}

// funannotate annotate with the NCBI settings: organism name, locus-tag prefix (--rename),
// submission template, isolate/strain, curated products (--fix/--remove).
process FUNANNOTATE_NCBI {
    tag "funannotate_ncbi"

    conda 'bioconda::funannotate=1.8.17'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/10_release_curation${params.release_folder_suffix}/funannotate_ncbi", mode: params.publish_mode

    input:
    path gff
    path assembly
    path eggnog
    path iprscan_xml
    path phobius_results
    path funannotate_db
    path fix_tsv
    path remove_txt
    path sbt

    output:
    path "anno", emit: dir
    path "versions.yml", emit: versions

    script:
    def args = []
    if (!iprscan_xml.baseName.startsWith('NO_'))     args << "--iprscan ${iprscan_xml}"
    if (!phobius_results.baseName.startsWith('NO_')) args << "--phobius ${phobius_results}"
    args << "--rename ${params.locus_tag}"
    args << "--sbt ${sbt}"
    if (params.ncbi_isolate)                     args << "--isolate '${params.ncbi_isolate}'"
    if (params.ncbi_strain)                      args << "--strain '${params.ncbi_strain}'"
    args << "--tbl2asn '-l ${linkage_evidence()}'"
    """
    export FUNANNOTATE_DB=\$(readlink -f ${funannotate_db})
    cp -r "\${AUGUSTUS_CONFIG_PATH:-/usr/local/config}" augustus_config
    export AUGUSTUS_CONFIG_PATH="\$PWD/augustus_config"
    fix=""; [ -s ${fix_tsv} ] && fix="--fix ${fix_tsv}"
    rem=""; [ -s ${remove_txt} ] && rem="--remove ${remove_txt}"
    funannotate annotate --gff ${gff} --fasta ${assembly} -s '${params.ncbi_species}' \\
      --busco_db ${params.funannotate_busco_db} --eggnog ${eggnog} \\
      ${args.join(' ')} \$fix \$rem --cpus ${task.cpus} -o anno

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "\$(funannotate --version 2>&1 | grep -m1 '^funannotate v' | sed 's/funannotate v//' || echo unknown)"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p anno/annotate_results
    for x in gff3 proteins.fa stats.json tbl scaffolds.fa annotations.txt cds-transcripts.fa mrna-transcripts.fa gbk; do
        touch anno/annotate_results/stub.\$x
    done
    printf '"${task.process}":\\n    funannotate: stub\\n' > versions.yml
    """
}

// Mechanical NCBI fixes to funannotate's .tbl (ncRNA_class, EC numbers, product punctuation,
// rRNA names and partial ends, gene cross-references). Run twice: 'preflight' without feedback, and
// 'final' with what the preflight table2asn run changed (.ecn, .fixedproducts) and the final curation.
process FIX_TBL {
    label 'res_light'
    tag   "fix_tbl_${pass}"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir path: { "${params.outdir}/10_release_curation${params.release_folder_suffix}/table2asn" + (pass == 'final' ? '' : "/${pass}") }, mode: params.publish_mode, pattern: '*.tsv'

    input:
    val  pass                                      // 'preflight' or 'final'
    val  strip                                     // '--strip-gene-suffix' (params.strip_gene_suffix) or '': an input, so changing it reruns this step
    path anno
    path refined_gff
    path ecn,           stageAs: 'feedback/ecn/*'   // preflight table2asn .ecn, or NO_ placeholder
    path curation,      stageAs: 'feedback/curation/*'   // final product_curation.tsv, or NO_ placeholder
    path fixedproducts, stageAs: 'feedback/fixed/*' // preflight table2asn .fixedproducts, or NO_ placeholder

    output:
    path "anno_fixed", emit: dir
    path "fix_tbl_report.tsv"
    path "rrna_decisions.tsv", emit: rrna     // rRNA fragments: partial ends or removed (RELEASE_SETS removes the same genes)

    script:
    """
    mkdir -p anno_fixed/annotate_results
    for f in ${anno}/annotate_results/*; do ln -s \$(readlink -f \$f) anno_fixed/annotate_results/; done
    base=\$(basename ${anno}/annotate_results/*.stats.json .stats.json)
    rm anno_fixed/annotate_results/\$base.tbl
    cat feedback/ecn/* > preflight.ecn
    cat feedback/fixed/* > preflight.fixedproducts
    fix_tbl.py ${strip} ${anno}/annotate_results/\$base.tbl ${refined_gff} anno_fixed/annotate_results/\$base.tbl fix_tbl_report.tsv \\
        ${anno}/annotate_results/\$base.scaffolds.fa rrna_decisions.tsv preflight.ecn ${curation} preflight.fixedproducts
    """

    stub:
    """
    mkdir -p anno_fixed/annotate_results
    cp -r ${anno}/annotate_results/. anno_fixed/annotate_results/
    touch fix_tbl_report.tsv rrna_decisions.tsv
    """
}

// Mitogenome FASTA (NCBI source modifiers) and feature table from the Liftoff annotation,
// for the same table2asn submission as the nuclear scaffolds.
process MITO_TBL {
    label 'res_light'
    tag   "mito_tbl"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/10_release_curation${params.release_folder_suffix}/table2asn/mitogenome", mode: params.publish_mode

    input:
    path mito_fasta
    path mito_gff
    path mito_check     // mito_gene_check.tsv: [completeness=complete] only when every protein is full length

    output:
    path "mt.fsa", emit: fsa
    path "mt.tbl", emit: tbl

    script:
    """
    mito_tbl.py ${mito_fasta} ${mito_gff} ${mito_check} ${params.locus_tag ?: 'MITO'} mt.fsa mt.tbl
    """

    stub:
    """
    touch mt.fsa mt.tbl
    """
}

// table2asn: .sqn, .val, .dr from funannotate's .fsa/.tbl (after FIX_TBL), plus the mitogenome.
// table2asn is NCBI's static Linux binary (params.table2asn: the file, or its .gz download).
// 'preflight': what table2asn changes (.ecn, .fixedproducts) goes back into FIX_TBL; 'final': the submission.
process TABLE2ASN {
    label 'res_mid_4h'
    tag   "table2asn_${pass}"

    publishDir path: { "${params.outdir}/10_release_curation${params.release_folder_suffix}/table2asn" + (pass == 'final' ? '' : "/${pass}") }, mode: params.publish_mode

    input:
    val  pass                             // 'preflight' or 'final'
    path anno
    path table2asn
    path sbt
    path cmt                              // structured comment for -w
    path mito_fsa, stageAs: 'mito/mt.fsa' // empty placeholder when there is no mitogenome
    path mito_tbl, stageAs: 'mito/mt.tbl'

    output:
    path "ncbi/*", emit: files
    path "versions.yml", emit: versions

    script:
    def lin = params.ncbi_lineage ? " [lineage=${params.ncbi_lineage}]" : ''
    def iso = params.ncbi_isolate ? " [isolate=${params.ncbi_isolate}]" : (params.ncbi_strain ? " [strain=${params.ncbi_strain}]" : '')
    def bp  = params.bioproject ? " [bioproject=${params.bioproject}]" : ''
    def bs  = params.biosample ? " [biosample=${params.biosample}]" : ''
    def t_arg = "-t ${sbt}"
    def w_arg = "-w ${cmt}"
    """
    if grep -q '<FILL' ${cmt}; then
        echo "ERROR: structured comment ${cmt} still has <FILL> placeholders" >&2; exit 1
    fi
    case ${table2asn} in
        *.gz) gunzip -c ${table2asn} > t2a ;;
        *)    cp -L ${table2asn} t2a ;;
    esac
    chmod +x t2a
    mkdir -p in ncbi
    base=\$(basename ${anno}/annotate_results/*.scaffolds.fa .scaffolds.fa)
    cp -L ${anno}/annotate_results/\$base.scaffolds.fa in/\$base.fsa
    cp -L ${anno}/annotate_results/\$base.tbl in/\$base.tbl
    # NCBI: organelle sequences belong in the genome's submission file
    if [ -s mito/mt.fsa ]; then
        cat mito/mt.fsa >> in/\$base.fsa
        cat mito/mt.tbl >> in/\$base.tbl
    fi
    rc=0
    ./t2a -indir in -outdir ncbi ${t_arg} -M n -J -euk -gaps-min 10 -gaps-unknown 100 -gap-type scaffold \\
      -l ${linkage_evidence()} \\
      -j "[organism=${params.ncbi_species}] [tech=wgs]${iso}${lin}${bp}${bs}" ${w_arg} \\
      -T -Z -split-dr -V vb || rc=\$?
    # table2asn exits 2 when validation reports errors; the QC gate judges .val/.dr (G10)
    [ "\$rc" -eq 0 ] || [ "\$rc" -eq 2 ]
    ls ncbi/*.sqn > /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        table2asn: "\$(./t2a -version 2>&1 | head -1)"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ncbi
    touch ncbi/stub.sqn ncbi/stub.val ncbi/stub.dr ncbi/stub.gbf
    printf '"${task.process}":\\n    table2asn: stub\\n' > versions.yml
    """
}

// The release annotation files, made to agree with the submission: rRNA fragments FIX_TBL removed are
// removed, then GFF3 and annotations.txt take gene spans, symbols, products and EC numbers from the
// final .tbl (sync_release_files.py), and the .gbk is table2asn's GenBank flatfile of the .sqn.
process RELEASE_SETS {
    label 'res_light'
    tag   "release_sets"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/10_release_curation${params.release_folder_suffix}/sets", mode: params.publish_mode

    input:
    path anno
    path fixed, stageAs: 'fixed/*'   // FIX_TBL 'final' output: the .tbl table2asn submits
    path gbf                         // table2asn's .gbf of the final .sqn
    path refined_gff
    path tags
    path splits
    path rrna           // rrna_decisions.tsv (FIX_TBL): genes removed from the submission are removed here too

    output:
    path "sets/*", emit: files
    path "sync_release_files.tsv"

    script:
    """
    base=\$(basename ${anno}/annotate_results/*.gff3 .gff3)
    mkdir -p sets final dropped
    for x in gff3 annotations.txt; do
        drop_loci.py ${rrna} ${anno}/annotate_results/\$base.\$x dropped/\$base.\$x
    done
    sync_release_files.py ${fixed}/annotate_results/\$base.tbl dropped/\$base.gff3 dropped/\$base.annotations.txt \\
      final/\$base.gff3 final/\$base.annotations.txt sync_release_files.tsv
    cp final/\$base.annotations.txt sets/
    cp ${gbf} sets/\$base.gbk
    [ -s ${rrna} ] && cp ${rrna} sets/rrna_decisions.tsv
    release_sets.py --final-gff final/\$base.gff3 --final-proteins ${anno}/annotate_results/\$base.proteins.fa \\
      --annotations final/\$base.annotations.txt --refined-gff ${refined_gff} --tags ${tags} \\
      --splits ${splits} --prefix \$base --outdir sets
    """

    stub:
    """
    mkdir -p sets
    touch sets/stub.gff3 sets/stub.hc.gff3 sets/stub.annotations.txt sets/stub.gbk sets/id_map.tsv sets/release_sets_report.tsv sync_release_files.tsv
    """
}

// BUSCO on the symbiont bin with a lineage that fits it (e.g. alveolata_odb12); the host lineage
// gives a meaningless 0. params.symbiont_busco_lineage: a lineage directory or its .tar.gz.
process SYMBIONT_BUSCO {
    tag "symbiont_busco"

    conda 'bioconda::busco=6.1.0'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/10_release_curation${params.release_folder_suffix}/symbiont_busco", mode: params.publish_mode, pattern: 'short_summary*'

    input:
    path fasta
    path lineage

    output:
    path "short_summary*.txt", emit: summary
    path "short_summary*.json", optional: true
    path "versions.yml", emit: versions

    script:
    """
    if [ -d ${lineage} ]; then
        L=\$(readlink -f ${lineage})
    else
        mkdir -p lineages && tar -xzf ${lineage} -C lineages
        L=\$(ls -d \$PWD/lineages/*/ | head -1)
    fi
    if [ ! -s ${fasta} ]; then
        echo "symbiont bin is empty" > short_summary.empty.txt
    else
        busco -i ${fasta} -m genome -l \$L -o symbiont_busco --offline -c ${task.cpus}
        cp symbiont_busco/short_summary* .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: "\$(busco --version 2>&1 | sed 's/BUSCO //')"
    END_VERSIONS
    """

    stub:
    """
    touch short_summary.stub.txt
    printf '"${task.process}":\\n    busco: stub\\n' > versions.yml
    """
}

// Release QC checks (bin/annotation_qc_gate.py). Any FAIL stops the run.
process ANNOTATION_QC_GATE {
    label 'res_mid_4h'
    tag   "qc_gate"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/10_release_curation${params.release_folder_suffix}/qc", mode: params.publish_mode

    input:
    path release_gff,      stageAs: 'g/release.gff3'
    path refined_gff,      stageAs: 'g/refined.gff3'
    path refine_report,    stageAs: 'g/refine_report.tsv'
    path protein_identity, stageAs: 'g/protein_identity.tsv'
    path sets_report,      stageAs: 'g/release_sets_report.tsv'
    path final_proteins,   stageAs: 'g/final.proteins.fa'
    path refined_proteins, stageAs: 'g/refined.proteins.faa'
    path id_map,           stageAs: 'g/id_map.tsv'
    path gt_reports,       stageAs: 'g/gt/*'
    path stats_json,       stageAs: 'g/final.stats.json'
    path summary_tsv,      stageAs: 'g/summary.tsv'
    path name_map,         stageAs: 'g/map/*'      // scaffold_name_map.tsv with RagTag, else NO_ placeholder
    path agp,              stageAs: 'g/agp/*'      // RagTag AGP, else NO_ placeholder
    path fasta,            stageAs: 'g/scaffolds.fa'
    path contigs,          stageAs: 'g/contigs/*'  // the contigs that entered scaffolding (after the screens)
    path repeat_lib,       stageAs: 'g/repeats_classified.lib'
    path relabel_report,   stageAs: 'g/relabel_report.tsv'
    path screen_files,     stageAs: 'g/screen/*'   // screen_actions.tsv, screen_review.tsv (or NO_ placeholder)
    path ncbi_files,       stageAs: 'g/ncbi/*'
    path sbt,              stageAs: 'g/sbt/*'

    output:
    path "annotation_qc.tsv", emit: tsv

    script:
    def ncbi = (params.ncbi_isolate ? "--isolate '${params.ncbi_isolate}'" : '') +
               (params.ncbi_strain ? " --strain '${params.ncbi_strain}'" : '')
    def chk = [params.n_chromosome_scaffolds ? "--n-placed ${params.n_chromosome_scaffolds}" : '',
               params.qc_chrom_prefix ? "--chrom-prefix '${params.qc_chrom_prefix}'" : '',
               params.qc_expected_lineage ? "--expected-lineage '${params.qc_expected_lineage}'" : '',
               params.qc_expected_division ? "--expected-division '${params.qc_expected_division}'" : ''].findAll().join(' ')
    def has_map = !name_map.baseName.startsWith('NO_')
    def scaffold_args = has_map
        ? "--name-map ${name_map} --agp ${agp} --gap-evidence ${linkage_evidence().replace('-', '_')} " +
          "--gap-length ${params.ragtag_infer_gaps.toString().toBoolean() ? 0 : 100}"
        : ''
    """
    gts=\$(for f in g/gt/*; do printf -- '--gt-report %s ' \$f; done)
    annotation_qc_gate.py --release-gff g/release.gff3 --refined-gff g/refined.gff3 --refine-report g/refine_report.tsv \\
      --protein-identity g/protein_identity.tsv --release-sets-report g/release_sets_report.tsv \\
      --final-proteins g/final.proteins.fa --refined-proteins g/refined.proteins.faa --id-map g/id_map.tsv \\
      \$gts --stats-json g/final.stats.json --summary-tsv g/summary.tsv ${scaffold_args} \\
      --fasta g/scaffolds.fa --contigs ${contigs} --repeat-lib g/repeats_classified.lib --relabel-report g/relabel_report.tsv \\
      --screen-actions g/screen/screen_actions.tsv --screen-review g/screen/screen_review.tsv --val-dir g/ncbi \\
      --sbt ${sbt} --locus-tag '${params.locus_tag}' ${ncbi} ${chk} \\
      ${params.run_screens.toString().toBoolean() ? '' : '--screens-not-run'} --out annotation_qc.tsv
    """

    stub:
    """
    printf 'check\\tlevel\\tstatus\\tvalue\\tdetail\\n' > annotation_qc.tsv
    """
}

// The release folder: <outdir>/release_curation/<release_name>/ with genome/, annotation/, repeats/,
// ncbi/, mitogenome/, qc/ (summary.tsv among the QC files). PUBLISH_MULTIQC adds the MultiQC report.
process RELEASE_ASSEMBLE {
    label 'res_light'
    tag   "release_assemble"

    publishDir "${params.outdir}/release_curation${params.release_folder_suffix}", mode: 'copy'

    input:
    path qc_gate,       stageAs: 'in/gate/*'   // ANNOTATION_QC_GATE result: the folder is built only after the gate passed
    path fasta,         stageAs: 'in/genome/scaffolds.fa'
    path name_map,      stageAs: 'in/genome/map/*'   // with RagTag; else NO_ placeholder
    path anno,          stageAs: 'in/anno'
    path sets,          stageAs: 'in/sets/*'
    path repeat_files,  stageAs: 'in/repeats/*'
    path ncbi_files,    stageAs: 'in/ncbi/*'
    path mito_files,    stageAs: 'in/mito/*'
    path qc_files,      stageAs: 'in/qc/*'

    output:
    path "${params.release_name}", emit: dir

    script:
    def D = params.release_name
    """
    set -euo pipefail
    mkdir -p ${D}/{genome,annotation/extra,repeats,ncbi,mitogenome,qc}
    base=\$(basename in/anno/annotate_results/*.proteins.fa .proteins.fa)
    cp -L in/genome/scaffolds.fa ${D}/genome/\$base.scaffolds.fa
    for f in in/genome/map/*; do case \$(basename \$f) in NO_*) ;; *) cp -L \$f ${D}/genome/ ;; esac; done
    cp -L in/sets/\$base.gff3 ${D}/annotation/
    for x in proteins.fa cds-transcripts.fa mrna-transcripts.fa; do
        cp -L in/anno/annotate_results/\$base.\$x ${D}/annotation/
    done
    cp -L in/sets/\$base.annotations.txt in/sets/\$base.gbk ${D}/annotation/   # rRNA fragments removed (RELEASE_SETS)
    for f in in/sets/*; do
        case \$(basename \$f) in \$base.gff3|\$base.annotations.txt|\$base.gbk|release_sets_report.tsv) ;; *) cp -L \$f ${D}/annotation/extra/ ;; esac
    done
    cp -L in/repeats/* ${D}/repeats/
    cp -L in/ncbi/* ${D}/ncbi/
    for f in in/mito/* in/qc/*; do
        [ -e "\$f" ] || continue                     # an empty input leaves the pattern unexpanded
        case \$(basename \$f) in NO_*) continue ;; esac
        case \$f in in/mito/*) cp -L \$f ${D}/mitogenome/ ;; *) cp -L \$f ${D}/qc/ ;; esac
    done
    """

    // no stub: the script only copies files, so -stub runs check the folder layout
}
