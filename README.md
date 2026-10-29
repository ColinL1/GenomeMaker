# GenomeMaker - Coral Genome Assembly & Annotation Pipeline 
**Not final**

[![Alpha Release](https://img.shields.io/badge/release-alpha-orange)](https://github.com/ColinL1/GenomeMaker/releases)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.04.0-%233099C1)](https://www.nextflow.io/)
[![DSL2](https://img.shields.io/badge/DSL2-Compatible-brightgreen)](https://www.nextflow.io/docs/latest/dsl2.html)

[Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) pipeline for *de novo* assembly and annotation of coral (Scleractinia) genomes from Oxford Nanopore Technologies (ONT) long reads, from raw reads to an NCBI-ready submission.

## Table of Contents

- [GenomeMaker - Coral Genome Assembly \& Annotation Pipeline](#genomemaker---coral-genome-assembly--annotation-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Pipeline Overview](#pipeline-overview)
  - [Installation Requirements](#installation-requirements)
    - [Tools set up outside the containers](#tools-set-up-outside-the-containers)
  - [Reference Data and Databases](#reference-data-and-databases)
  - [Quick Start](#quick-start)
    - [1. Clone the repository](#1-clone-the-repository)
    - [2. Stub test (no computation)](#2-stub-test-no-computation)
    - [3. De novo run](#3-de-novo-run)
    - [4. Starting from an existing assembly](#4-starting-from-an-existing-assembly)
    - [Profiles](#profiles)
  - [Parameter Reference](#parameter-reference)
    - [Core inputs](#core-inputs)
    - [Assembly and filtering](#assembly-and-filtering)
    - [Annotation](#annotation)
    - [Executor](#executor)
  - [Multi-Assembler Comparison](#multi-assembler-comparison)
  - [Stage Toggles and Dependencies](#stage-toggles-and-dependencies)
  - [Contamination Screens](#contamination-screens)
  - [NCBI Release](#ncbi-release)
  - [Output Structure](#output-structure)
    - [Release folder](#release-folder)
  - [License](#license)

---

## Pipeline Overview

| Stage | Module | Description |
|-------|--------|-------------|
| **00** | `00_preprocess_reads` | Porechop adapter removal, FastQC/NanoPlot/NanoStat QC, chopper read sets (assembly, polishing, high-quality k-mer reads), optional symbiont read removal (chopper `--contam`) |
| **01** | `01_kmer_genome_size` | Meryl k-mer database and histogram, GenomeScope 2 genome size estimate |
| **02** | `02_multi_assembler` | Parallel assembly with NECAT, Flye, hifiasm and/or Raven; BUSCO, gfastats and QUAST per assembly |
| **02b** | `02b_select_best_assembly` | Ranks the assemblies (BUSCO completeness with a duplication penalty, N50 within a tie band, genome-size outliers ranked last) and picks one |
| **02c** | `02c_input_assembly_qc` | Replaces 01-02b when an existing assembly is given with `--input_assembly` |
| **03b** | `03b_purge_dups` | Haplotig removal with purge_dups (ONT coverage), followed by a BUSCO + gfastats checkpoint |
| **03** | `03_assembly_filtering` | BlobToolKit host filter (BLASTN vs nt, DIAMOND vs UniProt, ONT coverage, BUSCO); splits off a symbiont bin |
| **03c** | `03c_screens` | NCBI contamination screens: FCS-adaptor, VecScreen (UniVec), FCS-GX; unclear hits are reviewed before actions are applied |
| **03d** | `03d_scaffold` | Optional reference-guided scaffolding with RagTag |
| **04** | `04_assembly_polish` | Medaka polishing, funannotate clean/sort; Merqury QV, BUSCO and gfastats before and after; the unpolished assembly is kept if polishing makes it worse |
| **05** | `05_mitogenome` | Mitogenome with MitoHiFi (ONT) or GetOrganelle (Illumina DNA), Medaka polishing, rotation, Liftoff and MITOS2 annotation, adaptor/vector screens |
| **06** | `06_repeats` | RepeatModeler + EDTA library, RepeatClassifier for unknown families, RepeatMasker soft-masking, barrnap rRNA |
| **07** | `07_structural_annotation` | RNA-seq trimming and STAR mapping, BRAKER4 gene prediction, tRNAscan-SE, Infernal/Rfam ncRNAs, StringTie, AGAT merge and GFF3 validation, compleasm and OMArk |
| **07c** | `07c_annotation_refine` | Gene-model refinement: UTRs, gene membership, tRNA/ncRNA filters, haplotig and gap-fusion scans, isoform filter, protein identity check |
| **08** | `08_functional_annotation` | InterProScan, eggNOG-mapper, Phobius, funannotate annotate |
| **10** | `10_release` | Optional: product-name curation, table2asn NCBI submission, release QC checks, release folder with README |
| **09** | `09_summary` | Summary tables, software versions and an integrated MultiQC report (runs last, after the release) |

---

## Installation Requirements

| Requirement | Version | Notes |
|-------------|---------|-------|
| [Nextflow](https://www.nextflow.io/) | >= 24.04.0 | `resourceLimits` in the binac2 profile needs 24.04 |
| [Apptainer](https://apptainer.org/) / Singularity | >= 1.1.0 | Used by the `binac2` and `local` profiles; also needed on `PATH` for the FCS images and BRAKER4 |
| [Conda](https://docs.conda.io/) / Mamba | >= 22.0 | `conda` profile, and the Snakemake environment that runs BRAKER4 |
| [Git](https://git-scm.com/) | any | For cloning the repository |

```bash
# Install Nextflow
conda install -c bioconda nextflow
# or
curl -s https://get.nextflow.io | bash

nextflow -v
```

### Tools set up outside the containers

| Tool | Parameter | Notes |
|------|-----------|-------|
| BRAKER4 | `--braker4_dir`, `--braker4_snakemake_env` | A checkout of [Gaius-Augustus/BRAKER4](https://github.com/Gaius-Augustus/BRAKER4) and a conda env with Snakemake 8 and pandas |
| NECAT | `--necat_path_bin` | Only when `necat` is in `--assemblers` |
| Phobius | `--phobius_dir` | Licence-restricted; install from the licensed tarball, or `--run_phobius false` |
| InterProScan data | `--interproscan_data` | Must match the release of `--interproscan_container` (default `interpro/interproscan:5.73-104.0`) |
| FCS-adaptor / FCS-GX | `--fcs_adaptor_image`, `--fcs_gx_image` | NCBI `.sif` images, run with `apptainer exec` |
| table2asn | `--table2asn` | NCBI Linux binary, for the release stage only |

---

## Reference Data and Databases

Most defaults in `nextflow.config` point into a `config/` folder in the pipeline directory. That folder holds site-specific data and is **not part of the repository**: create it, or pass the paths on the command line or in your own config file.

| Parameter | Content |
|-----------|---------|
| `--blast_db` | BLAST nucleotide database (e.g. `core_nt`) with `new_taxdump/` next to it |
| `--diamond_db` | DIAMOND database of UniProt reference proteomes |
| `--busco_offline_dir` | BUSCO lineage datasets (null = download) |
| `--protein_refs` | Protein evidence for BRAKER4 (e.g. OrthoDB Metazoa) |
| `--rnaseq_dir` | Directory of paired-end RNA-seq FASTQ files |
| `--trimmomatic_adapters` | Adapter FASTA for RNA-seq trimming |
| `--funannotate_db`, `--eggnog_db` | funannotate and eggNOG-mapper databases |
| `--omamer_db` | OMAmer database for OMArk |
| `--rfam_cm`, `--rfam_clanin`, `--rfam_family` | Rfam covariance models (cmpress'd), clan file and family table |
| `--mito_ref` | Reference mitogenome for read recruitment |
| `--mito_annot_ref_fasta`, `--mito_annot_ref_gff` | Annotated related mitogenome for Liftoff (`false` = MITOS2 only) |
| `--bin_path` | Folder holding `rank.py` and `combine_comparison.py`; set it to `<pipeline>/bin/` on a fresh clone |

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/ColinL1/GenomeMaker.git
cd GenomeMaker
```

### 2. Stub test (no computation)

```bash
nextflow run . -stub -profile test
```

Runs every process in stub mode, including the optional stages and the release, with the dummy inputs in `assets/test/`. It checks the workflow graph and channel wiring without tools, data or containers.

### 3. De novo run

```bash
nextflow run main.nf -profile binac2 \
    --ont_raw_reads data/ONT_raw_reads.fastq.gz \
    --outdir results \
    --species_name Porites_sp \
    --run_screens false
```

The contamination screens are on by default and need their databases (see [Contamination Screens](#contamination-screens)); turn them off as above, or provide them.

### 4. Starting from an existing assembly

```bash
nextflow run main.nf -profile binac2 \
    --input_assembly assembly.fasta \
    --ont_raw_reads data/ONT_raw_reads.fastq.gz \
    --outdir results
```

This skips stages 01-02b. ONT reads are still required for purge_dups, BlobToolKit coverage, polishing and the mitogenome.

### Profiles

| Profile | Executor | Software |
|---------|----------|----------|
| `binac2` | SLURM (`--project` sets the account) | Apptainer, cached in `--apptainer_cache_dir` |
| `local` | local | Apptainer |
| `conda` | local | conda/mamba environments (for quick local iteration) |
| `test` | local | none (use with `-stub`) |

Per-process CPUs, memory and time are set in `process_resource.config`; adjust them there for another cluster. Failed tasks are retried with more resources only for out-of-memory and walltime exit codes (`--max_retries`, default 3).

---

## Parameter Reference

The full list, with comments, is in `nextflow.config`. The most used ones:

### Core inputs

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ont_raw_reads` | | ONT reads (FASTQ, gzipped) |
| `--input_assembly` | `null` | Existing assembly; skips k-mer sizing and assembly |
| `--outdir` | `results` | Output directory |
| `--species_name` | `Porites_sp` | Name used in output file names |
| `--symbiont_refs` | `null` | Symbiont reference genome(s) for read removal (comma-separated) |
| `--genome_size` | `null` | Genome size in bp; null = GenomeScope 2 estimate |
| `--busco_lineage` | `anthozoa_odb12` | BUSCO lineage for all assembly QC and BRAKER4 |

### Assembly and filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assemblers` | `['flye','hifiasm','raven']` | Any of `necat`, `flye`, `hifiasm`, `raven` |
| `--busco_tie_band` | `0.5` | BUSCO difference (percentage points) treated as a tie, decided by N50 |
| `--assembly_size_tolerance` | `0.30` | Deviation from the size estimate before an assembly is ranked last |
| `--purge_dups_transition` | `null` | calcuts haploid/diploid transition; null = automatic |
| `--taxon_id` | `46719` | Taxon ID of the BlobToolKit dataset |
| `--host_taxon_name` / `--symbiont_taxon_name` | `Cnidaria` / `Symbiodiniaceae` | BlobToolKit taxa kept as host and symbiont |
| `--auto_filter` | `true` | Automatic best-hit taxon filter |
| `--filter_json` | `null` | Manual selection exported from the BlobToolKit Viewer, used instead of the automatic filter |
| `--scaffold_reference` | `null` | Reference genome for RagTag; null = no scaffolding |
| `--medaka_model` | `null` | Medaka model; null = read from the ONT read headers |

### Annotation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--rnaseq_strandedness` | `reverse` | `unstranded`, `forward` or `reverse` |
| `--funannotate_busco_db` | `metazoa` | BUSCO set inside the funannotate database |
| `--related_proteomes` | `null` | Folder of related-species protein FASTA files; enables the gap-fusion split and the isoform filter |
| `--n_chromosome_scaffolds` | `null` | Number of chromosome-scale scaffolds (`scaffold_1..N` after sorting) |
| `--transcript_edits` | `null` | TSV of manual transcript edits (drop a transcript or replace its CDS) |
| `--rm_libraries` | `null` | RepeatMasker `Libraries/` folder with Dfam partitions for RepeatClassifier |

### Executor

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--project` | `null` | SLURM account |
| `--publish_mode` | `copy` | Publish mode for results |
| `--publish_mode_intermediate` | `symlink` | Publish mode for large intermediates; use `copy` if the work dir is purged |
| `--max_retries` | `3` | Retries for out-of-memory and walltime failures |
| `--email` / `--email_on_failure` | `null` / `false` | Completion e-mail |

---

## Multi-Assembler Comparison

```bash
nextflow run main.nf -profile binac2 \
    --ont_raw_reads data/ONT_raw_reads.fastq.gz \
    --assemblers 'necat,flye,raven' \
    --outdir results
```

1. Each assembler runs independently.
2. Every assembly gets BUSCO, gfastats and QUAST.
3. A comparison table (TSV and JSON) and a ranking table are written to `02_assembly/`.
4. The top-ranked assembly continues downstream; all assemblies are published.

| Assembler | Notes | Default memory |
|-----------|-------|----------------|
| `necat` | ONT reads; needs a local NECAT install | 256 GB |
| `flye` | Runs with `--meta` by default (`--flye_meta`) | 256 GB |
| `hifiasm` | PacBio HiFi or high-accuracy ONT | 256 GB |
| `raven` | Fast, lower memory | 128 GB |

---

## Stage Toggles and Dependencies

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--run_genomescope` | `true` | GenomeScope 2 size estimate (ignored with `--genome_size`) |
| `--skip_porechop` | `false` | Use the raw reads without adapter removal |
| `--skip_symbiont_filter` | `false` | Keep reads matching `--symbiont_refs` |
| `--skip_blobtools` | `false` | No BlobToolKit host filter |
| `--run_assembly_checkpoints` | `true` | BUSCO + gfastats after the host filter |
| `--run_screens` | `true` | NCBI contamination screens |
| `--run_polish_qc` | `true` | QV/BUSCO comparison before and after polishing |
| `--run_mitogenome` | `true` | Mitogenome (also needs `--mito_ref`) |
| `--run_repeats` | `true` | Repeat library and soft-masking |
| `--run_structural_annotation` | `true` | Gene prediction (needs `--run_repeats`, `--rnaseq_dir` and `--protein_refs`) |
| `--run_functional_annotation` | `true` | Functional annotation (needs structural annotation) |
| `--run_release` | `false` | NCBI submission and release folder (needs functional annotation) |

```
00 Preprocessing
 ├── 01 K-mer size estimate ── 02 Assembly ── 02b Selection      (or 02c with --input_assembly)
 │                                               │
 │                                    03b purge_dups ── 03 BlobToolKit ── 03c Screens ── 03d RagTag
 │                                                                                         │
 │                                                                                  04 Polishing
 │                                                                                         │
 ├── 05 Mitogenome (needs --mito_ref)                                              06 Repeats
 │                                                                                         │
 │                                                     07 Structural annotation (+07c refinement)
 │                                                                                         │
 │                                                                  08 Functional annotation
 │                                                                                         │
 │                                                                         10 Release (optional)
 └──────────────────────────────────────────────────────────────────────── 09 Summary + MultiQC
```

Structural annotation is skipped when `--rnaseq_dir` or `--protein_refs` does not exist. The run stops at start-up if structural annotation is on and repeats are off, because BRAKER4 needs the soft-masked assembly.

---

## Contamination Screens

With `--run_screens true` (default) the host-filtered contigs go through FCS-adaptor, VecScreen and FCS-GX before scaffolding. Required: `--univec`, `--fcs_gx_db`, `--fcs_adaptor_image`, `--fcs_gx_image` and `--fcs_tax_id` (must be in the GX database). The run checks these at start-up.

Clear-cut hits are applied automatically. Unclear hits (FCS-GX REVIEW, VecScreen Strong/Moderate) stop the run and are listed in `03c_screens/contigs/review/screen_review.tsv`. Decide them in a TSV (`seq`, `start`, `end`, `decision` = `exclude` | `mask` | `keep`, `reason`), pass it with `--screen_decisions`, and resume with `-resume`.

---

## NCBI Release

`--run_release true` builds a table2asn submission and a release folder. Required: `--table2asn`, `--ncbi_species` (exactly as in the BioSample), `--locus_tag`, `--ncbi_sbt`, `--ncbi_cmt`, `--bioproject`, `--biosample`, and `--ncbi_isolate` or `--ncbi_strain`. These are checked at start-up.

Useful options:

| Parameter | Description |
|-----------|-------------|
| `--release_name`, `--release_version` | Release folder name and version |
| `--release_notes` | Markdown file added to the release README as "Known limitations" |
| `--product_draft` | TSV of checked product names (symbol to product), used before the UniProt lookup |
| `--product_cache`, `--uniprot_offline` | Reuse a UniProt lookup cache; offline mode uses only the cache |
| `--symbiont_busco_lineage` | BUSCO lineage for the symbiont bin |
| `--ncbi_linkage_evidence` | Scaffold-gap linkage evidence (default `align-genus` with RagTag) |
| `--qc_expected_lineage`, `--qc_expected_division`, `--qc_chrom_prefix` | Expectations for the release QC checks |

Before the folder is written, the release QC checks validate the models, protein identity, the protein-coding count across all reports, scaffold names, repeat labels, the screens and the table2asn reports. Any failed check stops the run; the results are in `qc/annotation_qc.tsv`.

---

## Output Structure

```
outdir/
├── 00_preprocessing/
│   ├── porechop/                    # adapter-trimmed reads
│   └── reads/                       # assembly, polishing and k-mer read sets, read QC
├── 01_kmer_genome_size/
│   ├── meryl/                       # k-mer database and histogram
│   └── genomescope/                 # genome size and heterozygosity model
├── 02_assembly/
│   ├── flye/ hifiasm/ raven/ necat/ # one folder per assembler
│   ├── qc/{busco,gfastats,quast}/   # per-assembly QC
│   ├── comparison/                  # comparison table
│   └── selection/                   # ranking and the selected assembly
├── 03b_purge_dups/                  # purged assembly, haplotigs, coverage cutoffs
├── 03_filtering/
│   ├── blobtools/                   # BlobDir (open with the BlobToolKit Viewer)
│   ├── filtered/                    # host and symbiont assemblies
│   └── plots/                       # blob, snail and cumulative plots
├── 03c_screens/contigs/             # FCS-adaptor, VecScreen, FCS-GX, review and actions
├── 03d_scaffolding/                 # RagTag scaffolds and AGP
├── 04_polishing/
│   ├── medaka/                      # polished assembly
│   ├── funannotate/                 # cleaned and sorted assembly
│   └── qc/                          # BUSCO, Merqury QV, polishing decision
├── 05_mitogenome/
│   ├── mitohifi/ getorganelle/ ont_polish/
│   └── final/                       # mitogenome FASTA and annotation
├── 06_repeats/
│   ├── repeatmodeler/ edta/ combined_library/ classify/
│   ├── repeatmasker/                # soft-masked assembly, .out, .tbl
│   └── barrnap/                     # rRNA predictions
├── 07_structural_annotation/
│   ├── rnaseq/                      # trimmed reads, STAR index and alignments
│   ├── braker/                      # BRAKER4 output
│   ├── trnascan/ ncrna/ stringtie/
│   ├── merged_annotation/           # merged and validated GFF3
│   ├── refine/                      # refined GFF3, reports, haplotig and gap-fusion scans
│   └── qc/{compleasm,omark}/
├── 08_functional_annotation/
│   ├── interproscan/ eggnog/ phobius/
│   └── funannotate/                 # annotated GFF3, proteins, GenBank, stats
├── 10_release/                      # intermediate release files (if --run_release)
├── release/<release_name>/          # final release folder (if --run_release)
├── summary/
│   ├── multiqc/multiqc_report.html  # integrated report
│   └── ...                          # summary tables and software versions
└── pipeline_info/                   # Nextflow timeline, report, trace and DAG
```

### Release folder

| Path | Content |
|------|---------|
| `README.md`, `summary.tsv` | Release summary, methods with tool versions, known limitations |
| `genome/` | Soft-masked, screened scaffolds; `scaffold_name_map.tsv` |
| `annotation/` | GFF3, proteins, CDS and mRNA transcripts, `annotations.txt`, `.gbk` (all consistent with the submission) |
| `annotation/extra/` | High-confidence GFF3 and proteins, `splits.tsv`, `id_map.tsv`, haplotig and gap candidate lists, `rrna_decisions.tsv` |
| `repeats/` | Repeat GFF3, RepeatMasker `.tbl` and `.out`, classified library |
| `ncbi/` | table2asn `.sqn`, `.val`, `.dr`, `.gbf` (the mitogenome is in the same `.sqn`) |
| `mitogenome/` | Mitogenome sequence and annotation |
| `qc/` | `annotation_qc.tsv`, refinement and protein-identity reports, product curation, screening reports, BUSCO, symbiont BUSCO, Merqury QV |

---

## License

See the [LICENSE](LICENSE) file for details.
