# GenomeMaker - Coral Genome Assembly & Annotation Pipeline

[![Alpha Release](https://img.shields.io/badge/release-alpha-orange)](https://github.com/your-username/your-repo/releases)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A523.04.0-%233099C1)](https://www.nextflow.io/)
[![DSL2](https://img.shields.io/badge/DSL2-Compatible-brightgreen)](https://www.nextflow.io/docs/latest/dsl2.html)

[Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) pipeline for *de novo* assembly and annotation of coral (Scleractinia) genomes from Oxford Nanopore Technologies (ONT) long reads.

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Installation Requirements](#installation-requirements)
- [Quick Start](#quick-start)
- [Parameter Reference](#parameter-reference)
- [Multi-Assembler Comparison](#multi-assembler-comparison)
- [Optional Stage Toggles](#optional-stage-toggles)
- [Output Structure](#output-structure)
<!-- - [Citation](#citation) -->

---

## Pipeline Overview

This pipeline automates the complete workflow for assembling and annotating coral genomes from ONT long-read sequencing data. The workflow consists of nine major stages:

| Stage | Module | Description |
|-------|--------|-------------|
| **00** | `00_preprocess_reads` | ONT read preprocessing: Porechop adapter removal, FastQC/MultiQC/NanoPlot QC, chopper read splitting (assembly vs. polishing reads), symbiont removal |
| **01** | `01_kmer_genome_size` | K-mer profiling with Meryl, histogram generation, and GenomeScope 2.0 genome size/coverage estimation |
| **02** | `02_multi_assembler` | Parallel *de novo* assembly with configurable assemblers (NECAT, Flye, hifiasm, Raven) with independent BUSCO, gfastats, and QUAST QC per assembler |
| **03** | `03_assembly_filtering` | Assembly filtering via BlobTools to remove contamination |
| **04** | `04_assembly_polish` | Polishing with Racon and Medaka, followed by funannotate clean/sort and final BUSCO assessment |
| **05** | `05_mitogenome` | Mitochondrial genome extraction: Canu assembly, Circlator circularization, Racon polishing |
| **06** | `06_repeats` | Repeat annotation with EDTA and RepeatModeler, soft-masking with RepeatMasker |
| **07** | `07_structural_annotation` | Gene prediction with BRAKER3 using RNA-Seq and protein evidence, tRNA prediction with tRNAscan-SE, UTR annotation, AGAT validation |
| **08** | `08_functional_annotation` | Functional annotation via InterProScan, EggNOG-mapper, Phobius, and funannotate synthesis |
| **09** | `09_summary` | Aggregated summary report combining all QC results |

---

## Installation Requirements

| Requirement | Version | Notes |
|-------------|---------|-------|
| [Nextflow](https://www.nextflow.io/) | >= 23.04.0 | Required for DSL2 support |
| [Apptainer](https://apptainer.org/) / Singularity | >= 1.1.0 | Recommended for production (binac2 profile) |
| [Docker](https://www.docker.com/) or [Podman](https://podman.io/) | Latest | Alternative for local profile |
| [Conda](https://docs.conda.io/) / Mamba | >= 22.0 | Fallback for tools without containers |
| [Git](https://git-scm.com/) | Latest | For cloning the repository |

### Getting Nextflow

```bash
# Install via Conda (recommended)
conda install -c bioconda nextflow

# Or download directly
curl -s https://get.nextflow.io | bash
```

Verify installation:

```bash
nextflow -v
# nextflow version 23.04.0 ...
```

---

## Container Setup

The pipeline uses containerized processes for reproducibility. Two container strategies are supported:

### Apptainer/Singularity (Recommended for HPC)

The **binac2** profile uses Apptainer by default. Images are auto-pulled on first run and cached to disk:

```bash
# Configure cache directory (set in nextflow.config params)
# Default: /pfs/10/project/apptainer_cache/${USER}
export APPTAINER_CACHEDIR=/path/to/apptainer_cache
```

**Important:** Some custom tools (e.g., **NECAT**) may require custom container builds if the Biocontainers image is unavailable or outdated. In such cases, build a local Apptainer definition:

```bash
# Example: Build a custom NECAT container
sudo apptainer build necat.sif \
    --bind /pfs/10/project/apptainer_cache:/apptainer_cache \
    apptainer/necat.def
```


---

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/your-org/ID-por-genome.git
cd ID-por-genome
```

### 2. Prepare Input Data

Place your ONT raw reads in a known location:

```bash
# Example directory structure
data/
└── ONT_raw_reads.fastq.gz    # Compressed ONT FASTQ
└── rnaseq/                   # Optional: RNA-Seq directory for structural annotation
```

### 3. Run the Pipeline

#### Default Single-Assembler Run (Local)

```bash
nextflow run main.nf -profile local \
    --ont_raw_reads data/ONT_raw_reads.fastq.gz \
    --outdir results \
    --species_name "coral_sp"
```



#### Stub Validation Run (No Computation)

```bash
nextflow run main.nf -profile test
```

This runs all processes in stub mode (echo only) to validate the workflow graph and input paths without executing any tools.

---

## Parameter Reference

### Core Inputs

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--ont_raw_reads` | path | `data/ONT_raw_reads.fastq.gz` | Path to compressed ONT FASTQ reads |
| `--symbiont_refs` | path(s) | `null` | Comma-separated paths to Symbiodiniaceae reference genomes (null = auto-download) |
| `--symbiont_taxids` | string | `'58267'` | NCBI taxon ID for Symbiodiniaceae |
| `--coverage_refs` | path(s) | `null` | Comma-separated paths to Porites reference genomes for coverage assessment |
| `--coverage_taxid` | int | `46719` | NCBI taxon ID for Porites |
| `--mito_ref` | path | `null` | Path to reference mitogenome (null = auto-download from taxon ID) |
| `--mito_taxid` | int | `31689` | NCBI taxon ID for coral mitochondrial reference |
| `--rnaseq_dir` | path | `data/rnaseq` | Directory containing RNA-Seq data for structural annotation |
| `--protein_refs` | path | `refs/Metazoa.fa` | Path to Metazoa protein references for gene prediction |
| `--genomescope_script` | path | `scripts/genomescope.R` | Path to GenomeScope 2.0 R script |
| `--outdir` | path | (required) | Output directory for all pipeline results |
| `--species_name` | string | `'coral_sp'` | Species name used in sample metadata and output paths |

### Assembly Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--assemblers` | list | `['necat']` | Comma-separated list of assemblers: `necat`, `flye`, `hifiasm`, `raven` |
| `--necat_config` | path | `config/necat_config.template` | NECAT configuration template |
| `--necat_read_list` | path | `config/necat_read_list.txt` | NECAT read list file |
| `--genome_size` | int | `null` | Expected genome size in bp (null = auto-detected from k-mer analysis) |
| `--kmer_size` | int | `21` | K-mer size for genome size estimation |
| `--min_coverage` | int | `2` | Minimum coverage threshold |

### Annotation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--busco_lineage_euk` | string | `'eukaryota_odb10'` | BUSCO eukaryotic lineage dataset |
| `--busco_lineage_meta` | string | `'metazoa_odb10'` | BUSCO metazoan lineage dataset |
| `--rnaseq_strandedness` | string | `'unstranded'` | RNA-Seq strandedness: `unstranded`, `forward`, `reverse`, `auto` |
| `--funannotate_db` | path | `funannotate_DB` | Path to funannotate database |

### Stage Toggles

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_genomescope` | bool | `true` | Run k-mer genome size estimation |
| `--run_mitogenome` | bool | `true` | Extract and assemble mitochondrial genome |
| `--run_repeats` | bool | `true` | Run repeat annotation (EDTA + RepeatModeler + RepeatMasker) |
| `--run_structural_annotation` | bool | `true` | Run gene prediction (BRAKER3) |
| `--run_functional_annotation` | bool | `true` | Run functional annotation (InterProScan, EggNOG, funannotate) |
| `--skip_symbiont_filter` | bool | `false` | Skip symbiont read removal during preprocessing |
| `--skip_blobtools` | bool | `false` | Skip BlobTools assembly filtering |

### Resource Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--threads_default` | int | `16` | Default CPU threads |
| `--threads_busco` | int | `16` | BUSCO threads |
| `--threads_star` | int | `32` | STAR alignment threads |
| `--threads_braker` | int | `32` | BRAKER3 threads |
| `--threads_funannotate` | int | `32` | funannotate threads |
| `--resource_necat_memory` | string | `'256 GB'` | NECAT memory limit |
| `--resource_necat_cpus` | int | `32` | NECAT CPU limit |
| `--resource_necat_time` | string | `'48h'` | NECAT time limit |
| `--resource_flye_memory` | string | `'128 GB'` | Flye memory limit |
| `--resource_flye_cpus` | int | `16` | Flye CPU limit |
| `--resource_flye_time` | string | `'24h'` | Flye time limit |
| `--resource_hifiasm_memory` | string | `'256 GB'` | hifiasm memory limit |
| `--resource_hifiasm_cpus` | int | `32` | hifiasm CPU limit |
| `--resource_hifiasm_time` | string | `'48h'` | hifiasm time limit |
| `--resource_raven_memory` | string | `'64 GB'` | Raven memory limit |
| `--resource_raven_cpus` | int | `16` | Raven CPU limit |
| `--resource_raven_time` | string | `'12h'` | Raven time limit |

### Executor Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--project` | string | `null` | SLURM account/project identifier |
| `--publish_mode` | string | `'copy'` | Nextflow publish mode: `copy`, `symlink`, `rellink`, `move` |
| `--enable_resume` | bool | `true` | Enable pipeline resume from last successful checkpoint |
| `--max_retries` | int | `3` | Maximum retry attempts per failed process |


---

## Multi-Assembler Comparison

The pipeline supports running multiple assemblers in parallel for direct comparison. This is useful for benchmarking assembly strategies or selecting the best assembler for a given dataset.

```bash
# Run NECAT, Flye, and Raven in parallel
nextflow run main.nf -profile binac2 \
    --ont_raw_reads data/ONT_raw_reads.fastq.gz \
    --assemblers 'necat,flye,raven' \
    --outdir results/multi_assembly
```

**How it works:**

1. Each assembler runs independently and in parallel using the specified resources
2. Every assembly receives independent QC: BUSCO completeness, gfastats statistics, and QUAST reports
3. A **comparison report** aggregates results across all assemblers (TSV, HTML, and JSON formats)
4. The **primary assembly** is chosen autoatically and carried downstream. All assembly are exported to output 

**Recommended assemblers for ONT data:**

| Assembler | Best For | Min. Memory |
|-----------|----------|-------------|
| `necat` | High-quality ONT long reads | 256 GB |
| `flye` | General-purpose ONT assembly | 128 GB |
| `raven` | Fast ONT assembly (lower memory) | 64 GB |
| `hifiasm` | PacBio HiFi or high-accuracy ONT | 256 GB |

---

## Optional Stage Toggles

Fine-tune the pipeline by enabling or disabling individual stages. This is useful for quick runs, troubleshooting, or when certain inputs (e.g., RNA-Seq) are unavailable.

### Skipping Optional Stages

```bash
# Minimal pipeline: assembly only (no mitogenome, repeats, or annotation)
nextflow run main.nf -profile local \
    --ont_raw_reads data/ONT_raw_reads.fastq.gz \
    --outdir results/minimal \
    --run_mitogenome false \
    --run_repeats false \
    --run_structural_annotation false \
    --run_functional_annotation false
```

### Stage Dependencies

```
Stage 00 (Preprocessing)
    └── Stage 01 (K-mer estimation)
            └── Stage 02 (Assembly)
                    └── Stage 03 (Filtering)
                            └── Stage 04 (Polishing)
                                    ├── Stage 05 (Mitogenome) ── requires --mito_ref
                                    ├── Stage 06 (Repeats)
                                    ├── Stage 07 (Structural annotation) ── requires --rnaseq_dir and --protein_refs
                                    │       └── Stage 08 (Functional annotation) ── requires Stage 07 output
                                    └── Stage 09 (Summary)
```

**Note:** Structural and functional annotation require both `--rnaseq_dir` and `--protein_refs` to be set. If either is missing, those stages are automatically skipped regardless of toggle settings.

---

## Output Structure

All outputs are written to the directory specified by `--outdir`. The structure is organized by stage:

```
outdir/
├── 00_reads_preprocessing/
│   ├── porechop/
│   │   └── ONT_reads_pc.fastq.gz          # Adapter-trimmed reads
│   ├── read_qc/
│   │   ├── fastqc/                         # FastQC reports
│   │   ├── nanoplot/                       # NanoPlot QC
│   │   └── multiqc/                        # MultiQC summary
│   ├── chopper/
│   │   ├── assembly_reads.fastq.gz         # Long reads for assembly
│   │   └── polishing_reads.fastq.gz        # Shorter reads for polishing
│   └── symbiont_removal/                   # Symbiont-filtered reads (if enabled)
│
├── 01_kmer_genome_size/
│   ├── meryl/
│   │   └── kmer_histogram.tabular          # K-mer frequency histogram
│   └── genomescope/
│       └── genomescope_output/             # Genome size & coverage estimates
│
├── assemblies/
│   ├── necat/
│   │   └── sample_assembly_necat.fasta     # NECAT assembly
│   ├── flye/
│   │   └── sample_assembly_flye.fasta      # Flye assembly
│   ├── comparison/
│   │   ├── assembly_comparison.tsv         # Multi-assembler comparison (TSV)
│   │   ├── assembly_comparison.html        # Multi-assembler comparison (HTML)
│   │   └── assembly_comparison.json        # Multi-assembler comparison (JSON)
│   └── qc/
│       ├── gfastats/                       # Per-assembler assembly statistics
│       ├── busco/                          # Per-assembler BUSCO results
│       └── quast/                          # Per-assembler QUAST reports
│
├── 03_assembly_filtering/
│   └── filtered_assembly.fasta             # BlobTools-filtered assembly
│
├── 04_assembly_polishing/
│   ├── racon/                              # Racon polished assembly
│   ├── medaka/                             # Medaka polished assembly
│   ├── funannotate_clean/                  # Cleaned & sorted assembly
│   └── busco/                              # Final BUSCO assessment
│
├── 05_mitogenome/                          # (if --run_mitogenome)
│   ├── canu/                               # Canu mitogenome assembly
│   ├── circlator/                          # Circularized mitogenome
│   ├── mito_bam.bam                        # Coverage alignment
│   ├── mito_bam.bai                        # BAM index
│   ├── coverage/                           # Coverage statistics
│   └── flagstat/                           # SAMtools flagstat
│
├── 06_repeats/                             # (if --run_repeats)
│   ├── edt/                                # EDTA repeat library
│   ├── repeatmodeler/                      # RepeatModeler repeat library
│   └── repeatmasker/                       # Soft-masked assembly
│
├── 07_structural_annotation/               # (if --run_structural_annotation)
│   ├── trnascan/                           # tRNA predictions
│   ├── star/                               # STAR alignments (RNA-Seq)
│   ├── braker3/                            # BRAKER3 gene predictions (GFF3)
│   └── proteins.faa                        # Predicted protein sequences
│
├── 08_functional_annotation/               # (if --run_functional_annotation)
│   ├── interproscan/                       # InterProScan annotations
│   ├── eggnog/                             # EggNOG annotations
│   ├── phobius/                            # Signal peptide predictions
│   └── funannotate/                        # Synthesized annotation directory
│       ├── fun_annotate.gff3               # Final GFF3 annotation
│       ├── fun_annotate.faa                # Final protein FASTA
│       └── ...                             # Additional annotation files
│
├── 09_summary/
│   └── summary_report.*                    # Aggregated summary across all stages
│
└── work/                                   # Nextflow work directory (intermediate files)
```

---


## License

See the [LICENSE](LICENSE) file for details.

---

<!-- **Developed for** [Project Name] — Coral Genome Assembly Initiative -->
