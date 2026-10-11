# NGS Automated Variant Calling Pipeline
### (Nextflow DSL2)

**Author:** Maryam Sattar
*Bioinformatics | Computational Biology*

---

## 1. Project Overview

A reproducible Nextflow DSL2 pipeline for paired-end human whole-genome sequencing (WGS) data. It downloads reads from the NCBI Sequence Read Archive (SRA), trims and quality-checks them with `fastp`, aligns them to the human GRCh38/hg38 reference with `BWA-MEM`, sorts and indexes the BAM with `SAMtools`, and calls variants with `BCFtools`. The final output is a compressed, indexed VCF of SNVs and indels.

The example run uses the paired-end SRA accession `SRR062634`.

## 2. Pipeline Architecture

```
SRA accession list (samples.csv)
            |
            v
   +------------------+
   |   FETCH_FASTQ    |  prefetch + fasterq-dump
   +--------+---------+
            |
            v
   +------------------+
   |    FASTP_QC      |  trimming + quality filtering + HTML/JSON reports
   +--------+---------+
            |
            v
   +------------------+
   |  ALIGN_AND_SORT  |  bwa mem + samtools sort + samtools index
   +--------+---------+
            |
            v
   +------------------+
   | VARIANT_CALLING  |  bcftools mpileup + call + index
   +--------+---------+
            |
            v
     Compressed VCF
     (.vcf.gz + .csi)
```

## 3. Workflow Processes

| Process | Tools | What it does | Output |
|---|---|---|---|
| `FETCH_FASTQ` | `prefetch`, `fasterq-dump` | Downloads the SRA run and converts it to paired FASTQ (`--split-files`). Fails if either FASTQ is empty. | `results/raw_fastq/` |
| `FASTP_QC` | `fastp` | Removes adapters, filters low-quality reads, writes QC reports. | `results/clean_fastq/`, `results/qc/` |
| `ALIGN_AND_SORT` | `bwa mem`, `samtools` | Maps cleaned reads to hg38, sorts by coordinate, indexes the BAM. | `results/bam/` |
| `VARIANT_CALLING` | `bcftools` | `mpileup` then `call -mv` (multiallelic caller, variant sites only), compressed with `-Oz` and indexed. | `results/vcf/` |

> The pipeline ends at variant calling. It does not yet filter variants, annotate them, or interpret them biologically.

## 4. Requirements

- Linux or WSL2 (Ubuntu 20.04/22.04)
- Conda/Miniconda installed at `~/miniconda` (the `safe`, `conda` and `parallel` profiles source `$HOME/miniconda/etc/profile.d/conda.sh`)
- Java 17 or later for Nextflow
- Enough disk space for the reference genome, FASTQ, cleaned FASTQ, BAM and Nextflow's `work/` folder. A human WGS run needs many gigabytes.
- Suggested resources per task: download 2 CPUs / 4 GB, QC 2 CPUs / 4 GB, alignment 4 CPUs / 12 GB, variant calling 2 CPUs / 8 GB

### Conda environments

```bash
# Nextflow environment
conda create -n nextflow_env -y openjdk=17 nextflow

# NGS tools environment (the profiles activate this one)
conda create -n ngs_tools -y -c conda-forge -c bioconda sra-tools fastp bwa samtools bcftools
```

## 5. Project Structure

```
NGS_PROJECT/
├── main.nf              # workflow and processes
├── nextflow.config      # profiles and resource limits
├── samples.csv          # SRA accessions, one per line
├── environment.yml      # tool list for Conda (reproducibility)
├── Dockerfile           # optional Docker image
├── README.md
├── ref/
│   ├── hg38.fasta
│   ├── hg38.fasta.fai
│   └── hg38.fasta.amb / .ann / .bwt / .pac / .sa   # BWA index
├── work/                # Nextflow task files (keep for -resume)
└── results/
    ├── raw_fastq/
    ├── clean_fastq/
    ├── qc/
    ├── bam/
    └── vcf/
```

## 6. Reference Genome Setup

The pipeline expects `ref/hg38.fasta` and its index files next to it. It does **not** download or rebuild the reference. If you already have `hg38.fasta` with `.amb`, `.ann`, `.bwt`, `.fai`, `.pac` and `.sa`, keep them as they are. Otherwise build them once:

```bash
cd ~/NGS_PROJECT/ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O hg38.fasta.gz
gunzip hg38.fasta.gz

samtools faidx hg38.fasta
bwa index hg38.fasta
```

The reference file must be named `hg38.fasta`, because `main.nf` looks for that exact name.

## 7. Input Samples

`samples.csv` holds one SRA run accession per line, with no header:

```
SRR062634
```

Use accessions that are paired-end human DNA sequencing data.

## 8. Running the Pipeline

```bash
cd ~/NGS_PROJECT
conda activate nextflow_env
nextflow run main.nf -profile safe
```

If a run stops, fix the reported problem and continue where it stopped:

```bash
nextflow run main.nf -profile safe -resume
```

This is a real run. The processes have no `stub:` blocks, so `-stub-run` is not supported.

### Profiles

| Profile | Use it when |
|---|---|
| `safe` | Recommended. Uses the installed `ngs_tools` Conda environment, one task at a time. |
| `conda` | Same as `safe` (kept for clarity). |
| `parallel` | You have spare RAM and disk: up to two tasks at a time. |
| `docker` | You want containers. Docker must be installed and running. Build with `docker build -t ngs-wgs:1.0 .`, then run `nextflow run main.nf -profile docker -resume`. |

You can also override the defaults:

```bash
nextflow run main.nf -profile safe --samples samples.csv --ref ref/hg38.fasta --outdir results
```

## 9. Output Directory

```
results/
├── raw_fastq/       # FASTQ files from SRA
├── clean_fastq/     # fastp-processed FASTQ files
├── qc/              # fastp HTML and JSON reports
├── bam/             # sorted BAM files and indexes
└── vcf/             # compressed VCF files and .csi indexes
```

Outputs are hard-linked into `results/` instead of copied, which saves disk space. Keep `work/`, because `-resume` depends on it.

## 10. Reproducibility and Resource Management

- Nextflow DSL2 splits the workflow into separate processes, each tagged by sample.
- The `safe` profile runs one task at a time, which suits a laptop or WSL.
- Resource limits are set per process label (`download`, `qc`, `alignment`, `variants`) in `nextflow.config`.
- `environment.yml` and the `Dockerfile` record the tool list for reproducibility. The recommended `safe` run does not need Docker or a new Conda environment.
- Do not commit the reference genome, FASTQ, BAM, `work/` or `results/` to GitHub.

## 11. Current Scope and Future Extensions

The current version ends at variant calling. Possible next steps:

- Variant quality filtering
- Annotation with VEP, SnpEff or ANNOVAR
- ClinVar and gnomAD annotation
- Functional prediction (SIFT, PolyPhen, CADD, REVEL)
- Downstream statistical or machine-learning analysis of annotated variants

**Maryam Sattar**
*Bioinformatics | Computational Biology*

