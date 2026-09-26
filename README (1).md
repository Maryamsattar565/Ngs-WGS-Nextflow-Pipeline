# NGS Automated Variant Calling Pipeline
### (Nextflow DSL2)

**Author:** Maryam Sattar
*Bioinformatics | Computational Biology*

---

## 1. Project Overview

An automated and reproducible Nextflow DSL2 pipeline for processing paired-end human whole-genome sequencing (WGS) or targeted resequencing data. The workflow retrieves sequencing data from the NCBI Sequence Read Archive (SRA), performs adapter trimming and quality filtering with `fastp`, aligns reads to the human GRCh38/hg38 reference genome using `BWA-MEM`, sorts and indexes BAM files with `SAMtools`, and performs variant calling using `BCFtools`. The final output is a compressed and indexed VCF file containing detected SNVs and indels.

## 2. Pipeline Architecture

```
SRA Accession List (samples.csv)
            |
            v
      +------------------+
      |   FETCH_FASTQ    |  <-- prefetch + fasterq-dump
      +--------+---------+
               |
               v
      +------------------+
      |    FASTP_QC      |  <-- trimming + quality filtering + reports
      +--------+---------+
               |
               v
      +------------------+
      | ALIGN_AND_SORT   |  <-- BWA-MEM + SAMtools
      +--------+---------+
               |
               v
      +------------------+
      | VARIANT_CALLING  |  <-- BCFtools mpileup + call
      +--------+---------+
               |
               v
        Compressed VCF
        (.vcf.gz + .csi)
```

## 3. Workflow Processes

### 3.1 FETCH_FASTQ
Downloads paired-end sequencing reads from the NCBI SRA using `sra-tools`.
- `prefetch` retrieves the SRA run.
- `fasterq-dump` converts the run into FASTQ files using `--split-files`.
- The process outputs paired FASTQ files for downstream analysis.

### 3.2 FASTP_QC
Performs read preprocessing and generates quality-control reports using `fastp`.
- Removes adapter sequences.
- Filters low-quality bases/reads according to fastp settings.
- Produces cleaned paired FASTQ files.
- Produces HTML and JSON QC reports.

### 3.3 ALIGN_AND_SORT
Maps cleaned reads to the human GRCh38/hg38 reference genome and prepares sorted BAM files.
- `BWA-MEM` performs read alignment.
- `SAMtools sort` creates coordinate-sorted BAM files.
- `SAMtools index` creates BAM index files.

### 3.4 VARIANT_CALLING
Identifies candidate genomic variants from the aligned BAM file using `BCFtools`.
- `bcftools mpileup` generates genotype likelihood information from the alignment.
- `bcftools call` performs variant calling with the multiallelic caller (`-m`) and outputs variant sites (`-v`).
- The VCF is compressed with bgzip-compatible BCFtools output (`-Oz`).
- The final VCF is indexed with `bcftools index`.

> The current pipeline performs variant calling but does not yet include downstream variant filtering, functional annotation, ClinVar annotation, or biological interpretation.

## 4. Environmental Setup

### 4.1 Prerequisites
- Operating system: Linux or WSL2 (Ubuntu 20.04/22.04)
- Conda/Miniconda or Mamba
- Java 17 or later for Nextflow
- Sufficient disk space for the reference genome, SRA/FASTQ data, BAM files, and VCF outputs

### 4.2 Conda Environments

```bash
# Nextflow environment
conda create -n nextflow_env -y openjdk=17 nextflow

# NGS tools environment
conda create -n ngs_tools -y sra-tools fastp bwa samtools bcftools
```

Nextflow can use Conda environments to provide isolated software dependencies for processes.

## 5. Project Structure

```
Ngs_project/
├── main.nf
├── nextflow.config
├── samples.csv
├── .gitignore
├── README.md
├── ref/
│   ├── hg38.fasta
│   ├── hg38.fasta.fai
│   └── BWA index files
└── results/
    ├── raw_fastq/
    ├── clean_fastq/
    ├── qc/
    ├── bam/
    └── vcf/
```

## 6. Reference Genome Setup

The pipeline uses the UCSC GRCh38/hg38 reference FASTA. The FASTA is indexed for SAMtools and BWA before pipeline execution.

```bash
cd ~/Ngs_project/ref

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O hg38.fasta.gz
gunzip hg38.fasta.gz

samtools faidx hg38.fasta
bwa index hg38.fasta
```

## 7. Input Samples

SRA run accessions are supplied through `samples.csv`. The current example contains one accession:

```
SRR062634
```

The selected accession should represent suitable paired-end human DNA sequencing data for the intended WGS or targeted resequencing analysis.

## 8. Running the Pipeline

### 8.1 Stub / Validation Run
The stub run tests the workflow structure without performing the full data-processing workload.

```bash
nextflow run main.nf -profile safe -stub-run
```

### 8.2 Production Run
Runs the complete workflow with one task at a time, which is suitable for limited WSL or laptop resources.

```bash
conda activate nextflow_env
nextflow run main.nf -profile safe
```

### 8.3 Parallel Run
Runs with a maximum of two concurrent tasks when additional CPU and memory resources are available.

```bash
nextflow run main.nf -profile parallel
```

## 9. Output Directory

```
results/
├── raw_fastq/       # Original FASTQ files retrieved from SRA
├── clean_fastq/     # fastp-processed FASTQ files
├── qc/              # fastp HTML and JSON reports
├── bam/             # Sorted and indexed BAM files
└── vcf/             # Compressed and indexed VCF files
```

## 10. Reproducibility and Resource Management

- Nextflow DSL2 separates the workflow into modular processes.
- Conda environments provide isolated software dependencies.
- The `safe` profile limits local execution to one concurrent task.
- The `parallel` profile allows up to two concurrent tasks.
- Nextflow work directories and generated outputs can be excluded from version control.
- Large genomic data and reference files should not be committed to GitHub.
- Nextflow is designed for reproducible and portable scientific workflows, and supports execution with software environments such as Conda.

## 11. Current Scope and Future Extensions

The current version ends at variant calling. Possible future extensions include:

- Variant quality filtering.
- Variant annotation with VEP, SnpEff, or ANNOVAR.
- Integration of ClinVar annotations.
- Population-frequency annotation such as gnomAD.
- Functional prediction features such as SIFT, PolyPhen, CADD, or REVEL.
- Downstream statistical or machine-learning analysis of annotated variants.

## 12. Important Execution Note

The pipeline includes a stub-mode validation command. A successful production run requires successful SRA download, sufficient storage, an appropriate paired-end dataset, and successful completion of every process. The project should only claim a completed VCF analysis after a real sample has successfully passed through the entire workflow.

## 13. Author

**Maryam Sattar**
*Bioinformatics | Computational Biology*

## 14. References

- [Nextflow documentation and training materials](https://nextflow.io/)
- [BCFtools documentation](https://samtools.github.io/bcftools/)
