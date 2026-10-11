nextflow.enable.dsl=2

params.samples = params.samples ?: "${projectDir}/samples.csv"
params.ref     = params.ref     ?: "${projectDir}/ref/hg38.fasta"
params.outdir  = params.outdir  ?: "${projectDir}/results"

process FETCH_FASTQ {
    tag "$sra_id"
    label 'download'
    cpus 2
    publishDir "${params.outdir}/raw_fastq", mode: 'link', pattern: '*.fastq'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}_1.fastq"), path("${sra_id}_2.fastq"), emit: reads

    script:
    """
    mkdir -p sra_tmp
    prefetch --max-size 100G --output-directory . ${sra_id}

    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --outdir . \\
        --temp ./sra_tmp \\
        ${sra_id}/${sra_id}.sra

    test -s ${sra_id}_1.fastq
    test -s ${sra_id}_2.fastq
    rm -rf ${sra_id} sra_tmp
    """
}

process FASTP_QC {
    tag "$sra_id"
    label 'qc'
    cpus 2
    publishDir "${params.outdir}/clean_fastq", mode: 'link', pattern: '*_clean.fastq'
    publishDir "${params.outdir}/qc", mode: 'link', pattern: '*.html'
    publishDir "${params.outdir}/qc", mode: 'link', pattern: '*.json'

    input:
    tuple val(sra_id), path(r1), path(r2)

    output:
    tuple val(sra_id), path("${sra_id}_1_clean.fastq"), path("${sra_id}_2_clean.fastq"), emit: reads
    path "${sra_id}_fastp.html"
    path "${sra_id}_fastp.json"

    script:
    """
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${sra_id}_1_clean.fastq \\
        --out2 ${sra_id}_2_clean.fastq \\
        --html ${sra_id}_fastp.html \\
        --json ${sra_id}_fastp.json \\
        --thread ${task.cpus}
    """
}

process ALIGN_AND_SORT {
    tag "$sra_id"
    label 'alignment'
    cpus 4
    publishDir "${params.outdir}/bam", mode: 'link', pattern: '*.bam*'

    input:
    tuple val(sra_id), path(r1), path(r2)
    path ref_files

    output:
    tuple val(sra_id), path("${sra_id}.sorted.bam"), path("${sra_id}.sorted.bam.bai"), emit: bam

    script:
    def ref_fasta = ref_files.find { it.name == 'hg38.fasta' }
    """
    set -o pipefail
    bwa mem -t ${task.cpus} ${ref_fasta} ${r1} ${r2} \\
        | samtools sort -@ ${task.cpus} -m 1G -T ${sra_id}.sorttmp -o ${sra_id}.sorted.bam -

    samtools index ${sra_id}.sorted.bam
    """
}

process VARIANT_CALLING {
    tag "$sra_id"
    label 'variants'
    cpus 2
    publishDir "${params.outdir}/vcf", mode: 'link', pattern: '*.vcf.gz*'

    input:
    tuple val(sra_id), path(bam), path(bai)
    path ref_files

    output:
    path "${sra_id}.vcf.gz"
    path "${sra_id}.vcf.gz.csi"

    script:
    def ref_fasta = ref_files.find { it.name == 'hg38.fasta' }
    """
    set -o pipefail
    bcftools mpileup -f ${ref_fasta} -Ou ${bam} \\
        | bcftools call -mv -Oz -o ${sra_id}.vcf.gz

    bcftools index -f ${sra_id}.vcf.gz
    """
}

workflow {
    samples_ch = Channel
        .fromPath(params.samples)
        .splitText()
        .map { it.trim() }
        .filter { it }

    // Stage the FASTA and every matching existing index file together.
    ref_files_ch = Channel.fromPath("${params.ref}*").collect()

    FETCH_FASTQ(samples_ch)
    FASTP_QC(FETCH_FASTQ.out.reads)
    ALIGN_AND_SORT(FASTP_QC.out.reads, ref_files_ch)
    VARIANT_CALLING(ALIGN_AND_SORT.out.bam, ref_files_ch)
}
