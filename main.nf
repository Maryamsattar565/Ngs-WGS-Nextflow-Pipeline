nextflow.enable.dsl=2

params.samples = "${projectDir}/samples.csv"
params.ref     = "${projectDir}/ref/hg38.fasta"
params.outdir  = "${projectDir}/results"

process FETCH_FASTQ {
    tag "$sra_id"
    conda "${System.getenv('HOME')}/miniconda/envs/ngs_tools"
    cpus 4
    publishDir "${params.outdir}/raw_fastq",
        mode: 'copy',
        pattern: "*.fastq"

    input:
    val sra_id

    output:
    tuple val(sra_id),
          path("${sra_id}_1.fastq"),
          path("${sra_id}_2.fastq")

    script:
    """
    prefetch ${sra_id}
    fasterq-dump \
        --split-files \
        --threads ${task.cpus} \
        ${sra_id}
    rm -rf ${sra_id}
    """

    stub:
    """
    touch ${sra_id}_1.fastq
    touch ${sra_id}_2.fastq
    """
}
process FASTP_QC {
    tag "$sra_id"
    conda "${System.getenv('HOME')}/miniconda/envs/ngs_tools"
    cpus 4
    publishDir "${params.outdir}/clean_fastq",
        mode: 'copy',
        pattern: "*_clean.fastq"
    publishDir "${params.outdir}/qc",
        mode: 'copy',
        pattern: "*.{html,json}"

    input:
    tuple val(sra_id),
          path(r1),
          path(r2)

    output:
    tuple val(sra_id),
          path("${sra_id}_1_clean.fastq"),
          path("${sra_id}_2_clean.fastq"),
          emit: reads
    path "${sra_id}_fastp.html"
    path "${sra_id}_fastp.json"

    script:
    """
    fastp \
        -i ${r1} \
        -I ${r2} \
        -o ${sra_id}_1_clean.fastq \
        -O ${sra_id}_2_clean.fastq \
        -h ${sra_id}_fastp.html \
        -j ${sra_id}_fastp.json \
        --thread ${task.cpus}
    """

    stub:
    """
    touch ${sra_id}_1_clean.fastq
    touch ${sra_id}_2_clean.fastq
    touch ${sra_id}_fastp.html
    touch ${sra_id}_fastp.json
    """
}
process ALIGN_AND_SORT {
    tag "$sra_id"
    conda "${System.getenv('HOME')}/miniconda/envs/ngs_tools"
    cpus 4
    publishDir "${params.outdir}/bam",
        mode: 'copy',
        pattern: "*.bam*"

    input:
    tuple val(sra_id),
          path(r1),
          path(r2)
    path ref

    output:
    tuple val(sra_id),
          path("${sra_id}.sorted.bam"),
          path("${sra_id}.sorted.bam.bai")

    script:
    """
    bwa mem \
        -t ${task.cpus} \
        ${ref} \
        ${r1} \
        ${r2} \
    | samtools sort \
        -@ ${task.cpus} \
        -o ${sra_id}.sorted.bam -
    samtools index ${sra_id}.sorted.bam
    """

    stub:
    """
    touch ${sra_id}.sorted.bam
    touch ${sra_id}.sorted.bam.bai
    """
}
process VARIANT_CALLING {
    tag "$sra_id"
    conda "${System.getenv('HOME')}/miniconda/envs/ngs_tools"
    cpus 2
    publishDir "${params.outdir}/vcf",
        mode: 'copy'

    input:
    tuple val(sra_id),
          path(bam),
          path(bai)
    path ref

    output:
    path "${sra_id}.vcf.gz"
    path "${sra_id}.vcf.gz.csi"

    script:
    """
    bcftools mpileup \
        -f ${ref} \
        -Ou \
        ${bam} \
    | bcftools call \
        -mv \
        -Oz \
        -o ${sra_id}.vcf.gz
    bcftools index ${sra_id}.vcf.gz
    """

    stub:
    """
    touch ${sra_id}.vcf.gz
    touch ${sra_id}.vcf.gz.csi
    """
}

workflow {
    sra_ch = Channel
        .fromPath(params.samples)
        .splitText()
        .map { it.trim() }
        .filter { it }

    ref = file(params.ref)

    FETCH_FASTQ(sra_ch)
    FASTP_QC(FETCH_FASTQ.out)
    ALIGN_AND_SORT(
        FASTP_QC.out.reads,
        ref
    )
    VARIANT_CALLING(
        ALIGN_AND_SORT.out,
        ref
    )
}
