process TRIM_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/02_trimmed_filtered", mode: 'copy'
    
    conda "bioconda::cutadapt=4.9"
    container "quay.io/biocontainers/cutadapt:4.9--py310h0941dc4_0"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq"), emit: trimmed
    path "${sample_id}_cutadapt.log", emit: log

    script:
    """
    cutadapt -j ${task.cpus} -a ${params.adapter} -m 24 -M 32 -o ${sample_id}_trimmed.fastq $reads > ${sample_id}_cutadapt.log
    """
}
