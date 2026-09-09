process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/03_fastqc_clean", mode: 'copy'

    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc -t ${task.cpus} -q $reads
    """
}
