process FILTER_NCRS {
    tag "$sample_id"
    publishDir "${params.outdir}/04_filtered_non_pirna", mode: 'copy'

    conda "bioconda::bowtie=1.3.1"
    container "quay.io/biocontainers/bowtie:1.3.1--py310h4b830d6_3"

    input:
    tuple val(sample_id), path(reads)
    path index_dir

    output:
    tuple val(sample_id), path("${sample_id}_putative_pirnas.fastq"), emit: putative_reads
    path "${sample_id}_bowtie_filter.log", emit: log

    script:
    """
    bowtie -S -p ${task.cpus} -v 1 --best --strata --un ${sample_id}_putative_pirnas.fastq ${index_dir}/${params.rfam_mirbase_prefix} $reads > discarded.sam 2> ${sample_id}_bowtie_filter.log
    rm discarded.sam
    """
}
