process MAP_GENOME {
    tag "$sample_id"
    publishDir "${params.outdir}/05_genome_mapping", mode: 'copy'

    conda "bioconda::bowtie=1.3.1"
    container "quay.io/biocontainers/bowtie:1.3.1--py310h4b830d6_3"

    input:
    tuple val(sample_id), path(reads)
    path index_dir

    output:
    tuple val(sample_id), path("${sample_id}_mapped_pirnas.sam"), emit: sam
    path "${sample_id}_bowtie_map.log", emit: log

    script:
    """
    bowtie -S -p ${task.cpus} -v 1 -m 50 -a --best --strata ${index_dir}/${params.genome_prefix} $reads > ${sample_id}_mapped_pirnas.sam 2> ${sample_id}_bowtie_map.log || { cat ${sample_id}_bowtie_map.log; exit 1; }
    """
}
