process PIRNA_CLUSTERS {
    tag "$sample_id"
    publishDir "${params.outdir}/08_clusters", mode: 'copy'

    conda "conda-forge::python=3.10"
    container "quay.io/biocontainers/python:3.10"

    input:
    tuple val(sample_id), path(bed)

    output:
    path "${sample_id}_clusters.bed"
    path "${sample_id}_cluster_report.txt"

    script:
    """
    python3 ${params.scripts}/pirna_clusters.py -i $bed --max_dist 5000 --min_reads 50 --min_len 1000 --out_bed ${sample_id}_clusters.bed --out_report ${sample_id}_cluster_report.txt
    """
}
