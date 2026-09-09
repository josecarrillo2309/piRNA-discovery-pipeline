process PINGPONG_SIG {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/07_signatures", mode: 'copy'

    conda "conda-forge::python=3.10"
    container "quay.io/biocontainers/python:3.10"

    input:
    tuple val(sample_id), path(bed)

    output:
    path "${sample_id}_pingpong_report.txt"

    script:
    """
    python3 ${params.scripts}/ping_pong_signature.py -i $bed -o ${sample_id}_pingpong_report.txt
    """
}
