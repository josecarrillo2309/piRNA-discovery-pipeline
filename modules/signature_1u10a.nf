process SIGNATURE_1U10A {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/07_signatures", mode: 'copy'

    conda "conda-forge::python=3.10"
    container "quay.io/biocontainers/python:3.10"

    input:
    tuple val(sample_id), path(bam), path(fasta)

    output:
    path "${sample_id}_bias_report.txt"

    script:
    """
    python3 ${params.scripts}/pirna_signatures.py -i $fasta -o ${sample_id}_bias_report.txt
    """
}
