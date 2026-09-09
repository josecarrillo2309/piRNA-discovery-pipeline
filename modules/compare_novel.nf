process COMPARE_NOVEL {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/09_pirbase_comparison", mode: 'copy'

    conda "conda-forge::python=3.10"
    container "quay.io/biocontainers/python:3.10"

    input:
    tuple val(sample_id), path(bam), path(fasta)
    path pirbase_fa

    output:
    path "${sample_id}_novelty_report.txt"

    script:
    """
    python3 ${params.scripts}/compare_pirbase.py -c $fasta -p $pirbase_fa -o ${sample_id}_novelty_report.txt
    """
}
