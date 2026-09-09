process MULTIQC {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    conda "bioconda::multiqc=1.21.0"
    container "quay.io/biocontainers/multiqc:1.21.0--pyhdfd78af_0"

    input:
    path multiqc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/", emit: data

    script:
    """
    multiqc .
    """
}
