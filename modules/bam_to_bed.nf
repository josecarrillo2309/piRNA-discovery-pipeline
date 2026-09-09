process BAM_TO_BED {
    tag "$sample_id"
    
    conda "bioconda::bedtools=2.31.1"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"

    input:
    tuple val(sample_id), path(bam), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_mapped.bed")

    script:
    """
    bedtools bamtobed -i $bam > ${sample_id}_mapped.bed
    """
}
