process EXTRACT_FASTA {
    tag "$sample_id"
    publishDir "${params.outdir}/06_final_candidates", mode: 'copy'

    conda "bioconda::samtools=1.18"
    container "quay.io/biocontainers/samtools:1.18--h50ea8bc_1"

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}_mapped_pirnas.bam"), path("${sample_id}_final_candidates.fasta")

    script:
    """
    samtools view -@ ${task.cpus} -bS $sam | samtools sort -m 500M -@ ${task.cpus} -o ${sample_id}_mapped_pirnas.bam
    samtools index ${sample_id}_mapped_pirnas.bam
    samtools fasta ${sample_id}_mapped_pirnas.bam > ${sample_id}_final_candidates.fasta
    """
}
