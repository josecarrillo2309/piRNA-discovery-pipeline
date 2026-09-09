#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Modules Import ---
include { FASTQC_RAW }      from './modules/fastqc_raw.nf'
include { TRIM_READS }      from './modules/trim_reads.nf'
include { FASTQC_TRIMMED }  from './modules/fastqc_trimmed.nf'
include { FILTER_NCRS }     from './modules/filter_ncrs.nf'
include { MAP_GENOME }      from './modules/map_genome.nf'
include { EXTRACT_FASTA }   from './modules/extract_fasta.nf'
include { SIGNATURE_1U10A } from './modules/signature_1u10a.nf'
include { BAM_TO_BED }      from './modules/bam_to_bed.nf'
include { PINGPONG_SIG }    from './modules/pingpong_sig.nf'
include { PIRNA_CLUSTERS }  from './modules/pirna_clusters.nf'
include { COMPARE_NOVEL }   from './modules/compare_novel.nf'
include { MULTIQC }         from './modules/multiqc.nf'


// --- Main Workflow ---

workflow {
    reads_ch = Channel.fromPath(params.reads)
        .map { file -> tuple(file.baseName, file) }

    // QC & Trimming
    FASTQC_RAW(reads_ch)
    trimmed_ch = TRIM_READS(reads_ch)
    FASTQC_TRIMMED(trimmed_ch.trimmed)

    // Filtering & Mapping
    putative_ch = FILTER_NCRS(trimmed_ch.trimmed, params.bowtie_index_dir)
    mapped_sam_ch = MAP_GENOME(putative_ch.putative_reads, params.bowtie_index_dir)
    
    // Formatting
    bam_fasta_ch = EXTRACT_FASTA(mapped_sam_ch.sam)

    // Signatures & Novelty
    SIGNATURE_1U10A(bam_fasta_ch)
    
    if (file(params.pirbase).exists()) {
        COMPARE_NOVEL(bam_fasta_ch, params.pirbase)
    } else {
        log.warn "piRBase file not found at ${params.pirbase}. Skipping known vs novel comparison."
    }

    // Ping-Pong & Clusters require BED format
    bed_ch = BAM_TO_BED(bam_fasta_ch)
    PINGPONG_SIG(bed_ch)
    PIRNA_CLUSTERS(bed_ch)

    // MultiQC Report
    ch_multiqc = Channel.empty()
    ch_multiqc = ch_multiqc.mix(FASTQC_RAW.out.reports)
                           .mix(TRIM_READS.out.log)
                           .mix(FASTQC_TRIMMED.out.reports)
                           .mix(FILTER_NCRS.out.log)
                           .mix(MAP_GENOME.out.log)
    MULTIQC(ch_multiqc.collect())
}
