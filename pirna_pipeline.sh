#!/bin/bash

# ==============================================================================
# piRNA Identification Pipeline (Mouse) - Version 1.0
# ==============================================================================
# Description: Pipeline to identify piRNAs from small RNA-seq data.
# Usage: ./pirna_pipeline.sh <input_fastq> <output_dir> <threads>
# ==============================================================================

set -eo pipefail # Exit immediately if a command exits with a non-zero status

# --- Input Arguments ---
INPUT_FASTQ=$1
OUTPUT_DIR=$2
THREADS=${3:-4} # Default to 4 threads if not provided

if [ -z "$INPUT_FASTQ" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_fastq> <output_dir> [threads]"
    exit 1
fi

# --- Configuration & Paths ---
# TODO: Update these paths to where the reference indexes are actually stored
BOWTIE_INDEX_DIR="/path/to/bowtie/indexes"
RFAM_MIRBASE_INDEX="${BOWTIE_INDEX_DIR}/mouse_rfam_mirbase"
MOUSE_GENOME_INDEX="${BOWTIE_INDEX_DIR}/mm39"

# 3' Adapter sequence (Example: Illumina TruSeq small RNA adapter. Change if necessary)
ADAPTER_SEQ="TGGAATTCTCGGGTGCCAAGG" 

# Log file
LOG_FILE="${OUTPUT_DIR}/pipeline.log"
mkdir -p "$OUTPUT_DIR"

echo "[$(date)] Starting piRNA Pipeline on $INPUT_FASTQ" | tee -a "$LOG_FILE"
echo "[$(date)] Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"

# ==============================================================================
# STEP 1: Quality Control (Raw)
# ==============================================================================
echo "[$(date)] STEP 1: Running FastQC on raw reads..." | tee -a "$LOG_FILE"
mkdir -p "${OUTPUT_DIR}/01_fastqc_raw"
fastqc -t "$THREADS" -o "${OUTPUT_DIR}/01_fastqc_raw" "$INPUT_FASTQ" >> "$LOG_FILE" 2>&1

# ==============================================================================
# STEP 2 & 3: Adapter Trimming and Length Filtering (24-32nt)
# ==============================================================================
echo "[$(date)] STEP 2 & 3: Trimming adapters and filtering length (24-32nt) using Cutadapt..." | tee -a "$LOG_FILE"
TRIMMED_FASTQ="${OUTPUT_DIR}/02_trimmed_filtered/trimmed_filtered.fastq"
mkdir -p "${OUTPUT_DIR}/02_trimmed_filtered"

# -a: 3' adapter
# -m 24: minimum length 24
# -M 32: maximum length 32
# --discard-untrimmed: optional, depends on library prep
cutadapt -j "$THREADS" -a "$ADAPTER_SEQ" -m 24 -M 32 \
    -o "$TRIMMED_FASTQ" "$INPUT_FASTQ" >> "$LOG_FILE" 2>&1

echo "[$(date)] Running FastQC on trimmed/filtered reads..." | tee -a "$LOG_FILE"
mkdir -p "${OUTPUT_DIR}/03_fastqc_clean"
fastqc -t "$THREADS" -o "${OUTPUT_DIR}/03_fastqc_clean" "$TRIMMED_FASTQ" >> "$LOG_FILE" 2>&1

# ==============================================================================
# STEP 4: Filtering Known Non-piRNA small RNAs (rRNA, tRNA, miRNA, etc.)
# ==============================================================================
echo "[$(date)] STEP 4: Filtering non-piRNAs (mapping to Rfam/miRBase)..." | tee -a "$LOG_FILE"
mkdir -p "${OUTPUT_DIR}/04_filtered_non_pirna"
UNMAPPED_FASTQ="${OUTPUT_DIR}/04_filtered_non_pirna/putative_pirnas.fastq"
MAPPED_DISCARD="${OUTPUT_DIR}/04_filtered_non_pirna/discarded_rnas.sam"

# Bowtie 1 parameters for small RNA:
# -v 1: allow max 1 mismatch
# --best --strata: report only the best alignments
# --un: save unmapped reads (these are our candidates)
bowtie -S -p "$THREADS" -v 1 --best --strata \
    --un "$UNMAPPED_FASTQ" \
    "$RFAM_MIRBASE_INDEX" "$TRIMMED_FASTQ" > "$MAPPED_DISCARD" 2>> "$LOG_FILE"

# Clean up SAM to save space
rm "$MAPPED_DISCARD"

# ==============================================================================
# STEP 5: Genome Mapping (Mouse mm39)
# ==============================================================================
echo "[$(date)] STEP 5: Mapping putative piRNAs to Mouse Genome (mm39)..." | tee -a "$LOG_FILE"
mkdir -p "${OUTPUT_DIR}/05_genome_mapping"
MAPPED_PIRNAS_SAM="${OUTPUT_DIR}/05_genome_mapping/mapped_pirnas.sam"
MAPPED_PIRNAS_BAM="${OUTPUT_DIR}/05_genome_mapping/mapped_pirnas.bam"

# -v 1: allow max 1 mismatch
# -m 50: suppress alignments if > 50 alignments exist (piRNAs can be repetitive, but we cap it)
# -a: report all valid alignments (up to the limit of -m)
bowtie -S -p "$THREADS" -v 1 -m 50 -a --best --strata \
    "$MOUSE_GENOME_INDEX" "$UNMAPPED_FASTQ" > "$MAPPED_PIRNAS_SAM" 2>> "$LOG_FILE"

# Convert SAM to BAM, sort, and index using samtools
echo "[$(date)] Converting SAM to sorted BAM..." | tee -a "$LOG_FILE"
samtools view -@ "$THREADS" -bS "$MAPPED_PIRNAS_SAM" | samtools sort -m 500M -@ "$THREADS" -o "$MAPPED_PIRNAS_BAM"
samtools index "$MAPPED_PIRNAS_BAM"
rm "$MAPPED_PIRNAS_SAM"

# ==============================================================================
# STEP 6: Extract Fasta for Downstream Signature Analysis
# ==============================================================================
echo "[$(date)] STEP 6: Extracting mapped sequences to FASTA..." | tee -a "$LOG_FILE"
FINAL_FASTA="${OUTPUT_DIR}/06_final_candidates/putative_pirnas_mapped.fasta"
mkdir -p "${OUTPUT_DIR}/06_final_candidates"

# Convert BAM back to FastQ then to Fasta to get the unique sequences for analysis
samtools fasta -@ "$THREADS" "$MAPPED_PIRNAS_BAM" > "$FINAL_FASTA"

# ==============================================================================
# STEP 7: Signature Analysis (1U/10A and Ping-Pong)
# ==============================================================================
echo "[$(date)] STEP 7: Running Signature Analysis..." | tee -a "$LOG_FILE"
SIGNATURES_DIR="${OUTPUT_DIR}/07_signatures"
mkdir -p "$SIGNATURES_DIR"

# 1U / 10A Bias Analysis
echo "[$(date)] Analyzing 1U/10A Sequence Bias..." | tee -a "$LOG_FILE"
BIAS_REPORT="${SIGNATURES_DIR}/sequence_bias_report.txt"
python3 "$(dirname "$0")/pirna_signatures.py" -i "$FINAL_FASTA" -o "$BIAS_REPORT" | tee -a "$LOG_FILE"

# Ping-Pong Signature Analysis
echo "[$(date)] Converting BAM to BED for Ping-Pong analysis..." | tee -a "$LOG_FILE"
MAPPED_PIRNAS_BED="${SIGNATURES_DIR}/mapped_pirnas.bed"
# Note: Requires bedtools installed
bedtools bamtobed -i "$MAPPED_PIRNAS_BAM" > "$MAPPED_PIRNAS_BED"

echo "[$(date)] Analyzing Ping-Pong Signature (10nt overlap)..." | tee -a "$LOG_FILE"
PINGPONG_REPORT="${SIGNATURES_DIR}/ping_pong_report.txt"
python3 "$(dirname "$0")/ping_pong_signature.py" -i "$MAPPED_PIRNAS_BED" -o "$PINGPONG_REPORT" | tee -a "$LOG_FILE"

# piRNA Cluster Identification
echo "[$(date)] Identifying piRNA Clusters..." | tee -a "$LOG_FILE"
CLUSTER_DIR="${OUTPUT_DIR}/08_clusters"
mkdir -p "$CLUSTER_DIR"
CLUSTER_BED="${CLUSTER_DIR}/pirna_clusters.bed"
CLUSTER_REPORT="${CLUSTER_DIR}/cluster_report.txt"

python3 "$(dirname "$0")/pirna_clusters.py" -i "$MAPPED_PIRNAS_BED" \
    --max_dist 5000 --min_reads 50 --min_len 1000 \
    --out_bed "$CLUSTER_BED" --out_report "$CLUSTER_REPORT" | tee -a "$LOG_FILE"



# piRNA Known vs Novel Comparison (piRBase)
echo "[$(date)] Comparing candidates against piRBase..." | tee -a "$LOG_FILE"
PIRBASE_FASTA="../reference_data/piR_mouse.fa.gz" # Assuming this relative path for now
COMPARISON_REPORT="${SIGNATURES_DIR}/pirbase_comparison_report.txt"

if [ -f "$PIRBASE_FASTA" ]; then
    python3 "$(dirname "$0")/compare_pirbase.py" -c "$FINAL_FASTA" -p "$PIRBASE_FASTA" -o "$COMPARISON_REPORT" | tee -a "$LOG_FILE"
else
    echo "WARNING: piRBase FASTA ($PIRBASE_FASTA) not found. Skipping known vs novel comparison." | tee -a "$LOG_FILE"
fi

echo "==============================================================================" | tee -a "$LOG_FILE"
echo "[$(date)] Pipeline completed successfully!" | tee -a "$LOG_FILE"
echo "[$(date)] Final candidates: $FINAL_FASTA" | tee -a "$LOG_FILE"
echo "[$(date)] Analysis reports: $SIGNATURES_DIR" | tee -a "$LOG_FILE"
echo "[$(date)] Clusters BED:     $CLUSTER_BED" | tee -a "$LOG_FILE"
echo "==============================================================================" | tee -a "$LOG_FILE"
