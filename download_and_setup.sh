#!/bin/bash

# ==============================================================================
# piRNA Pipeline Data Setup Script
# ==============================================================================
# Description: Downloads reference genomes, non-piRNA databases, and builds Bowtie indexes.
# It also attempts to download 3 test SRA datasets for pipeline validation.
# ==============================================================================

set -eo pipefail

# --- Configuration ---
DATA_DIR="./reference_data"
TEST_DATA_DIR="./test_datasets"
BOWTIE_INDEX_DIR="${DATA_DIR}/bowtie_indexes"
THREADS=4

mkdir -p "$DATA_DIR" "$TEST_DATA_DIR" "$BOWTIE_INDEX_DIR"

echo "============================================="
echo "1. Downloading Databases for Negative Filtering"
echo "============================================="

# 1. miRBase (mature miRNAs)
echo "Downloading miRBase..."
wget -c -O "${DATA_DIR}/mature.fa" "https://www.mirbase.org/download/mature.fa"
# Extract only mouse (mmu) miRNAs
grep -A 1 "^>mmu" "${DATA_DIR}/mature.fa" | grep -v "^--" > "${DATA_DIR}/mmu_mirbase.fa"

# 2. Rfam (rRNA, tRNA, etc.)
echo "Downloading Rfam (this might take a while)..."
wget -c -O "${DATA_DIR}/Rfam.fasta.gz" "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz"
zcat "${DATA_DIR}/Rfam.fasta.gz" > "${DATA_DIR}/Rfam.fasta"

# Combine Rfam and miRBase for the negative filter index
cat "${DATA_DIR}/mmu_mirbase.fa" "${DATA_DIR}/Rfam.fasta" > "${DATA_DIR}/mouse_rfam_mirbase.fa"

if [ ! -f "${BOWTIE_INDEX_DIR}/mouse_rfam_mirbase.1.ebwt" ]; then
    echo "Building Bowtie Index for Rfam/miRBase (with strict memory limits)..."
    bowtie-build --threads "$THREADS" --bmaxdivn 16 --dcv 256 "${DATA_DIR}/mouse_rfam_mirbase.fa" "${BOWTIE_INDEX_DIR}/mouse_rfam_mirbase"
else
    echo "Rfam/miRBase Bowtie Index already exists, skipping build."
fi


echo "============================================="
echo "2. Downloading Reference Genome (mm39)"
echo "============================================="
echo "Downloading Mouse genome mm39 from UCSC (Approx. 900MB compressed)..."
wget -c -O "${DATA_DIR}/mm39.fa.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"

if [ ! -f "${BOWTIE_INDEX_DIR}/mm39.1.ebwt" ]; then
    echo "Decompressing mm39 genome for index building..."
    zcat "${DATA_DIR}/mm39.fa.gz" > "${DATA_DIR}/mm39.fa"
    
    echo "Building Bowtie Index for mm39 (with strict memory limits)..."
    bowtie-build --threads "$THREADS" --bmaxdivn 16 --dcv 256 "${DATA_DIR}/mm39.fa" "${BOWTIE_INDEX_DIR}/mm39"
    
    echo "Removing temporary uncompressed genome to save space..."
    rm "${DATA_DIR}/mm39.fa"
else
    echo "mm39 Bowtie Index already exists, skipping build."
fi

echo "============================================="
echo "3. Downloading piRBase (Validation Database)"
echo "============================================="
echo "Downloading piRBase v3.0 for mouse..."
wget -c -O "${DATA_DIR}/piR_mouse.fa.gz" "http://www.regulatoryrna.org/database/piRBase/download/v3.0/fasta/piR_mouse_v3.0.fa.gz"


echo "============================================="
echo "4. Downloading Test Datasets (SRA)"
echo "============================================="
# Using 3 accessions from small RNA-seq experiments in Mouse Testis (high piRNA abundance)
# Note: Ensure you have SRA-Toolkit (fasterq-dump) installed.
SRA_ACCESSIONS=("SRR1197125" "SRR1197126" "SRR1197127") # Example accessions (Mouse testis small RNA)

if command -v fasterq-dump &> /dev/null; then
    for sra in "${SRA_ACCESSIONS[@]}"; do
        echo "Downloading dataset: $sra"
        fasterq-dump -e "$THREADS" -m 1G -O "$TEST_DATA_DIR" "$sra"
    done
else
    echo "WARNING: 'fasterq-dump' (sra-toolkit) not found. Skipping SRA test dataset downloads."
    echo "Please install sra-toolkit (e.g., sudo apt install sra-toolkit) and run this section manually,"
    echo "or download your datasets using an alternative method."
fi

echo "============================================="
echo "Setup Complete!"
echo "Indexes are located in: $BOWTIE_INDEX_DIR"
echo "Test datasets are in: $TEST_DATA_DIR"
echo "============================================="
