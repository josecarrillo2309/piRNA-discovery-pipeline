#!/bin/bash
set -e

echo "Generating synthetic test data for CI/CD..."

# Create directories
mkdir -p .github_test_data/reads
mkdir -p .github_test_data/reference

# 1. Generate a tiny synthetic genome (10kb)
echo ">chr1" > .github_test_data/reference/dummy_genome.fa
# Generate random DNA sequence (1000 lines of 100 bases)
tr -dc A-T-C-G < /dev/urandom | head -c 10000 | fold -w 100 >> .github_test_data/reference/dummy_genome.fa

# 2. Generate a tiny dummy Rfam/miRBase (negative filter)
echo ">mmu-mir-1" > .github_test_data/reference/dummy_rfam.fa
echo "TGAGGTAGTAGGTTGTATAGTT" >> .github_test_data/reference/dummy_rfam.fa

# 3. Generate a tiny dummy piRBase
echo ">piR-mmu-1" > .github_test_data/reference/dummy_pirbase.fa
echo "TGGGTGGGGGGCGTCCTAGGCGATCGA" >> .github_test_data/reference/dummy_pirbase.fa
gzip -f .github_test_data/reference/dummy_pirbase.fa

# 4. Generate dummy synthetic small RNA reads (FASTQ)
# Some reads matching genome, some matching miRNAs, some with adapter (TGGAATTCTCGGGTGCCAAGG)
cat << 'EOF' > .github_test_data/reads/test_sample.fastq
@seq1_adapter
TGGGTGGGGGGCGTCCTAGGCGATCGATGGAATTCTCGGGTGCCAAGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq2_mirna
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGTGCCAAGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq3_genomic
ATCGATCGATCGATCGATCGATCGTGGAATTCTCGGGTGCCAAGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# 5. Build Bowtie indexes
echo "Building Bowtie indexes (this will be instant)..."
bowtie-build .github_test_data/reference/dummy_genome.fa .github_test_data/reference/dummy_genome
bowtie-build .github_test_data/reference/dummy_rfam.fa .github_test_data/reference/dummy_rfam

echo "Dummy data generation complete."
