# piRNA Discovery Pipeline

A robust, memory-efficient, and fully automated bioinformatics pipeline for the discovery and characterization of PIWI-interacting RNAs (piRNAs) from small RNA sequencing data.

## Overview
This pipeline is designed to accurately filter, map, and identify piRNAs while minimizing memory usage (Strict 6GB RAM limit built-in). It performs the following steps:
1. **Quality Control & Trimming**: Raw read assessment (`FastQC`) and adapter removal (`Cutadapt`).
2. **Negative Filtering**: Removal of structural and non-coding RNAs (rRNA, tRNA, miRNA) against Rfam and miRBase using `Bowtie 1`.
3. **Genome Mapping**: Exact/near-exact mapping of surviving reads to the reference genome (e.g., mm39).
4. **piRBase Comparison**: Identification of known vs. novel piRNA candidates against the piRBase database.
5. **Signature Analysis**: Statistical profiling of canonical piRNA features (1U bias and 10A bias).
6. **Ping-Pong Amplification Signature**: O(1) memory streaming analysis to detect the 10nt overlap signature of the Ping-Pong cycle.
7. **Genomic Clustering**: Identification of genomic loci producing piRNA clusters.

## Prerequisites
- Linux / Unix environment
- **Conda / Mamba** (for automatic dependency resolution)
- **Nextflow** (>= 22.10.x)

## Data Preparation (Inputs)
To run this pipeline on your own biological samples, you must prepare the following input files:

### 1. Biological Data (Mandatory)
- **Raw reads:** Small RNA-seq data in FASTQ format (`.fastq` or `.fastq.gz`).

### 2. Reference Genome Index (Mandatory)
- Download the FASTA file for your organism's reference genome (e.g., `mm39` for mouse).
- Build the Bowtie 1 index: `bowtie-build genome.fa genome_index`

### 3. Negative Filter Index (Mandatory)
- Download structural RNAs (rRNA, tRNA, snRNA) from **Rfam** and microRNAs from **miRBase**.
- Concatenate them into a single FASTA file and build the Bowtie 1 index: `bowtie-build trash_rna.fa negative_filter_index`

### 4. piRBase Database (Optional)
- Download the known piRNAs FASTA file for your organism from [piRBase](http://bigdata.ibp.ac.cn/piRBase/). If provided, the pipeline will classify discoveries as *Known* vs *Novel*.

## Execution

The pipeline is orchestrated by Nextflow, which automatically handles parallelization, environment isolation (via Conda), and file staging.

To execute the pipeline, point it to your prepared inputs using the terminal:

```bash
nextflow run main.nf \
  -profile conda \
  --reads "path/to/data/*.fastq.gz" \
  --bowtie_index_dir "path/to/indexes/" \
  --genome_prefix "genome_index" \
  --rfam_mirbase_prefix "negative_filter_index" \
  --pirbase "path/to/pirbase.fa" \
  --adapter "TGGAATTCTCGGGTGCCAAGG" \
  --threads 8
```

## Output Structure
The `results/` directory will contain organized outputs per sample:
- `01_fastqc_raw/` & `03_fastqc_trimmed/`: Individual HTML quality reports.
- `06_final_candidates/`: Filtered and mapped sequences in FASTA format (Putative piRNAs).
- `07_signatures/`: Text reports detailing 1U/10A biases and Ping-Pong distances.
- `08_clusters/`: BED files of identified piRNA clusters and statistical reports.
- `09_pirbase_comparison/`: Text reports comparing your candidates against known piRNAs.

## Automated CI/CD
This repository includes a GitHub Actions workflow (`.github/workflows/ci.yml`) that automatically tests the integrity of the pipeline using synthetically generated data on every push to the `main` branch.

## License
MIT License.
