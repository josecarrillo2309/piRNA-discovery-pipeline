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
- **SRA-Toolkit** (`fasterq-dump`) (Optional, for downloading test datasets)

## Setup and Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/project_piwirnas.git
   cd project_piwirnas
   ```

2. **Download databases and build indexes:**
   Run the automated setup script. This will download the mouse genome (mm39), miRBase, Rfam, piRBase, and build the required Bowtie memory-optimized indexes.
   ```bash
   nohup ./download_and_setup.sh > setup.log 2>&1 &
   ```
   *Note: This step requires a stable internet connection and will download several gigabytes of data.*

## Execution
The pipeline is orchestrated by Nextflow, which automatically handles parallelization, environment isolation (via Conda), and file staging.

To run the pipeline on the test datasets downloaded during setup:
```bash
nextflow run main.nf
```

To run the pipeline on your own FASTQ data, modify the `reads` parameter in `nextflow.config`:
```groovy
params {
    reads = "/path/to/your/data/*.fastq"
    // ...
}
```

### Advanced: Standalone Bash Pipeline
If you do not wish to use Nextflow, a standalone Bash script is provided that executes the pipeline sequentially on a single FASTQ file:
```bash
./pirna_pipeline.sh <input_fastq> <output_dir> <threads>
```

## Output Structure
The `results/` directory will contain organized outputs per sample:
- `01_fastqc_raw/` & `03_fastqc_trimmed/`: HTML quality reports.
- `06_final_candidates/`: Filtered and mapped sequences in FASTA format.
- `07_signatures/`: Text reports detailing 1U/10A biases and Ping-Pong distances.
- `08_clusters/`: BED files of identified piRNA clusters and statistical reports.
- `09_pirbase_comparison/`: Text reports comparing your candidates against known piRNAs.

## License
MIT License.
