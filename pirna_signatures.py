#!/usr/bin/env python3

import argparse
from collections import Counter
import sys

def analyze_signatures(fasta_file, output_file=None):
    """
    Analyzes a FASTA file of small RNAs for piRNA signatures:
    - 1U bias (Uracil/Thymine at position 1)
    - 10A bias (Adenine at position 10)
    """
    pos1_counts = Counter()
    pos10_counts = Counter()
    total_reads = 0

    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    seq = line.upper().replace('U', 'T') # Normalize to DNA alphabet if needed
                    if len(seq) >= 10:
                        pos1_counts[seq[0]] += 1
                        pos10_counts[seq[9]] += 1
                        total_reads += 1
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {fasta_file}")
        sys.exit(1)

    if total_reads == 0:
        msg = "No valid sequences found (>= 10nt) to analyze."
        print(msg)
        if output_file:
            with open(output_file, 'w') as out:
                out.write(msg + "\n")
        sys.exit(0)

    # Calculate percentages
    pos1_u_pct = (pos1_counts.get('T', 0) / total_reads) * 100
    pos10_a_pct = (pos10_counts.get('A', 0) / total_reads) * 100

    report = []
    report.append(f"--- piRNA Signature Analysis ---")
    report.append(f"Total sequences analyzed: {total_reads}")
    report.append(f"Position 1 'U/T' Bias: {pos1_u_pct:.2f}% (Expected > ~70% for primary piRNAs)")
    report.append(f"Position 10 'A' Bias:  {pos10_a_pct:.2f}% (Expected > ~70% for secondary ping-pong piRNAs)")
    report.append(f"--------------------------------")
    
    report.append("\nPosition 1 Nucleotide Distribution:")
    for nt in ['A', 'C', 'G', 'T']:
        pct = (pos1_counts.get(nt, 0) / total_reads) * 100
        report.append(f"  {nt}: {pct:.2f}%")

    report.append("\nPosition 10 Nucleotide Distribution:")
    for nt in ['A', 'C', 'G', 'T']:
        pct = (pos10_counts.get(nt, 0) / total_reads) * 100
        report.append(f"  {nt}: {pct:.2f}%")

    report_str = "\n".join(report)
    
    print(report_str)
    
    if output_file:
        with open(output_file, 'w') as out:
            out.write(report_str)
        print(f"\nReport saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for 1U and 10A piRNA sequence biases.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file of putative piRNAs")
    parser.add_argument("-o", "--output", help="Output text file for the report")
    
    args = parser.parse_args()
    analyze_signatures(args.input, args.output)
