#!/usr/bin/env python3

import argparse
import sys
import gzip

def read_fasta_seqs(fasta_file):
    """Reads a FASTA file (can be gzipped) and returns a set of uppercase sequences."""
    seqs = set()
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    mode = 'rt' if fasta_file.endswith('.gz') else 'r'
    
    try:
        with open_func(fasta_file, mode) as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>') and line:
                    seqs.add(line.upper().replace('U', 'T'))
    except FileNotFoundError:
        print(f"Error: Could not find FASTA file {fasta_file}")
        sys.exit(1)
    return seqs

def compare_to_pirbase(candidates_fasta, pirbase_fasta, output_report):
    print("Loading putative candidates...")
    candidate_seqs = read_fasta_seqs(candidates_fasta)
    print(f"Loaded {len(candidate_seqs)} unique putative piRNA sequences from our pipeline.")
    
    if len(candidate_seqs) == 0:
        msg = "No candidates to compare."
        print(msg)
        with open(output_report, 'w') as f:
            f.write(msg + "\n")
        sys.exit(0)
        
    print("Streaming piRBase known sequences to find matches (this saves memory)...")
    known_candidates = set()
    
    open_func = gzip.open if pirbase_fasta.endswith('.gz') else open
    mode = 'rt' if pirbase_fasta.endswith('.gz') else 'r'
    
    try:
        with open_func(pirbase_fasta, mode) as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>') and line:
                    seq = line.upper().replace('U', 'T')
                    if seq in candidate_seqs:
                        known_candidates.add(seq)
                        
                        # Optimization: if all candidates are found, we can stop early
                        if len(known_candidates) == len(candidate_seqs):
                            break
    except FileNotFoundError:
        print(f"Error: Could not find piRBase FASTA file {pirbase_fasta}")
        sys.exit(1)
            
    known_hits = len(known_candidates)
    novel_hits = len(candidate_seqs) - known_hits
    
    pct_known = (known_hits / len(candidate_seqs)) * 100
    pct_novel = (novel_hits / len(candidate_seqs)) * 100
    
    report = []
    report.append("==================================================")
    report.append("          piRNA Known vs Novel Report             ")
    report.append("==================================================")
    report.append(f"Total Unique Candidates Analyzed: {len(candidate_seqs)}")
    report.append("--------------------------------------------------")
    report.append(f"Previously Reported (in piRBase):  {known_hits} ({pct_known:.2f}%)")
    report.append(f"Novel / Unreported Candidates:     {novel_hits} ({pct_novel:.2f}%)")
    report.append("==================================================")
    
    report_str = "\n".join(report)
    print("\n" + report_str)
    
    with open(output_report, 'w') as f:
        f.write(report_str)
    print(f"\nReport saved to: {output_report}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare identified piRNAs against piRBase.")
    parser.add_argument("-c", "--candidates", required=True, help="FASTA of putative piRNAs")
    parser.add_argument("-p", "--pirbase", required=True, help="FASTA of known piRBase sequences (can be .gz)")
    parser.add_argument("-o", "--output", required=True, help="Output report text file")
    
    args = parser.parse_args()
    compare_to_pirbase(args.candidates, args.pirbase, args.output)
