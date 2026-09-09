#!/usr/bin/env python3

import argparse
import sys

def find_clusters(bed_file, max_dist=5000, min_reads=50, min_len=1000, out_bed=None, out_report=None):
    """
    Identifies piRNA clusters using a sliding window / distance-based grouping approach.
    Assumes the input BED file is sorted by chromosome and start coordinate.
    """
    print(f"Reading and clustering piRNAs from: {bed_file}")
    
    clusters = []
    
    current_chrom = None
    current_start = -1
    current_end = -1
    current_reads = 0
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                
                # Initialize first cluster
                if current_chrom is None:
                    current_chrom = chrom
                    current_start = start
                    current_end = end
                    current_reads = 1
                    continue
                
                # If same chromosome and within max_dist of the current cluster's end
                if chrom == current_chrom and start <= (current_end + max_dist):
                    # Extend the cluster
                    current_end = max(current_end, end)
                    current_reads += 1
                else:
                    # Save the finalized cluster
                    clusters.append({
                        'chrom': current_chrom,
                        'start': current_start,
                        'end': current_end,
                        'reads': current_reads,
                        'length': current_end - current_start
                    })
                    
                    # Start a new cluster
                    current_chrom = chrom
                    current_start = start
                    current_end = end
                    current_reads = 1
                    
            # Don't forget the last cluster in the file
            if current_chrom is not None:
                clusters.append({
                    'chrom': current_chrom,
                    'start': current_start,
                    'end': current_end,
                    'reads': current_reads,
                    'length': current_end - current_start
                })
                
    except FileNotFoundError:
        print(f"Error: BED file not found at {bed_file}")
        sys.exit(1)

    # Filter clusters based on user-defined thresholds
    valid_clusters = [c for c in clusters if c['length'] >= min_len and c['reads'] >= min_reads]
    
    # Sort clusters by number of reads (density proxy) descending
    valid_clusters.sort(key=lambda x: x['reads'], reverse=True)

    # 1. Output BED file for genome browsers (IGV, UCSC)
    if out_bed:
        with open(out_bed, 'w') as bed_out:
            for i, c in enumerate(valid_clusters):
                cluster_name = f"piRNA_cluster_{i+1}_reads:{c['reads']}"
                # BED6 format
                bed_out.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{cluster_name}\t{c['reads']}\t.\n")
        print(f"Cluster coordinates saved to: {out_bed}")

    # 2. Output Text Report
    report = []
    report.append("==================================================")
    report.append("             piRNA Cluster Analysis               ")
    report.append("==================================================")
    report.append(f"Parameters:")
    report.append(f" - Max distance between reads: {max_dist} nt")
    report.append(f" - Minimum cluster length: {min_len} nt")
    report.append(f" - Minimum reads per cluster: {min_reads}")
    report.append("--------------------------------------------------")
    report.append(f"Total clusters found: {len(valid_clusters)}")
    report.append("--------------------------------------------------")
    
    if len(valid_clusters) > 0:
        report.append(f"{'Rank':<5} | {'Chromosome':<10} | {'Start':<12} | {'End':<12} | {'Length(nt)':<12} | {'Reads':<8}")
        report.append("-" * 75)
        for i, c in enumerate(valid_clusters):
            # Only print top 20 to console to avoid spam
            line = f"{i+1:<5} | {c['chrom']:<10} | {c['start']:<12} | {c['end']:<12} | {c['length']:<12} | {c['reads']:<8}"
            if i < 20:
                report.append(line)
            
            # Save all to file if output_report is provided
            if out_report and i >= 20:
                pass # Handled below
        
        if len(valid_clusters) > 20:
            report.append(f"... and {len(valid_clusters) - 20} more clusters (see full report).")
    
    report_str = "\n".join(report)
    print(report_str)
    
    if out_report:
        with open(out_report, 'w') as rep_out:
            # Write full report including all clusters
            rep_out.write("==================================================\n")
            rep_out.write("             piRNA Cluster Analysis               \n")
            rep_out.write("==================================================\n")
            rep_out.write(f"Total clusters found: {len(valid_clusters)}\n")
            rep_out.write("-" * 75 + "\n")
            rep_out.write(f"{'Rank':<5} | {'Chromosome':<10} | {'Start':<12} | {'End':<12} | {'Length(nt)':<12} | {'Reads':<8}\n")
            rep_out.write("-" * 75 + "\n")
            for i, c in enumerate(valid_clusters):
                rep_out.write(f"{i+1:<5} | {c['chrom']:<10} | {c['start']:<12} | {c['end']:<12} | {c['length']:<12} | {c['reads']:<8}\n")
        print(f"Full report saved to: {out_report}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify piRNA clusters from mapped BED file.")
    parser.add_argument("-i", "--input", required=True, help="Input BED file (sorted by coordinate)")
    parser.add_argument("--max_dist", type=int, default=5000, help="Max gap between reads to be in same cluster (default: 5000)")
    parser.add_argument("--min_reads", type=int, default=50, help="Minimum reads required to form a cluster (default: 50)")
    parser.add_argument("--min_len", type=int, default=1000, help="Minimum cluster length in nt (default: 1000)")
    parser.add_argument("--out_bed", help="Output BED file for genomic visualization")
    parser.add_argument("--out_report", help="Output text file with the cluster summary")
    
    args = parser.parse_args()
    find_clusters(args.input, args.max_dist, args.min_reads, args.min_len, args.out_bed, args.out_report)
