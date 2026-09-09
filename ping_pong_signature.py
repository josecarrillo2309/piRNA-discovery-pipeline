#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict, Counter

def calculate_ping_pong(bed_file, max_distance=30, output_file=None):
    """
    Calculates the distance between 5' ends of piRNAs on opposite strands.
    A sharp peak at 10nt indicates a strong ping-pong amplification signature.
    Input must be a BED6 file sorted by chromosome and start position.
    Memory efficient streaming implementation (O(1) RAM).
    """
    from collections import deque
    
    print(f"Streaming BED file for Ping-Pong analysis: {bed_file}")
    
    distances = {i: 0 for i in range(1, max_distance + 1)}
    buffer = deque()
    current_chrom = None
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6: continue
                
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[5]
                
                # If we switch chromosomes, clear the memory buffer
                if chrom != current_chrom:
                    buffer.clear()
                    current_chrom = chrom
                
                # Evict reads from buffer that are safely behind our sliding window
                # max distance is 30, max read length is ~32. A buffer of 100nt is mathematically safe
                while buffer and buffer[0]['start'] < start - 100:
                    buffer.popleft()
                
                # Calculate 5' end
                pos = start if strand == '+' else end
                
                # Compare current read against all recent reads in the sliding window buffer
                for b in buffer:
                    if strand != b['strand']: # Only compare opposite strands
                        if strand == '+':
                            p_pos = pos
                            m_pos = b['pos']
                        else:
                            m_pos = pos
                            p_pos = b['pos']
                            
                        dist = m_pos - p_pos
                        if 1 <= dist <= max_distance:
                            distances[dist] += 1
                            
                # Add current read to the sliding window
                buffer.append({'start': start, 'pos': pos, 'strand': strand})
                
    except FileNotFoundError:
        print(f"Error: BED file not found at {bed_file}")
        sys.exit(1)

    report = []
    report.append("--- Ping-Pong Signature Analysis (5' to 5' distance) ---")
    report.append("Distance (nt) | Pair Count")
    report.append("--------------------------")
    
    for i in range(1, max_distance + 1):
        count = distances.get(i, 0)
        report.append(f"{i:>12}  | {count}")
        
    ping_pong_score = distances.get(10, 0)
    background = sum(distances.values()) - ping_pong_score
    avg_bg = background / (max_distance - 1) if (max_distance - 1) > 0 else 1
    z_score = (ping_pong_score - avg_bg) / max(1, avg_bg**0.5) # Simplified Z-score approximation
    
    report.append("--------------------------")
    report.append(f"Z-score approximation for 10nt overlap: {z_score:.2f}")
    if z_score > 3:
        report.append("Result: Strong Ping-Pong signature detected!")
    else:
        report.append("Result: No significant Ping-Pong signature detected.")

    report_str = "\n".join(report)
    print("\n" + report_str)
    
    if output_file:
        with open(output_file, 'w') as out:
            out.write(report_str)
        print(f"\nReport saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Ping-Pong signature from BED file.")
    parser.add_argument("-i", "--input", required=True, help="Input BED file (e.g. from bedtools bamtobed)")
    parser.add_argument("-d", "--max_dist", type=int, default=30, help="Maximum distance to check (default: 30)")
    parser.add_argument("-o", "--output", help="Output text file for the report")
    
    args = parser.parse_args()
    calculate_ping_pong(args.input, args.max_dist, args.output)
