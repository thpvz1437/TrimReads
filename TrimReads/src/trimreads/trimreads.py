#!/usr/bin/env python3
"""
TrimReads - High-throughput sequencing data quality trimming tool
Performs base-by-base and window-based quality trimming of FASTQ files
"""

import argparse
from collections import deque
import gzip
import os
import sys
from typing import Tuple, Optional

def phred_to_score(char: str) -> int:
    """Convert Phred character to quality score (Sanger encoding)"""
    return ord(char) - 33

def score_to_phred(score: int) -> str:
    """Convert quality score to Phred character (Sanger encoding)"""
    return chr(score + 33)

def base_trim(sequence: str, quality: str, threshold: int) -> Tuple[str, str]:
    """
    Trim low-quality bases from both ends of a read
    
    Args:
        sequence: Nucleotide sequence
        quality: Quality score string
        threshold: Minimum quality score to keep
    
    Returns:
        Tuple of (trimmed_sequence, trimmed_quality)
    """
    # Convert to list for processing
    seq_chars = list(sequence)
    qual_chars = list(quality)
    
    # Find left trim position
    start_index = 0
    for i in range(len(seq_chars)):
        if phred_to_score(qual_chars[i]) >= threshold:
            start_index = i
            break
    
    # Find right trim position
    end_index = len(seq_chars) - 1
    for i in range(len(seq_chars)-1, -1, -1):
        if phred_to_score(qual_chars[i]) >= threshold:
            end_index = i
            break
    
    # Discard read if completely low quality
    if start_index > end_index:
        return "", ""
    
    # Return trimmed sequences
    trimmed_seq = "".join(seq_chars[start_index:end_index+1])
    trimmed_qual = "".join(qual_chars[start_index:end_index+1])
    return trimmed_seq, trimmed_qual

def window_trim(
    sequence: str, 
    quality: str, 
    window_size: int, 
    threshold: int
) -> Tuple[str, str]:
    """
    Trim reads using sliding window approach
    
    Args:
        sequence: Nucleotide sequence
        quality: Quality score string
        window_size: Size of sliding window
        threshold: Minimum average quality to keep
    
    Returns:
        Tuple of (trimmed_sequence, trimmed_quality)
    """
    # Handle edge cases
    if len(sequence) < window_size:
        return sequence, quality
    
    # Convert quality to scores
    scores = [phred_to_score(q) for q in quality]
    
    # Find left trim position
    start_index = 0
    window = deque(scores[:window_size], maxlen=window_size)
    window_avg = sum(window) / window_size
    
    if window_avg >= threshold:
        start_index = 0
    else:
        found = False
        for i in range(window_size, len(scores)):
            # Update sliding window
            window.append(scores[i])
            window_avg = sum(window) / window_size
            
            # Check if window meets threshold
            if window_avg >= threshold:
                start_index = i - window_size + 1
                found = True
                break
        
        # If no suitable window found, discard read
        if not found:
            return "", ""
    
    # Find right trim position
    end_index = len(sequence) - 1
    window = deque(scores[-window_size:], maxlen=window_size)
    window_avg = sum(window) / window_size
    
    if window_avg >= threshold:
        end_index = len(sequence) - 1
    else:
        found = False
        for i in range(len(scores)-window_size-1, -1, -1):
            # Update sliding window
            window.appendleft(scores[i])
            window_avg = sum(window) / window_size
            
            # Check if window meets threshold
            if window_avg >= threshold:
                end_index = i + window_size - 1
                found = True
                break
        
        # If no suitable window found, discard read
        if not found:
            return "", ""
    
    # Return trimmed sequences
    trimmed_seq = sequence[start_index:end_index+1]
    trimmed_qual = quality[start_index:end_index+1]
    return trimmed_seq, trimmed_qual

def process_fastq(
    input_file: str,
    output_file: str,
    base_threshold: Optional[int] = None,
    window_size: Optional[int] = None,
    window_threshold: Optional[int] = None,
    min_length: int = 30
) -> dict:
    """
    Process FASTQ file with quality trimming
    
    Args:
        input_file: Path to input FASTQ file
        output_file: Path to output FASTQ file
        base_threshold: Quality threshold for base trimming
        window_size: Window size for window-based trimming
        window_threshold: Quality threshold for window-based trimming
        min_length: Minimum read length to keep
    
    Returns:
        Dictionary with processing statistics
    """
    # Initialize counters
    stats = {
        'total_reads': 0,
        'passed_reads': 0,
        'discarded_reads': 0,
        'min_length': min_length,
        'base_threshold': base_threshold,
        'window_size': window_size,
        'window_threshold': window_threshold
    }
    
    # Open input file (handles gzipped files)
    open_func = gzip.open if input_file.endswith('.gz') else open
    in_mode = 'rt' if input_file.endswith('.gz') else 'r'
    
    # Open output file
    out_mode = 'wt' if output_file.endswith('.gz') else 'w'
    output_handle = gzip.open(output_file, out_mode) if output_file.endswith('.gz') else open(output_file, out_mode)
    
    with open_func(input_file, in_mode) as f:
        while True:
            # Read FASTQ record (4 lines per read)
            header = f.readline().rstrip()
            if not header:
                break  # End of file
                
            sequence = f.readline().rstrip()
            plus_line = f.readline().rstrip()
            quality = f.readline().rstrip()
            
            stats['total_reads'] += 1
            
            # Validate FASTQ record
            if not header.startswith('@'):
                sys.stderr.write(f"Warning: Invalid header at read {stats['total_reads']}: {header}\n")
                continue
                
            if len(sequence) != len(quality):
                sys.stderr.write(f"Warning: Length mismatch at read {stats['total_reads']}: "
                                 f"seq_len={len(sequence)}, qual_len={len(quality)}\n")
                continue
            
            # Apply base trimming if requested
            if base_threshold is not None:
                sequence, quality = base_trim(sequence, quality, base_threshold)
                if not sequence:  # Read discarded
                    stats['discarded_reads'] += 1
                    continue
            
            # Apply window trimming if requested
            if window_size is not None and window_threshold is not None:
                sequence, quality = window_trim(sequence, quality, window_size, window_threshold)
                if not sequence:  # Read discarded
                    stats['discarded_reads'] += 1
                    continue
            
            # Discard reads below minimum length
            if len(sequence) < min_length:
                stats['discarded_reads'] += 1
                continue
            
            # Write trimmed read to output
            output_handle.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
            stats['passed_reads'] += 1
    
    output_handle.close()
    return stats

def main():
    """Command-line interface for TrimReads"""
    parser = argparse.ArgumentParser(
        description='Quality trimming for FASTQ files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTQ file (gzip compressed allowed)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output trimmed FASTQ file (use .gz extension for compression)')
    parser.add_argument('--base_threshold', type=int, default=None,
                        help='Quality threshold for base-by-base trimming')
    parser.add_argument('--window_size', type=int, default=None,
                        help='Window size for window-based trimming')
    parser.add_argument('--window_threshold', type=int, default=None,
                        help='Quality threshold for window-based trimming')
    parser.add_argument('--min_length', type=int, default=30,
                        help='Minimum read length to keep after trimming')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not os.path.exists(args.input):
        sys.exit(f"Error: Input file not found: {args.input}")
    
    if args.window_size is not None and args.window_threshold is None:
        sys.exit("Error: Must specify --window_threshold when using --window_size")
    
    if args.window_threshold is not None and args.window_size is None:
        sys.exit("Error: Must specify --window_size when using --window_threshold")
    
    # Process the FASTQ file
    print("Starting quality trimming...")
    stats = process_fastq(
        input_file=args.input,
        output_file=args.output,
        base_threshold=args.base_threshold,
        window_size=args.window_size,
        window_threshold=args.window_threshold,
        min_length=args.min_length
    )
    
    # Print summary statistics
    print("\n=== Trimming Summary ===")
    print(f"Total reads processed: {stats['total_reads']}")
    print(f"Reads passing filters: {stats['passed_reads']} "
          f"({stats['passed_reads']/stats['total_reads']:.2%})")
    print(f"Reads discarded: {stats['discarded_reads']} "
          f"({stats['discarded_reads']/stats['total_reads']:.2%})")
    
    if args.base_threshold:
        print(f"Base-by-base threshold: Q{args.base_threshold}")
    
    if args.window_size and args.window_threshold:
        print(f"Window-based trimming: {args.window_size}bp window, Q{args.window_threshold} threshold")
    
    print(f"Minimum length kept: {args.min_length}bp")
    print(f"Output saved to: {args.output}")

if __name__ == '__main__':
    main()