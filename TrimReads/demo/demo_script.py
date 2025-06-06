#!/usr/bin/env python3
"""
TrimReads Demo Script

This script demonstrates the functionality of the TrimReads package
using both simulated and real-world sequencing data.
"""

import os
import gzip
import random
import tempfile
import subprocess
from trimreads import trimreads, fastq_parser

def print_header(title):
    """Print a formatted header"""
    print("\n" + "=" * 80)
    print(f" {title} ".center(80, "="))
    print("=" * 80 + "\n")

def create_simulated_fastq(file_path, num_reads=1000, read_length=150):
    """
    Create a simulated FASTQ file with quality degradation at ends
    
    Args:
        file_path: Path to output FASTQ file
        num_reads: Number of reads to generate
        read_length: Length of each read
    """
    bases = "ACGT"
    
    # Open file (handle gzip compression)
    if file_path.endswith('.gz'):
        f = gzip.open(file_path, 'wt')
    else:
        f = open(file_path, 'w')
    
    with f:
        for i in range(num_reads):
            # Generate random sequence
            sequence = ''.join(random.choices(bases, k=read_length))
            
            # Generate quality scores:
            # - High quality in the middle (Q30-40)
            # - Degraded at ends (Q0-20)
            quality = []
            for pos in range(read_length):
                # Calculate position factor (0 at ends, 1 in middle)
                pos_factor = min(pos, read_length - pos - 1) / (read_length / 2)
                pos_factor = min(1.0, pos_factor)  # Cap at 1.0
                
                # Base quality: high in middle, low at ends
                base_quality = int(20 + 20 * pos_factor)
                quality.append(trimreads.score_to_phred(base_quality))
            quality = ''.join(quality)
            
            # Write FASTQ record
            f.write(f"@SIMULATED_READ_{i}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality}\n")
    
    print(f"Created simulated FASTQ file: {file_path}")
    print(f"  Reads: {num_reads}, Length: {read_length} bp")

def run_trimming_demo(input_file, output_file, method="base", threshold=25, window_size=10):
    """
    Run a trimming demo and show results
    
    Args:
        input_file: Input FASTQ file
        output_file: Output trimmed FASTQ file
        method: "base" or "window" trimming
        threshold: Quality threshold
        window_size: Window size (for window method)
    """
    # Run trimming based on method
    if method == "base":
        print(f"Running base-by-base trimming with threshold Q{threshold}...")
        stats = trimreads.process_fastq(
            input_file,
            output_file,
            base_threshold=threshold,
            min_length=30
        )
    elif method == "window":
        print(f"Running window-based trimming (window={window_size}bp, threshold Q{threshold})...")
        stats = trimreads.process_fastq(
            input_file,
            output_file,
            window_size=window_size,
            window_threshold=threshold,
            min_length=30
        )
    else:
        print(f"Running combined base and window trimming...")
        stats = trimreads.process_fastq(
            input_file,
            output_file,
            base_threshold=threshold,
            window_size=window_size,
            window_threshold=threshold,
            min_length=30
        )
    
    # Print summary statistics
    print("\nTrimming Summary:")
    print(f"  Total reads processed: {stats['total_reads']}")
    print(f"  Reads passing filters: {stats['passed_reads']} "
          f"({stats['passed_reads']/stats['total_reads']:.2%})")
    print(f"  Reads discarded: {stats['discarded_reads']} "
          f"({stats['discarded_reads']/stats['total_reads']:.2%})")
    
    # Calculate quality improvement
    input_stats = fastq_parser.calculate_file_stats(input_file)
    output_stats = fastq_parser.calculate_file_stats(output_file)
    
    print("\nQuality Improvement:")
    print(f"  Average quality before: {input_stats['avg_quality']:.2f}")
    print(f"  Average quality after:  {output_stats['avg_quality']:.2f}")
    print(f"  Bases with Q<20 before: {input_stats['low_quality_bases']} "
          f"({input_stats['low_quality_percent']:.2%})")
    print(f"  Bases with Q<20 after:  {output_stats['low_quality_bases']} "
          f"({output_stats['low_quality_percent']:.2%})")
    
    return stats

def visualize_quality_profiles(input_file, output_file):
    """
    Visualize quality profiles before and after trimming
    
    Args:
        input_file: Original FASTQ file
        output_file: Trimmed FASTQ file
    """
    # Get random samples
    input_sample = fastq_parser.get_random_records(input_file, 5)
    output_sample = fastq_parser.get_random_records(output_file, 5)
    
    print("\nQuality Profile Visualization:")
    print("Position: 1 2 3 4 5 6 7 8 9 10 ...")
    print("          | | | | | | | | | |")
    
    for i, (in_rec, out_rec) in enumerate(zip(input_sample, output_sample)):
        # Only show first 50 positions for clarity
        in_qual = in_rec.quality[:50]
        out_qual = out_rec.quality[:50]
        
        # Convert to scores for visualization
        in_scores = [trimreads.phred_to_score(q) for q in in_qual]
        out_scores = [trimreads.phred_to_score(q) for q in out_qual]
        
        # Create visual representation
        in_visual = ''.join('▅' if s >= 30 else '▄' if s >= 20 else '▃' if s >= 10 else '▁' for s in in_scores)
        out_visual = ''.join('▅' if s >= 30 else '▄' if s >= 20 else '▃' if s >= 10 else '▁' for s in out_scores)
        
        print(f"\nRead {i+1} - Before Trimming:")
        print(in_visual)
        print(f"Read {i+1} - After Trimming:")
        print(out_visual)
    
    print("\nQuality Key:")
    print("▅ = Q≥30 (Excellent)  ▄ = Q20-29 (Good)  ▃ = Q10-19 (Marginal)  ▁ = Q<10 (Poor)")

def download_real_data(output_dir):
    """
    Download real sequencing data from ENA (European Nucleotide Archive)
    
    Args:
        output_dir: Directory to save downloaded files
        
    Returns:
        Path to downloaded FASTQ file
    """
    # Use a small COVID-19 sequencing dataset from ENA
    accession = "SRR15338479"
    url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{accession}/{accession}"
    
    # Create output path
    fastq_gz = os.path.join(output_dir, f"{accession}.fastq.gz")
    
    print(f"Downloading real sequencing data ({accession}) from ENA...")
    
    # Use fasterq-dump if available, otherwise use curl
    try:
        # Try using fasterq-dump (from SRA Toolkit)
        subprocess.run(["fasterq-dump", "--progress", accession], check=True, cwd=output_dir)
        # Compress to gzip
        subprocess.run(["gzip", f"{accession}.fastq"], check=True, cwd=output_dir)
        print("Download completed using fasterq-dump")
    except (FileNotFoundError, subprocess.CalledProcessError):
        # Fallback to direct download
        print("fasterq-dump not available, using direct download...")
        try:
            subprocess.run(["curl", "-o", fastq_gz, url], check=True)
            print("Direct download completed")
        except subprocess.CalledProcessError:
            print("Failed to download real data. Using simulated data for demo.")
            return None
    
    return fastq_gz

def main():
    """Main demo function"""
    # Create temporary directory for files
    with tempfile.TemporaryDirectory() as tmpdir:
        print_header("TrimReads - High-Throughput Sequencing Data Quality Trimming")
        print("This demo showcases the functionality of the TrimReads package.")
        
        # Part 1: Simulated Data Demo
        print_header("Part 1: Simulated Data Demo")
        sim_fastq = os.path.join(tmpdir, "simulated.fastq")
        create_simulated_fastq(sim_fastq, num_reads=10000, read_length=150)
        
        # Show file stats before trimming
        sim_stats = fastq_parser.calculate_file_stats(sim_fastq)
        print("\nSimulated Data Statistics Before Trimming:")
        print(f"  Total reads: {sim_stats['total_reads']}")
        print(f"  Total bases: {sim_stats['total_bases']}")
        print(f"  Average read length: {sim_stats['avg_length']:.1f} bp")
        print(f"  Average quality: {sim_stats['avg_quality']:.2f}")
        print(f"  Bases with Q<20: {sim_stats['low_quality_bases']} "
              f"({sim_stats['low_quality_percent']:.2%})")
        
        # Demo 1: Base-by-base trimming
        print_header("Demo 1: Base-by-Base Trimming")
        base_output = os.path.join(tmpdir, "base_trimmed.fastq")
        run_trimming_demo(sim_fastq, base_output, method="base", threshold=25)
        
        # Demo 2: Window-based trimming
        print_header("Demo 2: Window-Based Trimming")
        window_output = os.path.join(tmpdir, "window_trimmed.fastq")
        run_trimming_demo(sim_fastq, window_output, method="window", threshold=20, window_size=10)
        
        # Demo 3: Combined trimming
        print_header("Demo 3: Combined Base and Window Trimming")
        combined_output = os.path.join(tmpdir, "combined_trimmed.fastq")
        run_trimming_demo(sim_fastq, combined_output, method="combined", threshold=25, window_size=10)
        
        # Visualize quality profiles
        print_header("Quality Profile Visualization")
        visualize_quality_profiles(sim_fastq, combined_output)
        
        # Part 2: Real-World Data Demo
        print_header("Part 2: Real-World Data Demo")
        real_fastq = download_real_data(tmpdir)
        
        if real_fastq and os.path.exists(real_fastq):
            # Show file stats before trimming
            real_stats = fastq_parser.calculate_file_stats(real_fastq)
            print("\nReal Data Statistics Before Trimming:")
            print(f"  Total reads: {real_stats['total_reads']}")
            print(f"  Total bases: {real_stats['total_bases']}")
            print(f"  Average read length: {real_stats['avg_length']:.1f} bp")
            print(f"  Average quality: {real_stats['avg_quality']:.2f}")
            print(f"  Bases with Q<20: {real_stats['low_quality_bases']} "
                  f"({real_stats['low_quality_percent']:.2%})")
            
            # Trim real data
            print_header("Trimming Real Sequencing Data")
            real_output = os.path.join(tmpdir, "real_trimmed.fastq.gz")
            run_trimming_demo(
                real_fastq, 
                real_output, 
                method="combined", 
                threshold=25, 
                window_size=10
            )
        else:
            print("Skipping real data demo due to download issues.")
        
        print_header("Demo Complete")
        print("All demo files were created in temporary directory and have been cleaned up.")
        print("To use TrimReads with your own data:")
        print("  trimreads -i input.fastq -o output.fastq --base_threshold 25 --window_size 10 --window_threshold 20")
        print("\nFor more information, see the package documentation.")

if __name__ == "__main__":
    main()