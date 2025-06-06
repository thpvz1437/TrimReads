"""
FASTQ File Parser Module

This module provides functionality for reading, writing, and processing FASTQ files,
including support for compressed (.gz) files and streaming processing.
"""

import gzip
import os
import sys
from collections import namedtuple
from typing import Iterator, List, Optional, TextIO, Union

# Define a named tuple for FASTQ records
FastqRecord = namedtuple('FastqRecord', ['header', 'sequence', 'plus', 'quality'])

class FastqParser:
    """
    A parser for FASTQ files with support for both plain text and gzip compressed files.
    
    Attributes:
        file_path (str): Path to the FASTQ file
        compressed (bool): True if file is gzip compressed
    """
    
    def __init__(self, file_path: str):
        """
        Initialize the FASTQ parser.
        
        Args:
            file_path: Path to the FASTQ file (supports .gz extension for compressed files)
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"FASTQ file not found: {file_path}")
            
        self.file_path = file_path
        self.compressed = file_path.endswith('.gz')
    
    def __enter__(self):
        """Context manager entry: open the file"""
        if self.compressed:
            self.file_handle = gzip.open(self.file_path, 'rt')
        else:
            self.file_handle = open(self.file_path, 'r')
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: close the file"""
        self.file_handle.close()
    
    def parse(self) -> Iterator[FastqRecord]:
        """
        Parse the FASTQ file and yield records one by one.
        
        Yields:
            FastqRecord: Parsed FASTQ record
            
        Raises:
            ValueError: If an invalid FASTQ record is encountered
        """
        line_count = 0
        record_lines = []
        
        for line in self.file_handle:
            line = line.strip()
            line_count += 1
            
            if line:  # Skip empty lines
                record_lines.append(line)
                
                # Every 4 lines make a complete FASTQ record
                if len(record_lines) == 4:
                    # Validate FASTQ record
                    if not record_lines[0].startswith('@'):
                        raise ValueError(
                            f"Invalid FASTQ header at line {line_count-3}: {record_lines[0]}"
                        f"\nExpected header to start with '@'"
                        f"\nFull record: {record_lines}"
                        )
                    
                    if record_lines[2] != '+':
                        raise ValueError(
                            f"Invalid separator line at line {line_count-1}: {record_lines[2]}"
                            f"\nExpected '+'"
                            f"\nFull record: {record_lines}"
                        )
                    
                    if len(record_lines[1]) != len(record_lines[3]):
                        raise ValueError(
                            f"Sequence/quality length mismatch at line {line_count}"
                            f"\nSequence length: {len(record_lines[1])}"
                            f"\nQuality length: {len(record_lines[3])}"
                            f"\nFull record: {record_lines}"
                        )
                    
                    # Create record and reset buffer
                    yield FastqRecord(
                        header=record_lines[0],
                        sequence=record_lines[1],
                        plus=record_lines[2],
                        quality=record_lines[3]
                    )
                    record_lines = []
    
    def read_all(self) -> List[FastqRecord]:
        """
        Read all records from the FASTQ file into memory.
        
        Returns:
            List of FastqRecord objects
            
        Warning: Not recommended for large files!
        """
        return list(self.parse())

class FastqWriter:
    """
    A writer for FASTQ files with support for both plain text and gzip compression.
    """
    
    def __init__(self, output_path: str, compress: bool = False):
        """
        Initialize the FASTQ writer.
        
        Args:
            output_path: Path to the output FASTQ file
            compress: If True, compress the output with gzip
        """
        self.output_path = output_path
        self.compress = compress or output_path.endswith('.gz')
        self.file_handle = None
    
    def __enter__(self):
        """Context manager entry: open the file"""
        if self.compress:
            self.file_handle = gzip.open(self.output_path, 'wt')
        else:
            self.file_handle = open(self.output_path, 'w')
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: close the file"""
        if self.file_handle:
            self.file_handle.close()
    
    def write_record(self, record: FastqRecord):
        """
        Write a single FASTQ record to the file.
        
        Args:
            record: FastqRecord to write
        """
        self.file_handle.write(f"{record.header}\n")
        self.file_handle.write(f"{record.sequence}\n")
        self.file_handle.write(f"{record.plus}\n")
        self.file_handle.write(f"{record.quality}\n")
    
    def write_records(self, records: List[FastqRecord]):
        """
        Write multiple FASTQ records to the file.
        
        Args:
            records: List of FastqRecord objects to write
        """
        for record in records:
            self.write_record(record)

def stream_fastq_records(
    input_file: Union[str, TextIO],
    output_file: Union[str, TextIO] = None,
    process_func: callable = None
) -> Iterator[FastqRecord]:
    """
    Stream FASTQ records from input to output with optional processing.
    
    Args:
        input_file: Input FASTQ file path or file-like object
        output_file: Output FASTQ file path or file-like object (optional)
        process_func: Function to process each record (should return a FastqRecord)
    
    Yields:
        FastqRecord: Processed FASTQ records
        
    Example:
        def trim_quality(record):
            # Trimming logic here
            return FastqRecord(...)
            
        for record in stream_fastq_records('input.fq', 'output.fq', trim_quality):
            # Do something with each record
    """
    input_handle = None
    output_handle = None
    
    # Handle input
    if isinstance(input_file, str):
        if input_file.endswith('.gz'):
            input_handle = gzip.open(input_file, 'rt')
        else:
            input_handle = open(input_file, 'r')
    else:
        input_handle = input_file
    
    # Handle output
    if output_file:
        if isinstance(output_file, str):
            if output_file.endswith('.gz'):
                output_handle = gzip.open(output_file, 'wt')
            else:
                output_handle = open(output_file, 'w')
        else:
            output_handle = output_file
    
    try:
        record_lines = []
        line_count = 0
        
        for line in input_handle:
            line = line.strip()
            line_count += 1
            
            if line:  # Skip empty lines
                record_lines.append(line)
                
                # Every 4 lines make a complete FASTQ record
                if len(record_lines) == 4:
                    # Create record
                    record = FastqRecord(
                        header=record_lines[0],
                        sequence=record_lines[1],
                        plus=record_lines[2],
                        quality=record_lines[3]
                    )
                    
                    # Process record if function provided
                    if process_func:
                        processed_record = process_func(record)
                    else:
                        processed_record = record
                    
                    # Write to output if provided
                    if output_handle:
                        output_handle.write(f"{processed_record.header}\n")
                        output_handle.write(f"{processed_record.sequence}\n")
                        output_handle.write(f"{processed_record.plus}\n")
                        output_handle.write(f"{processed_record.quality}\n")
                    
                    # Yield processed record
                    yield processed_record
                    
                    # Reset buffer
                    record_lines = []
    
    finally:
        # Clean up handles if we opened them
        if isinstance(input_file, str) and input_handle:
            input_handle.close()
        if isinstance(output_file, str) and output_handle:
            output_handle.close()

def validate_fastq(file_path: str) -> bool:
    """
    Validate the structure of a FASTQ file.
    
    Args:
        file_path: Path to the FASTQ file
        
    Returns:
        True if valid, False otherwise
    """
    try:
        with FastqParser(file_path) as parser:
            for _ in parser.parse():
                pass  # Just iterate through all records
        return True
    except (ValueError, UnicodeDecodeError) as e:
        print(f"Validation error: {str(e)}", file=sys.stderr)
        return False

def extract_quality_scores(record: FastqRecord) -> List[int]:
    """
    Extract quality scores from a FASTQ record (Sanger/Phred+33 encoding).
    
    Args:
        record: FastqRecord object
        
    Returns:
        List of integer quality scores
    """
    return [ord(char) - 33 for char in record.quality]

def calculate_average_quality(record: FastqRecord) -> float:
    """
    Calculate the average quality score for a FASTQ record.
    
    Args:
        record: FastqRecord object
        
    Returns:
        Average quality score
    """
    scores = extract_quality_scores(record)
    return sum(scores) / len(scores) if scores else 0.0

def filter_records(
    input_file: str,
    output_file: str,
    min_quality: float = 20.0,
    min_length: int = 50,
    max_length: int = None
) -> dict:
    """
    Filter FASTQ records based on quality and length criteria.
    
    Args:
        input_file: Input FASTQ file path
        output_file: Output FASTQ file path
        min_quality: Minimum average quality to keep
        min_length: Minimum sequence length to keep
        max_length: Maximum sequence length to keep (optional)
        
    Returns:
        Dictionary with filtering statistics
    """
    stats = {
        'total': 0,
        'passed': 0,
        'low_quality': 0,
        'too_short': 0,
        'too_long': 0
    }
    
    with FastqParser(input_file) as parser, FastqWriter(output_file) as writer:
        for record in parser.parse():
            stats['total'] += 1
            
            # Calculate average quality
            avg_quality = calculate_average_quality(record)
            seq_length = len(record.sequence)
            
            # Apply filters
            if avg_quality < min_quality:
                stats['low_quality'] += 1
                continue
                
            if seq_length < min_length:
                stats['too_short'] += 1
                continue
                
            if max_length is not None and seq_length > max_length:
                stats['too_long'] += 1
                continue
                
            # Write passing record
            writer.write_record(record)
            stats['passed'] += 1
    
    return stats

if __name__ == '__main__':
    # Example usage
    import argparse
    
    parser = argparse.ArgumentParser(description='FASTQ File Utilities')
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate a FASTQ file')
    validate_parser.add_argument('file', help='FASTQ file to validate')
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter a FASTQ file')
    filter_parser.add_argument('input', help='Input FASTQ file')
    filter_parser.add_argument('output', help='Output FASTQ file')
    filter_parser.add_argument('--min_quality', type=float, default=20.0,
                              help='Minimum average quality score')
    filter_parser.add_argument('--min_length', type=int, default=50,
                              help='Minimum sequence length')
    filter_parser.add_argument('--max_length', type=int, 
                              help='Maximum sequence length')
    
    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Get statistics for a FASTQ file')
    stats_parser.add_argument('file', help='FASTQ file to analyze')
    
    args = parser.parse_args()
    
    if args.command == 'validate':
        result = validate_fastq(args.file)
        print(f"File {args.file} is {'VALID' if result else 'INVALID'}")
        sys.exit(0 if result else 1)
    
    elif args.command == 'filter':
        print(f"Filtering {args.input} to {args.output}")
        stats = filter_records(
            args.input,
            args.output,
            min_quality=args.min_quality,
            min_length=args.min_length,
            max_length=args.max_length
        )
        print("\nFiltering Statistics:")
        print(f"Total records: {stats['total']}")
        print(f"Passed records: {stats['passed']} ({stats['passed']/stats['total']:.1%})")
        print(f"Low quality: {stats['low_quality']}")
        print(f"Too short: {stats['too_short']}")
        if args.max_length:
            print(f"Too long: {stats['too_long']}")
    
    elif args.command == 'stats':
        total_records = 0
        total_bases = 0
        min_length = float('inf')
        max_length = 0
        quality_sum = 0.0
        
        with FastqParser(args.file) as parser:
            for record in parser.parse():
                total_records += 1
                seq_len = len(record.sequence)
                total_bases += seq_len
                min_length = min(min_length, seq_len)
                max_length = max(max_length, seq_len)
                quality_sum += calculate_average_quality(record)
        
        avg_quality = quality_sum / total_records if total_records else 0
        avg_length = total_bases / total_records if total_records else 0
        
        print(f"FASTQ File Statistics for {args.file}:")
        print(f"Total records: {total_records}")
        print(f"Total bases: {total_bases}")
        print(f"Average sequence length: {avg_length:.1f} bp")
        print(f"Min sequence length: {min_length} bp")
        print(f"Max sequence length: {max_length} bp")
        print(f"Average quality score: {avg_quality:.2f}")
    
    else:
        parser.print_help()
        sys.exit(1)