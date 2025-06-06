"""
TrimReads - High-throughput sequencing data quality trimming package

This package provides tools for quality trimming of FASTQ files from high-throughput sequencing data.
It includes both base-by-base and window-based trimming algorithms, as well as a comprehensive FASTQ
parser and utilities for quality control.

Key Components:
- trimreads: Core trimming algorithms and command-line interface
- fastq_parser: FASTQ file parsing, validation, and quality analysis

Usage:
    from trimreads import trimreads, fastq_parser
    
    # Perform trimming
    stats = trimreads.process_fastq('input.fastq', 'output.fastq', 
                                   base_threshold=25, 
                                   window_size=10, 
                                   window_threshold=20)
    
    # Parse and analyze FASTQ files
    with fastq_parser.FastqParser('input.fastq') as parser:
        for record in parser.parse():
            avg_quality = fastq_parser.calculate_average_quality(record)
            print(f"Read: {record.header}, Avg Quality: {avg_quality:.2f}")

Command-Line Tools:
    trimreads: Quality trimming for FASTQ files
    fastq_utils: FASTQ validation, filtering, and statistics
"""

from . import trimreads
from . import fastq_parser

# Import core functions for direct access
from .trimreads import (
    base_trim,
    window_trim,
    process_fastq,
    phred_to_score,
    score_to_phred
)

from .fastq_parser import (
    FastqRecord,
    FastqParser,
    FastqWriter,
    stream_fastq_records,
    validate_fastq,
    extract_quality_scores,
    calculate_average_quality,
    filter_records
)

# Package metadata
__version__ = "1.0.0"
__author__ = "Bioinformatics Team"
__email__ = "bioinfo@example.com"
__license__ = "MIT"
__all__ = [
    'trimreads',
    'fastq_parser',
    'base_trim',
    'window_trim',
    'process_fastq',
    'phred_to_score',
    'score_to_phred',
    'FastqRecord',
    'FastqParser',
    'FastqWriter',
    'stream_fastq_records',
    'validate_fastq',
    'extract_quality_scores',
    'calculate_average_quality',
    'filter_records'
]

# Package initialization message
print(f"Loaded TrimReads package v{__version__}")
print("Available modules: trimreads, fastq_parser")
print("Use help(trimreads) or help(fastq_parser) for more information")