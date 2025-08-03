"""
conda activate ngs

Split FASTQ data into multiple non-overlapping replicates
Supports both single-end and paired-end sequencing data
Usage: 
    Single-end: python split_fastq.py sample.fq.gz --split_ratio 0.5 0.5
    Paired-end: python split_fastq.py sample_1.fq.gz sample_2.fq.gz --split_ratio 0.33 0.33 0.34
"""

import sys
import gzip
import random
import argparse
import re
from pathlib import Path
import io

def detect_file_format(filename):
    """Detect if file is gzipped and return appropriate open function"""
    if filename.endswith('.gz'):
        return gzip.open, 'rb'
    else:
        return open, 'r'

def parse_fastq_batch(file_handle, batch_size=10000, is_gzipped=True):
    """Parse multiple FASTQ records in a batch for better performance"""
    records = []
    
    for _ in range(batch_size):
        if is_gzipped:
            header = file_handle.readline()
            if not header:
                break
            header = header.decode('utf-8').strip()
            sequence = file_handle.readline().decode('utf-8').strip()
            plus = file_handle.readline().decode('utf-8').strip()
            quality = file_handle.readline().decode('utf-8').strip()
        else:
            header = file_handle.readline().strip()
            if not header:
                break
            sequence = file_handle.readline().strip()
            plus = file_handle.readline().strip()
            quality = file_handle.readline().strip()
        
        records.append((header, sequence, plus, quality))
    
    return records

def write_fastq_batch(file_handle, records, is_gzipped=True):
    """Write multiple FASTQ records in a batch for better performance"""
    if is_gzipped:
        # Build content as string first, then encode once
        content = '\n'.join(f"{header}\n{sequence}\n{plus}\n{quality}" 
                           for header, sequence, plus, quality in records)
        content += '\n'
        file_handle.write(content.encode('utf-8'))
    else:
        # For text files, use join for better performance
        content = '\n'.join(f"{header}\n{sequence}\n{plus}\n{quality}" 
                           for header, sequence, plus, quality in records)
        content += '\n'
        file_handle.write(content)

def get_base_name(filename):
    """Extract base name from various FASTQ file formats"""
    path = Path(filename)
    name = path.name
    
    # Remove .gz extension if present
    if name.endswith('.gz'):
        name = name[:-3]
    
    # Remove common FASTQ extensions
    extensions = ['.fastq', '.fq']
    for ext in extensions:
        if name.endswith(ext):
            name = name[:-len(ext)]
            break
    
    return name

def get_output_filename(input_filename, replicate_num, pair_num=None):
    """Generate output filename based on input filename and replicate number"""
    base_name = get_base_name(input_filename)
    
    # Handle paired-end naming
    if pair_num is not None:
        # Remove existing pair indicators (_1, _2, _R1, _R2, etc.)
        patterns = [r'_[12]$', r'_R[12]$', r'\.R[12]$']
        for pattern in patterns:
            base_name = re.sub(pattern, '', base_name)
        
        output_name = f"{base_name}_rep{replicate_num}_{pair_num}"
    else:
        output_name = f"{base_name}_rep{replicate_num}"
    
    # Determine output extension
    if input_filename.endswith('.gz'):
        output_name += '.fq.gz'
    else:
        output_name += '.fq'
    
    return output_name

def generate_replicate_assignments(num_reads, split_ratios, seed=42):
    """Pre-generate all replicate assignments for better performance"""
    random.seed(seed)
    
    # Calculate cumulative ratios once
    cumulative_ratios = []
    cumulative = 0
    for ratio in split_ratios:
        cumulative += ratio
        cumulative_ratios.append(cumulative)
    
    # Generate all random values at once
    random_values = [random.random() for _ in range(num_reads)]
    
    # Assign replicates efficiently
    assignments = []
    for rand_val in random_values:
        for i, cum_ratio in enumerate(cumulative_ratios):
            if rand_val <= cum_ratio:
                assignments.append(i)
                break
        else:
            assignments.append(len(split_ratios) - 1)
    
    return assignments

def count_reads_in_file(file_path):
    """Count total number of reads in FASTQ file"""
    open_func, mode = detect_file_format(file_path)
    is_gzipped = file_path.endswith('.gz')
    
    count = 0
    with open_func(file_path, mode) as f:
        if is_gzipped:
            # For gzipped files, count lines and divide by 4
            for line in f:
                count += 1
        else:
            # For text files, count lines and divide by 4
            for line in f:
                count += 1
    
    return count // 4

def split_single_end_fastq(file_path, split_ratios, seed=42, batch_size=20000):
    """Split single-end FASTQ file into multiple replicates"""
    
    num_replicates = len(split_ratios)
    
    # Detect file format
    open_func, mode = detect_file_format(file_path)
    is_gzipped = file_path.endswith('.gz')
    
    # Generate output filenames
    output_files = []
    for i in range(num_replicates):
        output_files.append(get_output_filename(file_path, i + 1))
    
    print(f"Splitting single-end file: {file_path}")
    print(f"Split ratios: {split_ratios}")
    for i, out_file in enumerate(output_files):
        print(f"Replicate {i + 1}: {out_file} ({split_ratios[i]*100:.1f}%)")
    
    # Count total reads first (optional, for progress tracking)
    print("Counting reads...")
    total_reads = count_reads_in_file(file_path)
    print(f"Total reads to process: {total_reads:,}")
    
    # Pre-generate replicate assignments
    print("Generating replicate assignments...")
    assignments = generate_replicate_assignments(total_reads, split_ratios, seed)
    
    # Determine output file opening mode
    if is_gzipped:
        out_mode = 'wb'
        out_open = gzip.open
    else:
        out_mode = 'w'
        out_open = open
    
    # Open all output files
    output_handles = []
    for out_file in output_files:
        output_handles.append(out_open(out_file, out_mode))
    
    # Initialize batch containers for each replicate
    replicate_batches = [[] for _ in range(num_replicates)]
    
    try:
        # Process the file in batches
        with open_func(file_path, mode) as infile:
            read_count = 0
            replicate_counts = [0] * num_replicates
            
            while read_count < total_reads:
                # Read batch
                batch_records = parse_fastq_batch(infile, batch_size, is_gzipped)
                
                if not batch_records:
                    break
                
                # Distribute records to replicate batches
                for i, record in enumerate(batch_records):
                    replicate_idx = assignments[read_count + i]
                    replicate_batches[replicate_idx].append(record)
                    replicate_counts[replicate_idx] += 1
                
                read_count += len(batch_records)
                
                # Write batches when they get large enough
                for rep_idx in range(num_replicates):
                    if len(replicate_batches[rep_idx]) >= batch_size:
                        write_fastq_batch(output_handles[rep_idx], 
                                        replicate_batches[rep_idx], is_gzipped)
                        replicate_batches[rep_idx] = []
                
                # Progress report
                if read_count % batch_size == 0:
                    print(f"Processed {read_count:,} reads ({read_count/total_reads*100:.1f}%)")
            
            # Write remaining records in batches
            for rep_idx in range(num_replicates):
                if replicate_batches[rep_idx]:
                    write_fastq_batch(output_handles[rep_idx], 
                                    replicate_batches[rep_idx], is_gzipped)
    
    finally:
        # Close all output files
        for handle in output_handles:
            handle.close()
    
    print(f"\nSplitting completed:")
    print(f"Total reads processed: {read_count:,}")
    for i, count in enumerate(replicate_counts):
        percentage = count / read_count * 100 if read_count > 0 else 0
        print(f"Replicate {i + 1}: {count:,} reads ({percentage:.1f}%)")

def split_paired_end_fastq(file1_path, file2_path, split_ratios, seed=42, batch_size=20000):
    """Split paired-end FASTQ files into multiple replicates"""
    
    num_replicates = len(split_ratios)
    
    # Detect file formats
    open_func1, mode1 = detect_file_format(file1_path)
    open_func2, mode2 = detect_file_format(file2_path)
    is_gzipped1 = file1_path.endswith('.gz')
    is_gzipped2 = file2_path.endswith('.gz')
    
    # Generate output filenames
    output_files1 = []
    output_files2 = []
    for i in range(num_replicates):
        output_files1.append(get_output_filename(file1_path, i + 1, 1))
        output_files2.append(get_output_filename(file2_path, i + 1, 2))
    
    print(f"Splitting paired-end files: {file1_path}, {file2_path}")
    print(f"Split ratios: {split_ratios}")
    for i in range(num_replicates):
        print(f"Replicate {i + 1}: {output_files1[i]}, {output_files2[i]} ({split_ratios[i]*100:.1f}%)")
    
    # Count total read pairs
    print("Counting read pairs...")
    total_pairs = count_reads_in_file(file1_path)
    print(f"Total read pairs to process: {total_pairs:,}")
    
    # Pre-generate replicate assignments
    print("Generating replicate assignments...")
    assignments = generate_replicate_assignments(total_pairs, split_ratios, seed)
    
    # Determine output file opening modes
    out_mode1 = 'wb' if is_gzipped1 else 'w'
    out_mode2 = 'wb' if is_gzipped2 else 'w'
    out_open1 = gzip.open if is_gzipped1 else open
    out_open2 = gzip.open if is_gzipped2 else open
    
    # Open all output files
    output_handles1 = []
    output_handles2 = []
    for i in range(num_replicates):
        output_handles1.append(out_open1(output_files1[i], out_mode1))
        output_handles2.append(out_open2(output_files2[i], out_mode2))
    
    # Initialize batch containers for each replicate
    replicate_batches1 = [[] for _ in range(num_replicates)]
    replicate_batches2 = [[] for _ in range(num_replicates)]
    
    try:
        # Process the files in batches
        with open_func1(file1_path, mode1) as f1, open_func2(file2_path, mode2) as f2:
            read_count = 0
            replicate_counts = [0] * num_replicates
            
            while read_count < total_pairs:
                # Read batches from both files
                batch_records1 = parse_fastq_batch(f1, batch_size, is_gzipped1)
                batch_records2 = parse_fastq_batch(f2, batch_size, is_gzipped2)
                
                if not batch_records1 or not batch_records2:
                    break
                
                # Ensure both batches have the same size
                min_batch_size = min(len(batch_records1), len(batch_records2))
                
                # Distribute record pairs to replicate batches
                for i in range(min_batch_size):
                    replicate_idx = assignments[read_count + i]
                    replicate_batches1[replicate_idx].append(batch_records1[i])
                    replicate_batches2[replicate_idx].append(batch_records2[i])
                    replicate_counts[replicate_idx] += 1
                
                read_count += min_batch_size
                
                # Write batches when they get large enough
                for rep_idx in range(num_replicates):
                    if len(replicate_batches1[rep_idx]) >= batch_size:
                        write_fastq_batch(output_handles1[rep_idx], 
                                        replicate_batches1[rep_idx], is_gzipped1)
                        write_fastq_batch(output_handles2[rep_idx], 
                                        replicate_batches2[rep_idx], is_gzipped2)
                        replicate_batches1[rep_idx] = []
                        replicate_batches2[rep_idx] = []
                
                # Progress report
                if read_count % batch_size == 0:
                    print(f"Processed {read_count:,} read pairs ({read_count/total_pairs*100:.1f}%)")
            
            # Write remaining records in batches
            for rep_idx in range(num_replicates):
                if replicate_batches1[rep_idx]:
                    write_fastq_batch(output_handles1[rep_idx], 
                                    replicate_batches1[rep_idx], is_gzipped1)
                    write_fastq_batch(output_handles2[rep_idx], 
                                    replicate_batches2[rep_idx], is_gzipped2)
    
    finally:
        # Close all output files
        for handle in output_handles1:
            handle.close()
        for handle in output_handles2:
            handle.close()
    
    print(f"\nSplitting completed:")
    print(f"Total read pairs processed: {read_count:,}")
    for i, count in enumerate(replicate_counts):
        percentage = count / read_count * 100 if read_count > 0 else 0
        print(f"Replicate {i + 1}: {count:,} read pairs ({percentage:.1f}%)")

def detect_paired_end_files(file1, file2=None):
    """Detect if files are paired-end based on their names"""
    if file2 is None:
        return False
    
    base1 = get_base_name(file1)
    base2 = get_base_name(file2)
    
    # Remove pair indicators and check if base names are the same
    patterns = [r'_[12]$', r'_R[12]$', r'\.R[12]$']
    
    for pattern in patterns:
        clean_base1 = re.sub(pattern, '', base1)
        clean_base2 = re.sub(pattern, '', base2)
        
        if clean_base1 == clean_base2:
            return True
    
    return False

def validate_split_ratios(ratios):
    """Validate that split ratios sum to 1 and are all positive"""
    if not ratios:
        return False, "No split ratios provided"
    
    if any(ratio <= 0 for ratio in ratios):
        return False, "All split ratios must be positive"
    
    total = sum(ratios)
    if abs(total - 1.0) > 1e-6:  # Allow for small floating point errors
        return False, f"Split ratios must sum to 1.0, got {total}"
    
    return True, "Valid"

def main():
    parser = argparse.ArgumentParser(description='Split FASTQ data into multiple replicates')
    parser.add_argument('fastq1', help='First FASTQ file (or single-end file)')
    parser.add_argument('fastq2', nargs='?', help='Second FASTQ file (for paired-end)')
    parser.add_argument('--split_ratio', type=float, nargs='+', default=[0.5, 0.5],
                       help='Split ratios for each replicate (must sum to 1.0). '
                            'Example: --split_ratio 0.5 0.5 for 2 replicates, '
                            '--split_ratio 0.33 0.33 0.34 for 3 replicates')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not Path(args.fastq1).exists():
        print(f"Error: {args.fastq1} not found")
        sys.exit(1)
    
    if args.fastq2 and not Path(args.fastq2).exists():
        print(f"Error: {args.fastq2} not found")
        sys.exit(1)
    
    # Validate split ratios
    valid, message = validate_split_ratios(args.split_ratio)
    if not valid:
        print(f"Error: {message}")
        sys.exit(1)
    
    num_replicates = len(args.split_ratio)
    print(f"Splitting into {num_replicates} replicates")
    
    # Determine if single-end or paired-end
    if args.fastq2 is None:
        # Single-end
        print("Detected single-end sequencing data")
        split_single_end_fastq(args.fastq1, args.split_ratio, args.seed, batch_size=100000)
    else:
        # Paired-end
        if detect_paired_end_files(args.fastq1, args.fastq2):
            print("Detected paired-end sequencing data")
        else:
            print("Warning: File names don't follow typical paired-end naming convention")
            print("Proceeding as paired-end data...")

        split_paired_end_fastq(args.fastq1, args.fastq2, args.split_ratio, args.seed, batch_size=100000)

if __name__ == "__main__":
    main()
