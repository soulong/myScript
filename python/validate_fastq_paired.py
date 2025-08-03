#!/usr/bin/env python3
"""
Validate if two FASTQ files are properly paired-end sequencing data
Usage: python validate_paired_fastq.py file1.fq.gz file2.fq.gz
"""

import sys
import gzip
import argparse
import re
from pathlib import Path
from collections import Counter

def detect_file_format(filename):
    """Detect if file is gzipped and return appropriate open function"""
    if filename.endswith('.gz'):
        return gzip.open, 'rb'
    else:
        return open, 'r'

def parse_fastq_record(file_handle, is_gzipped=True):
    """Parse a single FASTQ record from file handle"""
    if is_gzipped:
        header = file_handle.readline().decode('utf-8').strip()
    else:
        header = file_handle.readline().strip()
    
    if not header:
        return None
    
    if is_gzipped:
        sequence = file_handle.readline().decode('utf-8').strip()
        plus = file_handle.readline().decode('utf-8').strip()
        quality = file_handle.readline().decode('utf-8').strip()
    else:
        sequence = file_handle.readline().strip()
        plus = file_handle.readline().strip()
        quality = file_handle.readline().strip()
    
    return (header, sequence, plus, quality)

def extract_read_id(header):
    """Extract read ID from FASTQ header, removing pair information"""
    # Remove @ symbol
    if header.startswith('@'):
        header = header[1:]
    
    # Common patterns for paired-end reads
    patterns = [
        r'/[12]$',      # /1, /2 at the end
        r'\.[12]$',     # .1, .2 at the end
        r'_[12]$',      # _1, _2 at the end
        r' [12]:.*$',   # space followed by 1: or 2: (Illumina format)
        r'\s+[12]:.*$', # whitespace followed by 1: or 2:
    ]
    
    clean_id = header
    for pattern in patterns:
        clean_id = re.sub(pattern, '', clean_id)
    
    # For Illumina format, keep only the part before the first space
    clean_id = clean_id.split()[0]
    
    return clean_id

def get_pair_info(header):
    """Extract pair information (1 or 2) from FASTQ header"""
    # Remove @ symbol
    if header.startswith('@'):
        header = header[1:]
    
    # Check for different pair indicators
    if re.search(r'/1$', header) or re.search(r'\.1$', header) or re.search(r'_1$', header):
        return 1
    elif re.search(r'/2$', header) or re.search(r'\.2$', header) or re.search(r'_2$', header):
        return 2
    elif re.search(r' 1:', header) or re.search(r'\s+1:', header):
        return 1
    elif re.search(r' 2:', header) or re.search(r'\s+2:', header):
        return 2
    else:
        return None

def validate_filename_pairing(file1, file2):
    """Check if filenames suggest paired-end data"""
    base1 = Path(file1).stem
    base2 = Path(file2).stem
    
    # Remove .gz extension if present
    if base1.endswith('.gz'):
        base1 = base1[:-3]
    if base2.endswith('.gz'):
        base2 = base2[:-3]
    
    # Remove common FASTQ extensions
    extensions = ['.fastq', '.fq']
    for ext in extensions:
        if base1.endswith(ext):
            base1 = base1[:-len(ext)]
        if base2.endswith(ext):
            base2 = base2[:-len(ext)]
    
    # Check for common paired-end patterns
    patterns = [
        (r'_1$', r'_2$'),
        (r'_R1$', r'_R2$'),
        (r'\.R1$', r'\.R2$'),
        (r'\.1$', r'\.2$'),
        (r'_forward$', r'_reverse$'),
        (r'_fwd$', r'_rev$'),
        (r'_f$', r'_r$'),
    ]
    
    for pattern1, pattern2 in patterns:
        if re.search(pattern1, base1) and re.search(pattern2, base2):
            clean_base1 = re.sub(pattern1, '', base1)
            clean_base2 = re.sub(pattern2, '', base2)
            if clean_base1 == clean_base2:
                return True, f"Filenames match paired-end pattern: {pattern1}/{pattern2}"
    
    return False, "Filenames don't follow common paired-end naming conventions"

def validate_read_count(file1, file2, max_reads=10000):
    """Validate that both files have the same number of reads"""
    open_func1, mode1 = detect_file_format(file1)
    open_func2, mode2 = detect_file_format(file2)
    is_gzipped1 = file1.endswith('.gz')
    is_gzipped2 = file2.endswith('.gz')
    
    count1 = 0
    count2 = 0
    
    with open_func1(file1, mode1) as f1, open_func2(file2, mode2) as f2:
        while count1 < max_reads and count2 < max_reads:
            record1 = parse_fastq_record(f1, is_gzipped1)
            record2 = parse_fastq_record(f2, is_gzipped2)
            
            if record1 is None and record2 is None:
                break
            elif record1 is None or record2 is None:
                return False, f"Files have different read counts (file1: {count1}, file2: {count2})"
            
            count1 += 1
            count2 += 1
    
    # Check if we reached the end of both files
    record1 = parse_fastq_record(f1, is_gzipped1)
    record2 = parse_fastq_record(f2, is_gzipped2)
    
    if record1 is not None or record2 is not None:
        if max_reads == count1:
            return True, f"Read counts match (checked first {max_reads} reads)"
        else:
            return False, f"Files have different read counts"
    
    return True, f"Read counts match ({count1} reads)"

def validate_read_pairing(file1, file2, max_reads=1000):
    """Validate that reads are properly paired by checking read IDs"""
    open_func1, mode1 = detect_file_format(file1)
    open_func2, mode2 = detect_file_format(file2)
    is_gzipped1 = file1.endswith('.gz')
    is_gzipped2 = file2.endswith('.gz')
    
    paired_reads = 0
    unpaired_reads = 0
    pair_info_counts = Counter()
    read_count = 0
    
    mismatched_examples = []
    
    with open_func1(file1, mode1) as f1, open_func2(file2, mode2) as f2:
        while read_count < max_reads:
            record1 = parse_fastq_record(f1, is_gzipped1)
            record2 = parse_fastq_record(f2, is_gzipped2)
            
            if record1 is None or record2 is None:
                break
            
            read_count += 1
            
            # Extract read IDs and pair info
            header1, header2 = record1[0], record2[0]
            read_id1 = extract_read_id(header1)
            read_id2 = extract_read_id(header2)
            pair_info1 = get_pair_info(header1)
            pair_info2 = get_pair_info(header2)
            
            # Count pair information
            if pair_info1:
                pair_info_counts[f"file1_pair_{pair_info1}"] += 1
            if pair_info2:
                pair_info_counts[f"file2_pair_{pair_info2}"] += 1
            
            # Check if reads are paired
            if read_id1 == read_id2:
                paired_reads += 1
                
                # Check if pair information is consistent
                if pair_info1 and pair_info2:
                    if not ((pair_info1 == 1 and pair_info2 == 2) or 
                            (pair_info1 == 2 and pair_info2 == 1)):
                        unpaired_reads += 1
                        if len(mismatched_examples) < 3:
                            mismatched_examples.append((header1, header2))
            else:
                unpaired_reads += 1
                if len(mismatched_examples) < 3:
                    mismatched_examples.append((header1, header2))
    
    pairing_rate = paired_reads / read_count if read_count > 0 else 0
    
    result = {
        'total_checked': read_count,
        'paired_reads': paired_reads,
        'unpaired_reads': unpaired_reads,
        'pairing_rate': pairing_rate,
        'pair_info_counts': dict(pair_info_counts),
        'mismatched_examples': mismatched_examples
    }
    
    return result

def validate_sequence_characteristics(file1, file2, max_reads=1000):
    """Validate sequence characteristics (length distribution, quality, etc.)"""
    open_func1, mode1 = detect_file_format(file1)
    open_func2, mode2 = detect_file_format(file2)
    is_gzipped1 = file1.endswith('.gz')
    is_gzipped2 = file2.endswith('.gz')
    
    lengths1, lengths2 = [], []
    read_count = 0
    
    with open_func1(file1, mode1) as f1, open_func2(file2, mode2) as f2:
        while read_count < max_reads:
            record1 = parse_fastq_record(f1, is_gzipped1)
            record2 = parse_fastq_record(f2, is_gzipped2)
            
            if record1 is None or record2 is None:
                break
            
            read_count += 1
            lengths1.append(len(record1[1]))  # sequence length
            lengths2.append(len(record2[1]))  # sequence length
    
    if lengths1 and lengths2:
        avg_len1 = sum(lengths1) / len(lengths1)
        avg_len2 = sum(lengths2) / len(lengths2)
        
        return {
            'reads_checked': read_count,
            'avg_length_file1': avg_len1,
            'avg_length_file2': avg_len2,
            'length_difference': abs(avg_len1 - avg_len2),
            'min_length_file1': min(lengths1),
            'max_length_file1': max(lengths1),
            'min_length_file2': min(lengths2),
            'max_length_file2': max(lengths2)
        }
    
    return None

def main():
    parser = argparse.ArgumentParser(description='Validate paired-end FASTQ files')
    parser.add_argument('file1', help='First FASTQ file')
    parser.add_argument('file2', help='Second FASTQ file')
    parser.add_argument('--max_reads', type=int, default=10000,
                       help='Maximum number of reads to check (default: 10000)')
    parser.add_argument('--quick', action='store_true',
                       help='Quick validation (check fewer reads)')
    
    args = parser.parse_args()
    
    if args.quick:
        args.max_reads = 1000
    
    # Check if files exist
    if not Path(args.file1).exists():
        print(f"Error: {args.file1} not found")
        sys.exit(1)
    
    if not Path(args.file2).exists():
        print(f"Error: {args.file2} not found")
        sys.exit(1)
    
    print(f"Validating paired-end FASTQ files:")
    print(f"File 1: {args.file1}")
    print(f"File 2: {args.file2}")
    print(f"Max reads to check: {args.max_reads}")
    print("-" * 60)
    
    # 1. Filename validation
    print("1. Filename validation:")
    filename_valid, filename_msg = validate_filename_pairing(args.file1, args.file2)
    print(f"   {filename_msg}")
    print(f"   Result: {'PASS' if filename_valid else 'WARNING'}")
    print()
    
    # 2. Read count validation
    print("2. Read count validation:")
    try:
        count_valid, count_msg = validate_read_count(args.file1, args.file2, args.max_reads)
        print(f"   {count_msg}")
        print(f"   Result: {'PASS' if count_valid else 'FAIL'}")
    except Exception as e:
        print(f"   Error: {e}")
        print(f"   Result: FAIL")
        count_valid = False
    print()
    
    # 3. Read pairing validation
    print("3. Read pairing validation:")
    try:
        pairing_result = validate_read_pairing(args.file1, args.file2, min(args.max_reads, 1000))
        
        print(f"   Reads checked: {pairing_result['total_checked']}")
        print(f"   Properly paired: {pairing_result['paired_reads']}")
        print(f"   Unpaired: {pairing_result['unpaired_reads']}")
        print(f"   Pairing rate: {pairing_result['pairing_rate']:.1%}")
        
        if pairing_result['pair_info_counts']:
            print(f"   Pair indicators found: {pairing_result['pair_info_counts']}")
        
        pairing_valid = pairing_result['pairing_rate'] > 0.95
        print(f"   Result: {'PASS' if pairing_valid else 'FAIL'}")
        
        if pairing_result['mismatched_examples'] and not pairing_valid:
            print("   Examples of mismatched reads:")
            for i, (h1, h2) in enumerate(pairing_result['mismatched_examples'][:3]):
                print(f"     {i+1}. {h1}")
                print(f"        {h2}")
    
    except Exception as e:
        print(f"   Error: {e}")
        print(f"   Result: FAIL")
        pairing_valid = False
    print()
    
    # 4. Sequence characteristics validation
    print("4. Sequence characteristics:")
    try:
        seq_result = validate_sequence_characteristics(args.file1, args.file2, min(args.max_reads, 1000))
        
        if seq_result:
            print(f"   Reads checked: {seq_result['reads_checked']}")
            print(f"   Average length file1: {seq_result['avg_length_file1']:.1f}")
            print(f"   Average length file2: {seq_result['avg_length_file2']:.1f}")
            print(f"   Length difference: {seq_result['length_difference']:.1f}")
            print(f"   Length range file1: {seq_result['min_length_file1']}-{seq_result['max_length_file1']}")
            print(f"   Length range file2: {seq_result['min_length_file2']}-{seq_result['max_length_file2']}")
            
            length_valid = seq_result['length_difference'] < 10  # Allow some variation
            print(f"   Result: {'PASS' if length_valid else 'WARNING'}")
        else:
            print("   Could not analyze sequence characteristics")
            print("   Result: FAIL")
            length_valid = False
    
    except Exception as e:
        print(f"   Error: {e}")
        print(f"   Result: FAIL")
        length_valid = False
    print()
    
    # Overall result
    print("=" * 60)
    print("OVERALL VALIDATION RESULT:")
    
    if count_valid and pairing_valid:
        print("✓ Files appear to be valid paired-end FASTQ data")
        sys.exit(0)
    elif count_valid and pairing_result['pairing_rate'] > 0.8:
        print("⚠ Files are likely paired-end but may have some issues")
        sys.exit(1)
    else:
        print("✗ Files do not appear to be properly paired-end FASTQ data")
        sys.exit(2)

if __name__ == "__main__":
    main()