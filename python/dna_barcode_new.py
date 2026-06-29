#!/usr/bin/env python3
# MIT License
# Copyright (c) 2020 David Jacob Feldman
# Modified for performance optimization

import argparse
from collections import defaultdict
import logging
import os
import re
import sys
import time
import random
import gc as garbage_collector
import traceback

logger = logging.getLogger(__name__)

# ============================================================================
# Fast Levenshtein distance with early termination
# ============================================================================
def levenshtein_distance_cutoff(s1, s2, cutoff):
    """Calculate Levenshtein distance with early termination."""
    len1, len2 = len(s1), len(s2)
    
    if abs(len1 - len2) > cutoff:
        return cutoff + 1
    
    # Use arrays for better memory efficiency
    prev_row = list(range(len2 + 1))
    curr_row = [0] * (len2 + 1)
    
    for i in range(1, len1 + 1):
        curr_row[0] = i
        min_val = i
        
        for j in range(1, len2 + 1):
            if s1[i-1] == s2[j-1]:
                curr_row[j] = prev_row[j-1]
            else:
                curr_row[j] = 1 + min(prev_row[j-1], prev_row[j], curr_row[j-1])
            
            min_val = min(min_val, curr_row[j])
        
        if min_val > cutoff:
            return cutoff + 1
        
        # Swap rows
        prev_row, curr_row = curr_row, prev_row
    
    return prev_row[len2]

# ============================================================================
# Sequence generation with filtering
# ============================================================================
INT_TO_NUC = ['A', 'C', 'G', 'T']

def count_gc_int(num, length):
    """Count GC content from integer representation."""
    gc_count = 0
    for _ in range(length):
        bits = num & 3
        if bits == 1 or bits == 2:
            gc_count += 1
        num >>= 2
    return gc_count

def has_homopolymer_int(num, length, homopolymer):
    """Check for homopolymers."""
    if length < homopolymer:
        return False
    
    consecutive = 1
    prev_bits = num & 3
    num >>= 2
    
    for _ in range(1, length):
        curr_bits = num & 3
        if curr_bits == prev_bits:
            consecutive += 1
            if consecutive >= homopolymer:
                return True
        else:
            consecutive = 1
        prev_bits = curr_bits
        num >>= 2
    
    return False

def int_to_seq(num, length):
    """Convert integer to DNA sequence."""
    seq = []
    for _ in range(length):
        seq.append(INT_TO_NUC[num & 3])
        num >>= 2
    return ''.join(reversed(seq))

def generate_sequences(length, gc_min, gc_max, homopolymer, target_count, exclude_pattern=None):
    """Generate sequences with all filters applied."""
    random.seed(42)
    total = 4 ** length
    
    gc_min_count = int(gc_min * length)
    gc_max_count = int(gc_max * length)
    
    sequences = []
    seen = set()
    attempts = 0
    max_attempts = target_count * 200
    
    pattern = re.compile(exclude_pattern) if exclude_pattern else None
    
    logger.info(f"Generating {target_count:,} sequences (length={length}, GC={gc_min}-{gc_max}, hp<{homopolymer})...")
    
    while len(sequences) < target_count and attempts < max_attempts:
        num = random.randint(0, total - 1)
        
        # Quick GC check
        gc_count = count_gc_int(num, length)
        if gc_count < gc_min_count or gc_count > gc_max_count:
            attempts += 1
            continue
        
        # Homopolymer check
        if has_homopolymer_int(num, length, homopolymer):
            attempts += 1
            continue
        
        # Convert to sequence
        seq = int_to_seq(num, length)
        
        # Exclude pattern check
        if pattern and pattern.search(seq):
            attempts += 1
            continue
        
        # Duplicate check
        if seq in seen:
            attempts += 1
            continue
        
        seen.add(seq)
        sequences.append(seq)
        
        if (attempts + 1) % 50000 == 0:
            logger.info(f"  Progress: {len(sequences):,}/{target_count:,} sequences found ({100*len(sequences)/max_attempts:.1f}% of max attempts)")
    
    logger.info(f"Generated {len(sequences):,} sequences from {attempts:,} attempts")
    return sequences

# ============================================================================
# Chunked greedy selection - memory efficient
# ============================================================================
def greedy_select_chunked(sequences, k, chunk_size=10000):
    """Greedy selection in chunks to manage memory better."""
    if not sequences:
        return []
    
    random.seed(42)
    candidates = sequences.copy()
    random.shuffle(candidates)
    
    selected = []
    total = len(candidates)
    
    for idx, candidate in enumerate(candidates):
        # Check distance to all selected sequences
        valid = True
        for sel in selected:
            d = levenshtein_distance_cutoff(candidate, sel, k - 1)
            if d < k:
                valid = False
                break
        
        if valid:
            selected.append(candidate)
        
        # Progress update every 2000 sequences
        if (idx + 1) % 2000 == 0:
            pct = 100 * (idx + 1) / total
            logger.info(f"  Progress: {pct:.1f}% - selected {len(selected):,} barcodes so far")
            
            # Garbage collection every 10% to free memory
            if (idx + 1) % (total // 10) == 0:
                garbage_collector.collect()
                logger.info(f"  [GC] Memory cleanup performed")
    
    return selected

# ============================================================================
# Main function
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description='DNA Barcode Generator (Optimized)')
    parser.add_argument('--distance', '-d', type=int, default=4, help='Minimum edit distance')
    parser.add_argument('--length', '-l', type=int, default=16, help='Barcode length')
    parser.add_argument('--gc_min', type=int, default=40, help='Min GC %')
    parser.add_argument('--gc_max', type=int, default=60, help='Max GC %')
    parser.add_argument('--homopolymer', type=int, default=3, help='Max homopolymer length')
    parser.add_argument('--limit', type=int, default=100000, help='Max input sequences')
    parser.add_argument('--exclude', type=str, default='ATG|TAA|TGA|TAG', help='Exclude pattern')
    parser.add_argument('--verbosity', '-v', type=int, default=2, help='Verbosity level')
    parser.add_argument('--cores', type=int, default=1, help='Cores (not used)')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        format='%(asctime)s -- %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=args.verbosity * 10
    )
    
    logger.info("=" * 60)
    logger.info("DNA Barcode Generator - Optimized Version")
    logger.info("=" * 60)
    
    total_start = time.time()
    
    try:
        # Generate sequences
        sequences = generate_sequences(
            length=args.length,
            gc_min=args.gc_min / 100,
            gc_max=args.gc_max / 100,
            homopolymer=args.homopolymer,
            target_count=args.limit,
            exclude_pattern=args.exclude
        )
        
        if not sequences:
            logger.error("No sequences generated!")
            return
        
        # Clear memory
        garbage_collector.collect()
        
        # Greedy selection
        logger.info(f"Starting greedy selection (distance >= {args.distance})...")
        select_start = time.time()
        barcodes = greedy_select_chunked(sequences, args.distance)
        select_time = time.time() - select_start
        
        logger.info(f"Selection completed in {select_time:.1f}s")
        logger.info(f"Selected {len(barcodes):,} barcodes from {len(sequences):,} candidates")
        logger.info(f"Selection rate: {100*len(barcodes)/len(sequences):.1f}%")
        
        # Validation (sample for large sets)
        validate_count = min(len(barcodes), 2000)
        logger.info(f"Validating {validate_count:,} barcodes...")
        validate_sample = barcodes[:validate_count]
        failures = 0
        checked = 0
        
        for i in range(len(validate_sample)):
            for j in range(i + 1, len(validate_sample)):
                d = levenshtein_distance_cutoff(validate_sample[i], validate_sample[j], args.distance - 1)
                checked += 1
                if d < args.distance:
                    failures += 1
                    if failures <= 5:
                        logger.warning(f"  Failure: {validate_sample[i]} vs {validate_sample[j]} (d={d})")
        
        if failures > 0:
            logger.error(f"Validation FAILED: {failures} errors in {checked:,} checks")
        else:
            logger.info(f"Validation PASSED: {checked:,} pairs checked")
        
        # Save results
        import csv
        filename = f"barcodes_n{args.length}_k{args.distance}_{len(barcodes)}.csv"
        
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['barcode', 'gc_content'])
            for bc in barcodes:
                gc = (bc.count('G') + bc.count('C')) / len(bc)
                writer.writerow([bc, f"{gc:.3f}"])
        
        total_time = time.time() - total_start
        logger.info(f"Saved to {filename}")
        logger.info(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"Error occurred: {e}")
        logger.error(traceback.format_exc())
        raise

if __name__ == '__main__':
    main()
