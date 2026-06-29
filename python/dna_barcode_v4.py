#!/usr/bin/env python3
# MIT License
# Copyright (c) 2020 David Jacob Feldman
# Modified for performance optimization and target-oriented design
# Version 3: Adaptive strategy + target count semantics

"""
DNA Barcode Generator v3

Key features:
- Target-oriented: limit parameter specifies desired OUTPUT count (not input sampling)
- Adaptive strategy: automatically selects optimal algorithm based on barcode length
- Multi-round iteration: keeps sampling until target count is reached
- Optimized filtering: k-mer hashing for short sequences, direct comparison for long ones

Usage:
    python dna_barcode_v3.py --length 12 --distance 3 --limit 1000
"""

import argparse
from collections import defaultdict, Counter
import logging
import re
import time
import random
import csv
import numpy as np
import scipy.sparse
from typing import List, Dict, Tuple, Optional

logger = logging.getLogger(__name__)

# ============================================================================
# Constants and Configuration
# ============================================================================

INT_TO_NUC = ['A', 'C', 'G', 'T']
NUC_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# Strategy thresholds
EXHAUSTIVE_MAX_LENGTH = 10      # Use exhaustive generation for n <= 10
KHASH_MIN_SELECTED = 100        # Use k-mer hashing when selected count >= 100
DEFAULT_BATCH_SIZE = 10000      # Default sequences per sampling round
MAX_ROUNDS = 100                # Maximum sampling rounds before giving up

# ============================================================================
# Strategy Selection
# ============================================================================

def get_strategy_for_length(n: int) -> Dict:
    """
    Select optimal strategy based on barcode length.
    
    Short sequences (n<=10): exhaustive enumeration is feasible and gives optimal results
    Medium sequences (11-16): sampling with k-mer hashing for efficiency
    Long sequences (n>16): sampling with direct comparison (k-mer overhead not worth it)
    """
    if n <= EXHAUSTIVE_MAX_LENGTH:
        return {
            'generation': 'exhaustive',
            'selection': 'maxy_clique',
            'use_khash': True,
            'batch_size': None  # Not used for exhaustive
        }
    elif n <= 16:
        return {
            'generation': 'sampling',
            'selection': 'greedy_khash',
            'use_khash': True,
            'batch_size': DEFAULT_BATCH_SIZE
        }
    else:
        return {
            'generation': 'sampling',
            'selection': 'greedy_direct',
            'use_khash': False,
            'batch_size': 100000
        }


def validate_params(n: int, k: int, gc_min: float, gc_max: float) -> List[str]:
    """
    Validate input parameters and return list of warnings.
    Raises ValueError for parameters that make barcode generation impossible.
    """
    warnings = []

    if n < 4:
        warnings.append(f"Warning: length={n} is very short, may not find enough barcodes")

    if k < 1:
        raise ValueError(f"distance={k} must be at least 1")

    if k > n:
        raise ValueError(f"distance={k} > length={n}: impossible to satisfy, no valid barcodes exist")

    if gc_min < 0 or gc_min > 1:
        raise ValueError(f"gc_min={gc_min} must be in range [0, 1]")

    if gc_max < 0 or gc_max > 1:
        raise ValueError(f"gc_max={gc_max} must be in range [0, 1]")

    if gc_min > gc_max:
        raise ValueError(f"gc_min={gc_min} > gc_max={gc_max}: invalid GC range")

    return warnings

# ============================================================================
# Fast Levenshtein Distance with Early Termination
# ============================================================================

def levenshtein_distance_cutoff(s1: str, s2: str, cutoff: int) -> int:
    """
    Calculate Levenshtein distance with early termination.
    Returns immediately if distance exceeds cutoff, saving computation time.
    """
    len1, len2 = len(s1), len(s2)
    
    # Quick length check
    if abs(len1 - len2) > cutoff:
        return cutoff + 1
    
    # Use two rows for space efficiency
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
        
        # Early termination: if minimum value in row exceeds cutoff, stop
        if min_val > cutoff:
            return cutoff + 1
        
        prev_row, curr_row = curr_row, prev_row
    
    return prev_row[len2]

# ============================================================================
# Sequence Generation - Integer Encoding
# ============================================================================

def count_gc_int(num: int, length: int) -> int:
    """
    Count GC content from integer representation.
    Each nucleotide is encoded as 2 bits: A=00, C=01, G=10, T=11
    C and G have bit pattern *1, so we count bits where (bits & 1) == 1
    """
    gc_count = 0
    for _ in range(length):
        bits = num & 3
        if bits == 1 or bits == 2:  # C=01, G=10
            gc_count += 1
        num >>= 2
    return gc_count


def has_homopolymer_int(num: int, length: int, homopolymer: int) -> bool:
    """
    Check for homopolymers (runs of identical nucleotides of length >= homopolymer).
    Converts integer to sequence first for clarity, then does a single-pass run scan.
    """
    if length < homopolymer:
        return False

    seq = int_to_seq(num, length)

    run = 1
    for i in range(1, length):
        if seq[i] == seq[i - 1]:
            run += 1
            if run >= homopolymer:
                return True
        else:
            run = 1

    return False


def int_to_seq(num: int, length: int) -> str:
    """
    Convert integer to DNA sequence.
    Least significant 2 bits represent the first nucleotide (matching seq_to_int encoding).
    """
    seq = []
    for _ in range(length):
        seq.append(INT_TO_NUC[num & 3])
        num >>= 2
    return ''.join(seq)


def seq_to_int(seq: str) -> int:
    """
    Convert DNA sequence to integer representation.
    """
    num = 0
    for nuc in reversed(seq):
        num = (num << 2) | NUC_TO_INT[nuc]
    return num


def generate_all_barcodes_exhaustive(n: int, gc_min: float, gc_max: float, 
                                      homopolymer: int, exclude_pattern: Optional[re.Pattern] = None) -> List[str]:
    """
    Generate ALL possible barcodes of length n (exhaustive enumeration).
    Only feasible for n <= 10 (4^10 = 1,048,576 sequences).
    
    Applies all filters during generation to reduce memory usage.
    """
    total = 4 ** n
    gc_min_count = round(gc_min * n)
    gc_max_count = round(gc_max * n)
    barcodes = []

    for num in range(total):
        # Quick GC check
        gc_count = count_gc_int(num, n)
        if gc_count < gc_min_count or gc_count > gc_max_count:
            continue
        
        # Homopolymer check
        if has_homopolymer_int(num, n, homopolymer):
            continue
        
        # Convert to sequence
        seq = int_to_seq(num, n)
        
        # Exclude pattern check
        if exclude_pattern and exclude_pattern.search(seq):
            continue
        
        barcodes.append(seq)
    
    logger.info(f"Exhaustive generation complete: {len(barcodes):,} sequences passed filters")
    return barcodes


def generate_sequences_batch(n: int, gc_min: float, gc_max: float, homopolymer: int,
                             target_count: int, exclude_pattern: Optional[re.Pattern] = None,
                             seed: Optional[int] = None) -> List[str]:
    """
    Generate a batch of sequences via random sampling.
    Used for longer barcodes where exhaustive enumeration is infeasible.
    """
    if seed is not None:
        random.seed(seed)
    
    total = 4 ** n
    gc_min_count = round(gc_min * n)
    gc_max_count = round(gc_max * n)
    seen = set()
    sequences = []
    attempts = 0
    max_attempts = target_count * 200  # Safety limit
    
    logger.info(f"Sampling {target_count:,} sequences (length={n}, GC={gc_min}-{gc_max}, hp<{homopolymer})...")
    
    while len(sequences) < target_count and attempts < max_attempts:
        attempts += 1
        num = random.randint(0, total - 1)

        # Quick GC check
        gc_count = count_gc_int(num, n)
        if gc_count < gc_min_count or gc_count > gc_max_count:
            continue

        # Homopolymer check
        if has_homopolymer_int(num, n, homopolymer):
            continue

        # Convert to sequence
        seq = int_to_seq(num, n)

        # Exclude pattern check
        if exclude_pattern and exclude_pattern.search(seq):
            continue

        # Duplicate check
        if seq in seen:
            continue

        seen.add(seq)
        sequences.append(seq)

        if attempts % 50000 == 0:
            logger.info(f"  Progress: {len(sequences):,}/{target_count:,} found after {attempts:,} attempts")
    
    logger.info(f"Sampling complete: {len(sequences):,} sequences from {attempts:,} attempts")
    return sequences

# ============================================================================
# K-mer Hashing for Efficient Filtering
# ============================================================================

def khash(s: str, k: int) -> List[Tuple[int, str]]:
    """
    Divide sequence into substrings for k-mer hashing.

    Mathematical property: Two strings of the same length with Levenshtein
    distance < k will share at least one substring of length ceil((n-k)/k).

    Note: Uses LINEAR (non-circular) hashing for correct Levenshtein distance filtering.
    """
    n = len(s)
    window = int(np.ceil((n - k) / float(k)))
    seen = set()
    arr = []

    # Linear hashing: only consider substrings within the sequence
    for i in range(n):
        for j in (0, 1):
            pos = i + j
            if pos < n:
                kmer = s[pos:pos + window]
                if len(kmer) == window and (pos, kmer) not in seen:
                    seen.add((pos, kmer))
                    arr.append((pos, kmer))

    return arr


def build_kmer_index(sequences: List[str], k: int) -> Dict:
    """
    Build k-mer index for a set of sequences.
    Returns a dictionary mapping k-mer hashes to sequence indices.
    """
    index = defaultdict(set)
    
    for idx, seq in enumerate(sequences):
        for pos, kmer in khash(seq, k):
            index[(pos, kmer)].add(idx)
    
    return index


def conflicts_with_index(candidate: str, kmer_index: Dict, sequences: List[str], k: int) -> bool:
    """
    Check if a candidate sequence conflicts with any indexed sequence.
    Uses k-mer hashing to find candidates, then verifies with actual distance.
    """
    neighbor_indices = set()
    for pos, kmer in khash(candidate, k):
        if (pos, kmer) in kmer_index:
            neighbor_indices.update(kmer_index[(pos, kmer)])

    for idx in neighbor_indices:
        d = levenshtein_distance_cutoff(candidate, sequences[idx], k - 1)
        if d < k:
            return True

    return False

# ============================================================================
# Filtering and Selection
# ============================================================================

def filter_direct(candidates: List[str], selected: List[str], k: int) -> List[str]:
    """
    Filter candidates by checking distance against all selected sequences.
    Simple O(n*m) approach, efficient for small selected sets.
    """
    result = []
    
    for candidate in candidates:
        valid = True
        for sel in selected:
            d = levenshtein_distance_cutoff(candidate, sel, k - 1)
            if d < k:
                valid = False
                break
        
        if valid:
            result.append(candidate)
    
    return result


def filter_with_khash(candidates: List[str], selected: List[str], k: int) -> List[str]:
    """
    Filter candidates using k-mer hashing optimization.
    More efficient when selected set is large.
    """
    # Build k-mer index for selected sequences
    kmer_index = build_kmer_index(selected, k)

    result = []
    for candidate in candidates:
        if not conflicts_with_index(candidate, kmer_index, selected, k):
            result.append(candidate)

    return result


def filter_compatible(candidates: List[str], selected: List[str], k: int, 
                      use_khash: bool = True) -> List[str]:
    """
    Filter candidates to keep only those compatible with selected sequences.
    Automatically chooses between direct comparison and k-mer hashing.
    """
    if not selected:
        return candidates
    
    # Use k-mer hashing only when beneficial
    if use_khash and len(selected) >= KHASH_MIN_SELECTED:
        logger.debug(f"Using k-mer hashing for filtering (selected={len(selected)})")
        return filter_with_khash(candidates, selected, k)
    else:
        logger.debug(f"Using direct comparison for filtering (selected={len(selected)})")
        return filter_direct(candidates, selected, k)

# ============================================================================
# Maximum Clique Selection (for exhaustive mode)
# ============================================================================

def build_khash_dict(xs: List[str], k: int) -> Dict:
    """Build k-mer hash buckets for distance calculation."""
    D = defaultdict(list)
    for x in xs:
        for h in khash(x, k):
            D[h].append(x)
    
    # Use 'key' instead of 'k' to avoid shadowing the parameter
    return {key: sorted(set(v)) for key, v in D.items()}


def sparse_dist(hash_buckets: List[List[str]], k: int, 
                progress: callable = None) -> Dict:
    """
    Calculate pairwise distances within hash buckets.
    Only records distances < k (conflicts).
    """
    if progress is None:
        progress = lambda x: x
    
    D = {}
    for xs in progress(hash_buckets):
        for i, a in enumerate(xs):
            for b in xs[i+1:]:
                d = levenshtein_distance_cutoff(a, b, k - 1)
                if d < k:
                    key = tuple(sorted((a, b)))
                    D[key] = d
    
    return D


def maxy_clique_groups(cm, group_ids: List[int], verbose: bool = False) -> List[int]:
    """
    Select maximum clique using greedy algorithm.
    Prioritizes groups with fewest remaining candidates.
    """
    
    # counts => group_id
    d1 = defaultdict(set)
    for id_, counts in Counter(group_ids).items():
        d1[counts].add(id_)
    
    # group_id => indices
    d2 = defaultdict(list)
    for i, id_ in enumerate(group_ids):
        d2[id_].append(i)
    d2 = {k: v[::-1] for k, v in d2.items()}  # Reverse for efficient pop()
    
    # group_id => # selected
    d3 = Counter()
    
    selected = []
    available = np.array(range(len(group_ids)))
    available_set = set(range(len(group_ids)))

    while d1:
        if verbose and (len(selected) % 1000) == 0:
            logger.info(f"Selected {len(selected)} barcodes...")

        # Pick group with fewest remaining
        count = min(d1.keys())
        id_ = d1[count].pop()

        if len(d1[count]) == 0:
            d1.pop(count)

        # Find available index from this group
        index = None
        while d2[id_]:
            index = d2[id_].pop()
            if index in available_set:
                break
        else:
            index = None

        # Add to selection if compatible
        if index is not None:
            selected.append(index)
            d3[id_] += 1
            available_set.discard(index)
            available = available[available != index]

            # Remove incompatible barcodes
            remove = cm[index, available].indices
            mask = np.ones(len(available), dtype=bool)
            mask[remove] = False
            removed_indices = available[~mask]
            for ri in removed_indices:
                available_set.discard(int(ri))
            available = available[mask]

        # Move group to new bin
        n = len(d2[id_])
        if n > 0:
            d1[n].add(id_)
    
    return selected

# ============================================================================
# Target-Oriented Main Algorithm
# ============================================================================

def create_barcode_set_target(n: int, k: int, gc_min: float, gc_max: float,
                               homopolymer: int, target_count: int,
                               exclude: Optional[str] = None,
                               max_rounds: int = MAX_ROUNDS,
                               cores: int = 1) -> List[str]:
    """
    Generate barcodes with target output count.
    
    This is the main entry point. Uses adaptive strategy based on barcode length.
    
    Args:
        n: Barcode length
        k: Minimum edit distance
        gc_min: Minimum GC content (0-1 scale)
        gc_max: Maximum GC content (0-1 scale)
        homopolymer: Maximum allowed homopolymer length
        target_count: Desired number of output barcodes (KEY PARAMETER)
        exclude: Regex pattern to exclude
        max_rounds: Maximum sampling rounds
        cores: Number of parallel processes (not yet implemented)
    
    Returns:
        List of barcodes meeting all criteria (up to target_count)
    """
    # Get strategy for this length
    strategy = get_strategy_for_length(n)
    logger.info(f"Using strategy for length {n}: {strategy}")
    
    # Compile exclude pattern
    exclude_pattern = re.compile(exclude) if exclude else None
    
    selected = []
    
    if strategy['generation'] == 'exhaustive':
        # === Exhaustive mode (n <= 10) ===
        logger.info("=== EXHAUSTIVE MODE ===")
        
        # Generate all possible barcodes
        all_barcodes = generate_all_barcodes_exhaustive(
            n, gc_min, gc_max, homopolymer, exclude_pattern
        )
        
        if len(all_barcodes) <= target_count * 2:
            # Small set: use maximum clique algorithm
            logger.info(f"Using maximum clique selection on {len(all_barcodes):,} barcodes...")
            
            # Build distance matrix
            hash_buckets_dict = build_khash_dict(all_barcodes, k)
            hash_buckets = list(hash_buckets_dict.values())
            D = sparse_dist(hash_buckets, k)
            
            # Create sparse conflict matrix
            cm = sparse_view(all_barcodes, D)
            
            # Select maximum clique
            group_ids = [0] * len(all_barcodes)
            selected_indices = maxy_clique_groups(cm, group_ids)
            selected = [all_barcodes[i] for i in selected_indices]
            
        else:
            # Large set: use greedy selection
            logger.info(f"Using greedy selection on {len(all_barcodes):,} barcodes...")
            selected = greedy_select_simple(all_barcodes, k)
        
        logger.info(f"Exhaustive mode complete: {len(selected):,} barcodes selected")
    
    else:
        # === Sampling mode (n > 10) ===
        logger.info("=== SAMPLING MODE ===")
        
        base_batch_size = strategy.get('batch_size', DEFAULT_BATCH_SIZE)
        use_khash = strategy.get('use_khash', True)

        for round_idx in range(max_rounds):
            if len(selected) >= target_count:
                break

            # Calculate remaining needed
            remaining = target_count - len(selected)

            # Scale batch to remaining need, bounded by the strategy's base size
            batch_size = min(base_batch_size, remaining * 5)
            
            # Sample new candidates
            candidates = generate_sequences_batch(
                n, gc_min, gc_max, homopolymer,
                batch_size, exclude_pattern,
                seed=round_idx
            )
            
            if not candidates:
                logger.warning(f"Round {round_idx + 1}: No candidates generated")
                break
            
            # Step 1: Filter candidates internally (ensure new candidates don't conflict with each other)
            # This is critical: without this, a batch of 50k sequences has no internal filtering
            # Limit output to 1.5x remaining to avoid wasting time on sequences we won't use
            if len(candidates) > 1:
                max_internal = min(int(remaining * 1.5), 2000)  # Cap at 2000 for performance
                logger.info(f"  Internal filtering: {len(candidates):,} candidates -> max {max_internal:,}")
                candidates = greedy_filter_batch(candidates, k, max_output=max_internal, seed=round_idx)
                logger.info(f"  After internal filter: {len(candidates):,} remain")
            
            # Step 2: Filter compatible with already selected
            new_selected = filter_compatible(candidates, selected, k, use_khash)
            
            selected.extend(new_selected)
            
            logger.info(f"Round {round_idx + 1}/{max_rounds}: "
                       f"added {len(new_selected):,}, "
                       f"total {len(selected):,}/{target_count:,}")
        
        logger.info(f"Sampling mode complete: {len(selected):,} barcodes selected")
    
    # Post-process: trim to target count if exceeded
    if len(selected) > target_count:
        random.seed(42)
        selected = random.sample(selected, target_count)
        logger.info(f"Trimmed to target count: {target_count:,}")
    
    # Warn if target not met
    if len(selected) < target_count:
        logger.warning(f"Warning: Only generated {len(selected):,}/{target_count:,} barcodes")
        logger.warning("Possible causes: parameters too strict (narrow GC range, "
                      "high distance, short length)")
    
    return selected

# ============================================================================
# Greedy Selection (Simplified)
# ============================================================================

def greedy_select_simple(sequences: List[str], k: int) -> List[str]:
    """
    Simple greedy selection: shuffle and add if compatible.
    Used as fallback for large exhaustive sets.
    """
    if not sequences:
        return []
    
    random.seed(42)
    candidates = sequences.copy()
    random.shuffle(candidates)
    
    selected = []
    total = len(candidates)
    
    for idx, candidate in enumerate(candidates):
        valid = True
        for sel in selected:
            d = levenshtein_distance_cutoff(candidate, sel, k - 1)
            if d < k:
                valid = False
                break
        
        if valid:
            selected.append(candidate)
        
        if (idx + 1) % 2000 == 0:
            logger.info(f"  Progress: {100*(idx+1)/total:.1f}% - "
                       f"selected {len(selected):,} barcodes")
    
    return selected


def greedy_filter_batch(candidates: List[str], k: int, max_output: int = None,
                        seed: Optional[int] = None) -> List[str]:
    """
    Filter a batch of candidates to ensure mutual compatibility.
    Uses greedy selection: shuffle and add if compatible with already selected.
    
    This is critical for sampling mode: without internal filtering, a batch of
    sequences has no distance guarantees within the batch itself.
    
    Args:
        candidates: List of candidate sequences
        k: Minimum distance threshold
        max_output: Stop when this many sequences are selected (optional)
    """
    if len(candidates) <= 1:
        return candidates
    
    # Shuffle for randomness
    if seed is not None:
        random.seed(seed)
    shuffled = candidates.copy()
    random.shuffle(shuffled)
    
    selected = []
    total = len(shuffled)
    
    for idx, candidate in enumerate(shuffled):
        valid = True
        for sel in selected:
            d = levenshtein_distance_cutoff(candidate, sel, k - 1)
            if d < k:
                valid = False
                break
        
        if valid:
            selected.append(candidate)
            # Early exit if we have enough
            if max_output and len(selected) >= max_output:
                logger.info(f"  Greedy filter: reached max_output {max_output:,}")
                break
        
        # Progress logging for large batches
        if total > 10000 and (idx + 1) % 10000 == 0:
            logger.info(f"  Greedy filter progress: {100*(idx+1)/total:.1f}% - "
                        f"selected {len(selected):,}/{total:,}")
    
    return selected

# ============================================================================
# Sparse Matrix View
# ============================================================================

def sparse_view(xs: List[str], D: Dict, symmetric: bool = True):
    """
    Create sparse conflict matrix from distance dictionary.
    """
    if len(xs) != len(set(xs)):
        raise ValueError("sparse_view requires unique sequences in xs")

    mapper = {x: i for i, x in enumerate(xs)}

    if len(D) == 0:
        i, j, data = [], [], []
    else:
        i, j, data = zip(*[(mapper[a], mapper[b], True) for (a, b), v in D.items()])
        i = np.array(i)
        j = np.array(j)
        data = np.array(data)

    n = len(xs)
    cm = scipy.sparse.coo_matrix((data, (i, j)), shape=(n, n))

    if symmetric:
        cm = (cm + cm.T).tocsr()

    return cm

# ============================================================================
# Validation
# ============================================================================

def validate_barcode_set(barcodes: List[str], k: int,
                         max_to_check: int = int(1e6)) -> List[Tuple[str, str, int]]:
    """
    Validate that all barcode pairs have distance >= k.
    Returns list of failures (pairs with distance < k).
    """
    from itertools import combinations

    failures = []
    checked = 0

    logger.info(f"Validating {len(barcodes):,} barcodes...")

    for b1, b2 in combinations(barcodes, 2):
        d = levenshtein_distance_cutoff(b1, b2, k - 1)
        checked += 1

        if d < k:
            failures.append((b1, b2, d))

        if checked >= max_to_check:
            logger.info(f"Checked {max_to_check:,} pairs (limit reached)")
            return failures

    logger.info(f"Validation complete: {checked:,} pairs checked")
    return failures

# ============================================================================
# Argument Parsing
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='DNA Barcode Generator v3 - Target-oriented design',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--distance', '-d', type=int, default=3,
                        help='Minimum edit distance (k)')
    parser.add_argument('--length', '-l', type=int, default=12,
                        help='Barcode length (n)')
    parser.add_argument('--gc_min', type=int, default=40,
                        help='Minimum GC percentage (0-100)')
    parser.add_argument('--gc_max', type=int, default=60,
                        help='Maximum GC percentage (0-100)')
    parser.add_argument('--homopolymer', type=int, default=3,
                        help='Maximum homopolymer length')
    parser.add_argument('--limit', type=int, default=1000,
                        help='Target number of output barcodes (KEY PARAMETER)')
    parser.add_argument('--exclude', type=str, default='ATG|TAA|TGA|TAG',
                        help='Regex pattern to exclude in barcodes')
    parser.add_argument('--max_rounds', type=int, default=100,
                        help='Maximum sampling rounds (for n > 10)')
    parser.add_argument('--verbosity', '-v', type=int, default=2,
                        help='Logging verbosity level')
    parser.add_argument('--cores', type=int, default=1,
                        help='Number of parallel processes (not yet implemented)')
    parser.add_argument('--validate', action='store_true',
                        help='Run full validation after generation')
    
    return parser.parse_args()

# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    args = parse_args()
    
    # Setup logging
    verbosity_map = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    log_level = verbosity_map.get(args.verbosity, logging.DEBUG)
    logging.basicConfig(
        format='%(asctime)s -- %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=log_level
    )
    
    logger.info("=" * 70)
    logger.info("DNA Barcode Generator v3")
    logger.info("=" * 70)
    logger.info(f"Parameters:")
    logger.info(f"  Length: {args.length}")
    logger.info(f"  Distance: {args.distance}")
    logger.info(f"  GC range: {args.gc_min}% - {args.gc_max}%")
    logger.info(f"  Max homopolymer: {args.homopolymer}")
    logger.info(f"  Target count: {args.limit:,}")
    logger.info(f"  Exclude pattern: {args.exclude}")
    logger.info("=" * 70)
    
    # Validate parameters
    warnings = validate_params(
        args.length, args.distance,
        args.gc_min / 100, args.gc_max / 100
    )
    for warning in warnings:
        logger.warning(warning)

    if args.cores > 1:
        logger.warning(f"--cores={args.cores} requested but parallelism is not yet implemented; running single-threaded")
    
    total_start = time.time()
    
    # Generate barcodes
    barcodes = create_barcode_set_target(
        n=args.length,
        k=args.distance,
        gc_min=args.gc_min / 100,
        gc_max=args.gc_max / 100,
        homopolymer=args.homopolymer,
        target_count=args.limit,
        exclude=args.exclude if args.exclude else None,
        max_rounds=args.max_rounds,
        cores=args.cores
    )
    
    generation_time = time.time() - total_start
    logger.info(f"Generation completed in {generation_time:.1f}s")
    
    # Validation
    if args.validate or len(barcodes) <= 1000:
        logger.info(f"Running full validation...")
        failures = validate_barcode_set(barcodes, args.distance)
        
        if failures:
            logger.error(f"Validation FAILED: {len(failures)} errors found")
            for failure in failures[:5]:
                logger.error(f"  {failure[0]} vs {failure[1]} (d={failure[2]})")
        else:
            logger.info("Validation PASSED: all pairs have distance >= k")
    else:
        logger.info("Skipping full validation (use --validate to enable)")
        
        # Quick sample validation
        sample_size = min(100, len(barcodes))
        sample = barcodes[:sample_size]
        failures = validate_barcode_set(sample, args.distance, max_to_check=10000)
        
        if failures:
            logger.warning(f"Sample validation found {len(failures)} errors")
        else:
            logger.info(f"Sample validation PASSED ({sample_size} barcodes)")
    
    # Save results - use parameter names matching CLI args (-l -> L, -d -> D)
    filename = f"barcodes_L{args.length}_D{args.distance}_{len(barcodes)}.csv"
    
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['barcode', 'gc_content'])
        for bc in barcodes:
            gc = (bc.count('G') + bc.count('C')) / len(bc)
            writer.writerow([bc, f"{gc:.3f}"])
    
    total_time = time.time() - total_start
    logger.info("=" * 70)
    logger.info(f"Saved to {filename}")
    logger.info(f"Total time: {total_time:.1f}s")
    logger.info(f"Output: {len(barcodes):,} barcodes")
    logger.info("=" * 70)

if __name__ == '__main__':
    main()
