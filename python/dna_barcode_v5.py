#!/usr/bin/env python3

import argparse
import logging
import random
from multiprocessing import Manager, Process, cpu_count

import regex as re
import numpy as np
import pandas as pd
import Levenshtein
from tqdm import tqdm

logger = logging.getLogger(__name__)

COMPLEMENT = str.maketrans("ACTG", "TGAC")


# ----------------------------
# Basic utilities
# ----------------------------

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))


def calculate_gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def has_homopolymer(seq, max_run):
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False


def is_valid_barcode(seq, gc_min, gc_max, homopolymer, exclude_pat=None):
    gc = calculate_gc(seq)

    if not (gc_min <= gc <= gc_max):
        return False

    if has_homopolymer(seq, homopolymer):
        return False

    if exclude_pat and exclude_pat.search(seq):
        return False

    return True


# ----------------------------
# Fast stem detection (O(n²))
# ----------------------------

def longest_complementary_match_fast(seq, threshold=None):
    rc = revcomp(seq)
    n = len(seq)

    prev = [0] * (n + 1)
    max_len = 0

    for i in range(1, n + 1):
        curr = [0] * (n + 1)
        for j in range(1, n + 1):
            if seq[i - 1] == rc[j - 1]:
                curr[j] = prev[j - 1] + 1

                if curr[j] > max_len:
                    max_len = curr[j]

                    if threshold and max_len > threshold:
                        return max_len
            else:
                curr[j] = 0

        prev = curr

    return max_len


# ----------------------------
# k-hash
# ----------------------------

def khash(s, k):
    n = len(s)
    window = max(1, int(np.ceil((n - k) / float(k))))

    s2 = s + s
    arr = []

    for i in range(n):
        for j in (0, 1):
            arr.append(((i + j) % n, s2[i:i + window]))

    return arr


# ----------------------------
# Sampling
# ----------------------------

def random_barcode_gc_balanced(n, target_gc):
    probs = [
        (1 - target_gc) / 2,
        target_gc / 2,
        (1 - target_gc) / 2,
        target_gc / 2,
    ]
    return ''.join(np.random.choice(list("ACTG"), p=probs, size=n))


# ----------------------------
# Worker
# ----------------------------

def worker(
    wid,
    args,
    accepted,
    accepted_set,
    index,
    lock,
    counter,
    stop_flag,
):
    np.random.seed(args.seed + wid)
    random.seed(args.seed + wid)

    dist_func = Levenshtein.distance if args.use_levenshtein else hamming_distance
    k_eff = args.distance if args.use_levenshtein else max(1, args.distance - 1)

    exclude_pat = re.compile(args.exclude, re.IGNORECASE) if args.exclude else None
    target_gc = (args.gc_min + args.gc_max) / 2

    while True:
        if stop_flag.value:
            break

        if counter.value >= args.limit:
            stop_flag.value = 1
            break

        bc = random_barcode_gc_balanced(args.length, target_gc)

        if not is_valid_barcode(bc, args.gc_min, args.gc_max, args.homopolymer, exclude_pat):
            continue

        # early hairpin pruning
        if longest_complementary_match_fast(bc, threshold=args.hairpin_max) > args.hairpin_max:
            continue

        hashes = khash(bc, k_eff)

        with lock:
            if bc in accepted_set:
                continue

            neighbors = set()
            for h in hashes:
                neighbors.update(index.get(h, []))

            valid = True
            for existing in neighbors:
                if dist_func(bc, existing) < args.distance:
                    valid = False
                    break

            if not valid:
                continue

            accepted.append(bc)
            accepted_set[bc] = 1

            for h in hashes:
                if h not in index:
                    index[h] = []
                index[h].append(bc)

            counter.value += 1


# ----------------------------
# Parallel driver
# ----------------------------

def parallel_barcode_generation(args):

    manager = Manager()

    accepted = manager.list()
    accepted_set = manager.dict()
    index = manager.dict()

    lock = manager.Lock()
    counter = manager.Value("i", 0)
    stop_flag = manager.Value("i", 0)

    n_workers = args.workers or cpu_count()

    procs = []
    for i in range(n_workers):
        p = Process(
            target=worker,
            args=(i, args, accepted, accepted_set, index, lock, counter, stop_flag),
        )
        p.start()
        procs.append(p)

    pbar = tqdm(total=args.limit, desc="Accepted", unit="barcode")

    last = 0
    while counter.value < args.limit:
        current = counter.value
        pbar.update(current - last)
        last = current

    pbar.close()

    stop_flag.value = 1

    for p in procs:
        p.join()

    return list(accepted)


# ----------------------------
# Scoring (cleaned)
# ----------------------------

def compute_scores(barcodes):
    records = []

    for seq in barcodes:
        gc = calculate_gc(seq)
        hp = longest_complementary_match_fast(seq)

        entropy = -sum(
            p * np.log2(p)
            for p in [(seq.count(b) / len(seq)) for b in "ACTG"]
            if p > 0
        )

        records.append({
            "barcode": seq,
            "gc_content": round(gc, 3),
            "gc_deviation": round(abs(gc - 0.5), 3),
            "shannon_entropy": round(entropy, 3),
            "hairpin_stem": hp,  # integer
        })

    return pd.DataFrame(records)


# ----------------------------
# CLI
# ----------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parallel DNA barcode generator (clean metrics)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""

Output Metrics:

gc_content:
    Fraction of G/C bases (recommended: 0.40–0.60)

gc_deviation:
    |GC - 0.5| (lower is better)

shannon_entropy:
    Sequence complexity (max ≈ 2.0)
    Recommended: > 1.8

hairpin_stem:
    Longest reverse-complement match
    Recommended: ≤ 3 (strict), ≤ 4 acceptable

NOTE:
All numeric outputs are rounded to 3 decimal places.
"""
    )

    parser.add_argument("--length", type=int, default=12,
                            help="Barcode length")

    parser.add_argument("--distance", type=int, default=4,
                        help="Minimum edit/Hamming distance")

    parser.add_argument("--limit", type=int, default=1000,
                        help="Number of barcodes to generate")

    parser.add_argument("--gc_min", type=float, default=0.35,
                        help="Minimum GC fraction")

    parser.add_argument("--gc_max", type=float, default=0.65,
                        help="Maximum GC fraction")

    parser.add_argument("--homopolymer", type=int, default=3,
                        help="Max allowed homopolymer length")

    parser.add_argument("--exclude", type=str,
                        default="ATG|TAA|TGA|TAG",
                        help="Regex patterns to exclude")

    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")

    parser.add_argument("--use_levenshtein", action="store_true",
                        help="Use Levenshtein distance instead of Hamming")

    parser.add_argument("--hairpin_max", type=int, default=5,
                        help="Max allowed hairpin stem length")

    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel workers (default: CPU cores)")

    return parser.parse_args()


# ----------------------------
# Main
# ----------------------------

def main():
    args = parse_args()

    logging.basicConfig(level=logging.INFO)

    logger.info("Running with parameters:")
    for k, v in vars(args).items():
        logger.info(f"{k} = {v}")

    barcodes = parallel_barcode_generation(args)

    df = compute_scores(barcodes)

    filename = f"barcodes_n{args.length}_k{args.distance}.csv"
    df.to_csv(filename, index=False)

    logger.info(f"Saved {len(df)} barcodes to {filename}")


if __name__ == "__main__":
    main()