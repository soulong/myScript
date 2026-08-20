#!/usr/bin/env python3
# MIT License
# Copyright (c) 2020 David Jacob Feldman
# v6: merged rewrite of the dna_barcode family.
#   - exhaustive mode (n <= 10)  : full enumeration, khash buckets, sparse
#     conflict graph, greedy max-clique  (from dna_barcode.py, v4 fixes)
#   - sampling mode  (n > 10)    : parallel workers generate & hard-filter
#     candidates into a queue; the main process owns a central khash index
#     and validates each candidate exactly, so --limit is the number of
#     validated barcodes in the output
#   - mandatory validation pass, C Levenshtein when available
#
# Guarantee: every barcode in the output is at Levenshtein/Hamming distance
# >= k from every other barcode in the output. See the khash() docstring for
# the hash-correctness argument (windows are registered under enough position
# labels to cover any alignment offset up to k-1).

"""
DNA Barcode Generator v6

Designs DNA barcodes with a guaranteed pairwise edit distance >= k.

Modes (selected automatically by length):
    n <= 10  exhaustive enumeration -> maximum-density greedy clique
    n > 10   parallel random sampling into a queue; the main process holds
             a central khash index and validates each candidate exactly

Distance metric:
    Hamming by default (fast); add --use_levenshtein for full edit distance.

Validation:
    A validation pass always runs after generation. If the set is small
    enough (n*(n-1)/2 <= max_to_check pairs) it is exhaustive; otherwise a
    random subset of barcodes is checked pairwise and the log says so.

Usage:
    python dna_barcode_v6.py --length 12 --distance 3 --limit 1000
    python dna_barcode_v6.py --length 8 --distance 2 --workers 4
    python dna_barcode_v6.py --length 16 --distance 4 --limit 5000 \\
        --use_levenshtein --hairpin_max 4

Output columns (rows are sorted by mean_levenshtein desc, then mean_hamming desc):
    barcode          the sequence
    gc_content       fraction of G/C
    gc_center_diff   |GC - 0.5| (lower is better)
    shannon_entropy  sequence complexity (max ~2.0)
    hairpin_stem     longest reverse-complement match
    min_hamming      min pairwise Hamming distance to any other barcode
    min_levenshtein  min pairwise Levenshtein distance to any other barcode
    mean_hamming     mean pairwise Hamming distance to all other barcodes
    mean_levenshtein mean pairwise Levenshtein distance to all other barcodes
"""

import argparse
import logging
import math
import multiprocessing
import os
import random
import re
import time
from collections import Counter, defaultdict
from queue import Empty

import numpy as np

logger = logging.getLogger(__name__)

# ============================================================================
# Optional dependencies
# ============================================================================

try:
    import Levenshtein as _c_lev
except ImportError:
    _c_lev = None

try:
    import scipy.sparse as _sparse
except ImportError:
    _sparse = None

try:
    from tqdm import tqdm as _tqdm
except ImportError:
    _tqdm = None

# ============================================================================
# Constants
# ============================================================================

EXHAUSTIVE_MAX_LENGTH = 10      # 4^10 = 1,048,576 sequences
DEFAULT_MAX_TO_CHECK = int(1e6)
DEFAULT_MAX_ATTEMPTS = int(2e6)
INT_TO_NUC = ['A', 'C', 'G', 'T']
COMPLEMENT = str.maketrans("ACTG", "TGAC")

# ============================================================================
# Distance functions (cutoff semantics: returns exact d if d <= cutoff,
# otherwise cutoff + 1)
# ============================================================================


def hamming_cutoff(a, b, cutoff):
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > cutoff:
                return cutoff + 1
    return d


def levenshtein_cutoff(s1, s2, cutoff):
    len1, len2 = len(s1), len(s2)
    if abs(len1 - len2) > cutoff:
        return cutoff + 1
    prev = list(range(len2 + 1))
    curr = [0] * (len2 + 1)
    for i in range(1, len1 + 1):
        curr[0] = i
        row_min = i
        for j in range(1, len2 + 1):
            if s1[i - 1] == s2[j - 1]:
                v = prev[j - 1]
            else:
                v = 1 + min(prev[j - 1], prev[j], curr[j - 1])
            curr[j] = v
            if v < row_min:
                row_min = v
        # every alignment path crosses row i, so if the whole row exceeds the
        # cutoff the final distance must too
        if row_min > cutoff:
            return cutoff + 1
        prev, curr = curr, prev
    return prev[len2] if prev[len2] <= cutoff else cutoff + 1


def make_distance(use_levenshtein, k):
    """Return dist(a, b) with cutoff semantics for threshold k (exact below k)."""
    cutoff = k - 1
    if use_levenshtein:
        if _c_lev is not None:
            return lambda a, b: min(_c_lev.distance(a, b), k)
        return lambda a, b: levenshtein_cutoff(a, b, cutoff)
    return lambda a, b: hamming_cutoff(a, b, cutoff)


# ============================================================================
# Sequence filters
# ============================================================================


def calculate_gc(s):
    return (s.count('G') + s.count('C')) / len(s)


def gc_bounds(n, gc_min, gc_max):
    """Integer GC-count bounds that exactly mirror the float comparisons
    `gc_min <= gc <= gc_max` on rational GC fractions."""
    return (math.ceil(gc_min * n - 1e-9),
            math.floor(gc_max * n + 1e-9))


def count_gc_int(num, length):
    """GC count from 2-bit encoding (C=01, G=10 -> bits 1 and 2)."""
    gc = 0
    for _ in range(length):
        bits = num & 3
        if bits == 1 or bits == 2:
            gc += 1
        num >>= 2
    return gc


def has_homopolymer_int(num, length, homopolymer):
    """True if any run of identical nucleotides has length >= homopolymer."""
    if length < homopolymer:
        return False
    run = 1
    prev = num & 3
    num >>= 2
    for _ in range(1, length):
        cur = num & 3
        if cur == prev:
            run += 1
            if run >= homopolymer:
                return True
        else:
            run = 1
        prev = cur
        num >>= 2
    return False


def int_to_seq(num, length):
    """Least significant 2 bits encode the first nucleotide."""
    out = [None] * length
    for i in range(length):
        out[i] = INT_TO_NUC[num & 3]
        num >>= 2
    return ''.join(out)


def random_barcode_int(rng, length, target_gc):
    """Random barcode biased toward target_gc, as 2-bit integer."""
    p_at = (1 - target_gc) / 2
    p_gc = target_gc / 2
    weights = (p_at, p_gc, p_gc, p_at)
    num = 0
    for i in range(length):
        num |= rng.choices((0, 1, 2, 3), weights=weights)[0] << (2 * i)
    return num


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def longest_complementary_match_fast(seq, threshold=None):
    """Longest reverse-complement match (hairpin stem) via DP, O(n^2)."""
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
                    if threshold is not None and max_len > threshold:
                        return max_len
            else:
                curr[j] = 0
        prev = curr
    return max_len


# ============================================================================
# k-hash
# ============================================================================

def window_size(n, k):
    """Window length that guarantees no false negatives (see docstring below)."""
    return max(1, int(np.ceil((n - k + 1) / float(k))))


def khash(s, k, spread=1):
    """
    Return hash keys (position, substring) for s.

    Soundness: two same-length strings with edit distance < k share a
    substring of length W = ceil((n-k+1)/k). Proof: an alignment with
    <= k-1 edit columns leaves >= n-k+1 matching columns, split by the
    edits into at most k contiguous blocks, so the longest block has
    length >= ceil((n-k+1)/k) = W. Within a block both strings contain
    the same substring, and the block starts at positions differing by
    the alignment offset o, bounded by the number of indels before the
    block, i.e. |o| <= k-1. Each window is therefore registered under
    the 2*spread+1 labels around its start, so a window shared at offset
    o is keyed at the common label (i+o) whenever |o| <= spread; the
    circular extension (s+s) covers windows that wrap around the end.

    Use spread = k-1 for edit-distance (Levenshtein) filtering so every
    pair at distance < k shares a key (|o| <= k-1 always holds); the
    default spread = 1, i.e. labels {i, i+1}, suffices for Hamming
    distance, where the offset is always 0.
    """
    n = len(s)
    w = window_size(n, k)
    s2 = s + s
    keys = []
    seen = set()
    for i in range(n):
        sub = s2[i:i + w]
        for t in range(-spread, spread + 1):
            key = ((i + t) % n, sub)
            if key not in seen:
                seen.add(key)
                keys.append(key)
    return keys


# ============================================================================
# Exhaustive mode (n <= 10)
# ============================================================================

def generate_all_barcodes_exhaustive(n, gc_min, gc_max, homopolymer,
                                     exclude_pattern):
    total = 4 ** n
    gc_min_count, gc_max_count = gc_bounds(n, gc_min, gc_max)
    barcodes = []
    for num in range(total):
        gc = count_gc_int(num, n)
        if gc < gc_min_count or gc > gc_max_count:
            continue
        if has_homopolymer_int(num, n, homopolymer):
            continue
        seq = int_to_seq(num, n)
        if exclude_pattern is not None and exclude_pattern.search(seq):
            continue
        barcodes.append(seq)
    return barcodes


def build_khash_dict(xs, k, spread=1):
    D = defaultdict(list)
    for x in xs:
        for h in khash(x, k, spread=spread):
            D[h].append(x)
    return {key: sorted(set(v)) for key, v in D.items()}


def sparse_dist(hash_buckets, k, use_levenshtein, progress=None):
    """Pairwise distances within buckets; keeps only conflicts (d < k)."""
    if progress is None:
        progress = lambda x: x
    dist = make_distance(use_levenshtein, k)
    D = {}
    for xs in progress(hash_buckets):
        for i, a in enumerate(xs):
            for b in xs[i + 1:]:
                if dist(a, b) < k:
                    key = tuple(sorted((a, b)))
                    D[key] = True
    return D


def _bucket_worker(args):
    bucket, k, use_levenshtein = args
    return sparse_dist([bucket], k, use_levenshtein)


def sparse_view(xs, D, symmetric=True):
    if len(xs) != len(set(xs)):
        raise ValueError("sparse_view requires unique sequences in xs")
    mapper = {x: i for i, x in enumerate(xs)}
    if len(D) == 0:
        i, j, data = [], [], []
    else:
        i, j, data = zip(*[(mapper[a], mapper[b], True) for (a, b) in D])
        i = np.array(i)
        j = np.array(j)
        data = np.array(data)
    n = len(xs)
    if _sparse is None:
        raise RuntimeError("scipy is required for exhaustive mode "
                           "(pip install scipy)")
    cm = _sparse.coo_matrix((data, (i, j)), shape=(n, n))
    if symmetric:
        cm = (cm + cm.T).tocsr()
    return cm


def maxy_clique_groups(cm, group_ids, verbose=False, rng=None):
    """
    Greedy clique on the conflict graph, prioritizing groups with the
    fewest remaining candidates. Returns selected indices.

    Passing an rng randomizes the order in which equal-sized candidates
    are considered, which lets callers restart the greedy search and keep
    the best result.
    """
    # counts => group_id
    d1 = defaultdict(set)
    for id_, counts in Counter(group_ids).items():
        d1[counts].add(id_)

    # group_id => indices (reversed for efficient pop from the end)
    d2 = defaultdict(list)
    for i, id_ in enumerate(group_ids):
        d2[id_].append(i)
    d2 = {k: v[::-1] for k, v in d2.items()}
    if rng is not None:
        for k in d2:
            rng.shuffle(d2[k])

    selected = []
    available = np.array(range(len(group_ids)))
    available_set = set(range(len(group_ids)))

    while d1:
        if verbose and (len(selected) % 1000) == 0:
            logger.info(f"Selected {len(selected)} barcodes...")

        count = min(d1.keys())
        id_ = d1[count].pop()
        if len(d1[count]) == 0:
            d1.pop(count)

        index = None
        while d2[id_]:
            index = d2[id_].pop()
            if index in available_set:
                break
        else:
            index = None

        if index is not None:
            selected.append(index)
            available_set.discard(index)
            available = available[available != index]

            remove = cm[index, available].indices
            mask = np.ones(len(available), dtype=bool)
            mask[remove] = False
            removed_indices = available[~mask]
            for ri in removed_indices:
                available_set.discard(int(ri))
            available = available[mask]

        n = len(d2[id_])
        if n > 0:
            d1[n].add(id_)

    return selected


def run_exhaustive(n, k, gc_min, gc_max, homopolymer, exclude_pattern,
                   limit, workers, seed, use_levenshtein, progress=None):
    logger.info(f"=== EXHAUSTIVE MODE (n={n}) ===")
    t0 = time.time()
    all_barcodes = generate_all_barcodes_exhaustive(
        n, gc_min, gc_max, homopolymer, exclude_pattern)
    logger.info(f"Generated {len(all_barcodes):,} barcodes after filters "
                f"({time.time() - t0:.1f}s)")

    spread = k - 1 if use_levenshtein else 1
    buckets_dict = build_khash_dict(all_barcodes, k, spread=spread)
    buckets = list(buckets_dict.values())
    logger.info(f"Built {len(buckets):,} hash buckets")

    t0 = time.time()
    if workers > 1:
        with multiprocessing.Pool(processes=workers) as pool:
            D = {}
            it = pool.imap_unordered(
                _bucket_worker,
                ((bucket, k, use_levenshtein) for bucket in buckets))
            bar = progress(total=len(buckets)) if progress else None
            for i, d in enumerate(it):
                D.update(d)
                if bar is not None:
                    bar.update()
            if bar is not None:
                bar.close()
    else:
        D = sparse_dist(buckets, k, use_levenshtein, progress=progress)
    logger.info(f"Computed conflicts within buckets "
                f"({len(D):,} conflicts, {time.time() - t0:.1f}s)")

    cm = sparse_view(all_barcodes, D)
    logger.info("Selecting maximum clique...")
    group_ids = [0] * len(all_barcodes)
    indices = maxy_clique_groups(cm, group_ids)
    selected = [all_barcodes[i] for i in indices]
    logger.info(f"Clique found {len(selected):,} barcodes")

    rng = random.Random(seed)
    if len(selected) > limit:
        selected = rng.sample(selected, limit)
        logger.info(f"Trimmed to limit: {limit:,}")
    return selected


# ============================================================================
# Sampling mode (n > 10), parallel workers
# ============================================================================

def _worker_loop(wid, args, out_q, attempts_lock, attempts, stop):
    try:
        _worker_loop_core(wid, args, out_q, attempts_lock, attempts, stop)
    except BaseException:
        import traceback
        traceback.print_exc()
        try:
            with attempts_lock:
                stop.value = 1
        except BaseException:
            pass


def _worker_loop_core(wid, args, out_q, attempts_lock, attempts, stop):
    rng = random.Random(args.seed + wid)
    n = args.length
    hp = args.homopolymer
    gc_min_count, gc_max_count = gc_bounds(n, args.gc_min, args.gc_max)
    target_gc = (args.gc_min + args.gc_max) / 2
    exclude_pat = re.compile(args.exclude) if args.exclude else None
    cap = args.max_attempts

    # Workers only generate and hard-filter candidates, then push them onto
    # a multiprocessing.Queue. All distance bookkeeping lives in the main
    # process, so the queue path stays cheap (no Manager, no shared index).
    while not stop.value:
        num = random_barcode_int(rng, n, target_gc)
        gc = count_gc_int(num, n)
        if gc < gc_min_count or gc > gc_max_count:
            continue
        if has_homopolymer_int(num, n, hp):
            continue
        seq = int_to_seq(num, n)
        if exclude_pat is not None and exclude_pat.search(seq):
            continue

        # sequence passed all hard filters: consume from the shared budget
        with attempts_lock:
            attempts.value += 1
            if attempts.value >= cap:
                stop.value = 1
                break

        if args.hairpin_max and \
                longest_complementary_match_fast(seq, args.hairpin_max) > args.hairpin_max:
            continue

        try:
            out_q.put(seq)
        except (BrokenPipeError, EOFError, OSError):
            break


def run_sampling(args, progress=None):
    logger.info(f"=== SAMPLING MODE (n={args.length}) ===")

    ctx = multiprocessing.get_context()
    out_q = ctx.Queue(maxsize=max(100, args.workers * 8))
    attempts_lock = ctx.Lock()
    attempts = ctx.Value('i', 0)
    stop = ctx.Value('i', 0)

    procs = []
    for wid in range(args.workers):
        p = ctx.Process(
            target=_worker_loop,
            args=(wid, args, out_q, attempts_lock, attempts, stop))
        p.start()
        procs.append(p)

    n = args.length
    k = args.distance
    spread = k - 1 if args.use_levenshtein else 1
    dist = make_distance(args.use_levenshtein, k)

    # central index, plain in-process Python (no Manager needed)
    index = {}
    accepted = []
    accepted_set = set()

    bar = progress(total=args.limit) if progress else None
    last = 0

    def _all_workers_dead():
        return not any(p.is_alive() for p in procs)

    def _accept(seq):
        """Validate a candidate against the central index; accept it if it is
        at distance >= k from everything already accepted. Returns True if the
        candidate was added."""
        hashes = khash(seq, k, spread=spread)
        if seq in accepted_set:
            return False

        # the candidate is rejected on the first conflicting neighbor, so
        # iterate buckets directly (dist is deterministic; no need to
        # materialize a deduplicated neighbor set)
        for key in hashes:
            for nb in index.get(key, ()):
                if dist(seq, nb) < k:
                    return False

        accepted.append(seq)
        accepted_set.add(seq)
        for key in hashes:
            index.setdefault(key, []).append(seq)
        return True

    try:
        while len(accepted) < args.limit:
            if bar is not None:
                bar.update(len(accepted) - last)
                last = len(accepted)

            try:
                seq = out_q.get(timeout=0.2)
            except Empty:
                if not _all_workers_dead():
                    time.sleep(0.01)
                    continue
                # every worker has exited; a worker may still have enqueued
                # its final candidate right before dying, so drain whatever
                # is still in the queue before giving up
                while True:
                    try:
                        seq = out_q.get_nowait()
                    except Empty:
                        break
                    _accept(seq)
                break

            _accept(seq)
    except (BrokenPipeError, EOFError, OSError):
        logger.error("Sampling terminated abnormally: worker queue broke.")
        with attempts_lock:
            stop.value = 1
        for p in procs:
            if p.is_alive():
                p.terminate()
        for p in procs:
            p.join()
        return None

    with attempts_lock:
        stop.value = 1
    for p in procs:
        if p.is_alive():
            p.terminate()
        p.join()

    if bar is not None:
        bar.update(len(accepted) - last)
        bar.close()

    if len(accepted) < args.limit:
        logger.warning(f"Only generated {len(accepted):,}/{args.limit:,} "
                       f"barcodes (attempt cap of {args.max_attempts:,} "
                       f"reached, or the parameters are too strict / target "
                       f"infeasible; pass a larger --max_attempts)")
    return accepted


# ============================================================================
# Validation
# ============================================================================

def validate_barcode_set(barcodes, k, use_levenshtein,
                         max_to_check=DEFAULT_MAX_TO_CHECK, seed=42):
    """
    Verify all pairwise distances >= k. Full check when the number of pairs
    fits in max_to_check, otherwise a random subset of barcodes.
    """
    n = len(barcodes)
    total_pairs = n * (n - 1) // 2
    dist = make_distance(use_levenshtein, k)

    if total_pairs <= max_to_check:
        subset = barcodes
        mode = f"full ({total_pairs:,} pairs)"
    else:
        s = int((1 + math.sqrt(1 + 8 * max_to_check)) // 2)
        subset = random.Random(seed).sample(barcodes, s)
        mode = f"random subset of {s:,} barcodes ({s * (s - 1) // 2:,} pairs)"

    failures = []
    max_failures = 100
    for i, a in enumerate(subset):
        for b in subset[i + 1:]:
            d = dist(a, b)
            if d < k:
                failures.append((a, b, d))
                if len(failures) > max_failures:
                    logger.warning(f"Aborting validation: more than "
                                   f"{max_failures} failures")
                    return failures
    return failures


def pairwise_distance_stats(barcodes):
    """Return (min_hamming, min_levenshtein, mean_hamming, mean_levenshtein):
    per barcode the minimum and mean pairwise distance to every other barcode
    in the set, for both metrics independently. A singleton set yields zeros."""
    n = len(barcodes)
    if n < 2:
        return [0] * n, [0] * n, [0] * n, [0] * n

    length = len(barcodes[0])

    # Hamming: vectorized full distance matrix
    seqs = np.array([list(bc) for bc in barcodes], dtype='<U1')  # (N, L)
    diff = seqs[:, None, :] != seqs[None, :, :]                  # (N, N, L)
    hamming_mat = diff.sum(axis=2)                               # (N, N)

    tmp = hamming_mat.copy()
    np.fill_diagonal(tmp, length + 1)
    min_hamming = tmp.min(axis=1).tolist()

    mask = ~np.eye(n, dtype=bool)
    mean_hamming = (hamming_mat[mask].reshape(n, n - 1).sum(axis=1)
                    / (n - 1)).tolist()

    # Levenshtein: symmetric loop, C library when available
    if _c_lev is not None:
        lev_dist = _c_lev.distance
    else:
        max_len = max(len(b) for b in barcodes)
        lev_dist = lambda a, b: levenshtein_cutoff(a, b, max_len)
    min_lev = [length] * n
    sum_lev = [0] * n
    for i in range(n):
        best = min_lev[i]
        ai = barcodes[i]
        for j in range(i + 1, n):
            d = lev_dist(ai, barcodes[j])
            sum_lev[i] += d
            sum_lev[j] += d
            if d < best:
                best = d
            if d < min_lev[j]:
                min_lev[j] = d
        min_lev[i] = best
    mean_lev = [s / (n - 1) for s in sum_lev]

    return min_hamming, min_lev, mean_hamming, mean_lev
    return failures


# ============================================================================
# Output
# ============================================================================

def compute_metrics(barcode):
    gc = calculate_gc(barcode)
    entropy = -sum(
        p * math.log2(p)
        for p in (barcode.count(b) / len(barcode) for b in 'ACTG')
        if p > 0)
    return (round(gc, 3),
            round(abs(gc - 0.5), 3),
            round(entropy, 3),
            longest_complementary_match_fast(barcode))


def write_output(barcodes, min_hamming, min_lev, mean_hamming, mean_lev,
                 filename):
    import csv
    header = ['barcode', 'gc_content', 'gc_center_diff',
              'shannon_entropy', 'hairpin_stem',
              'min_hamming', 'min_levenshtein',
              'mean_hamming', 'mean_levenshtein']
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for bc, h, lv, mh, mlv in zip(barcodes, min_hamming, min_lev,
                                      mean_hamming, mean_lev):
            writer.writerow([bc] + list(compute_metrics(bc)) +
                            [h, lv, f"{mh:.3f}", f"{mlv:.3f}"])
    return filename


# ============================================================================
# CLI
# ============================================================================

def validate_params(n, k, gc_min, gc_max, homopolymer,
                    limit=1, workers=1, hairpin_max=0, max_to_check=1,
                    max_attempts=DEFAULT_MAX_ATTEMPTS):
    if n < 4:
        logger.warning(f"length={n} is very short, may not find enough barcodes")
    if k < 1:
        raise ValueError(f"distance={k} must be at least 1")
    if k > n:
        raise ValueError(f"distance={k} > length={n}: no valid barcodes exist")
    if gc_min < 0 or gc_min > 1:
        raise ValueError(f"gc_min={gc_min} must be in [0, 1]")
    if gc_max < 0 or gc_max > 1:
        raise ValueError(f"gc_max={gc_max} must be in [0, 1]")
    if gc_min > gc_max:
        raise ValueError(f"gc_min={gc_min} > gc_max={gc_max}")
    if homopolymer < 2:
        raise ValueError(f"homopolymer={homopolymer} < 2 rejects every barcode")
    if limit < 1:
        raise ValueError(f"limit={limit} must be at least 1")
    if workers < 1:
        raise ValueError(f"workers={workers} must be at least 1")
    if hairpin_max < 0:
        raise ValueError(f"hairpin_max={hairpin_max} must be >= 0")
    if max_to_check < 1:
        raise ValueError(f"max_to_check={max_to_check} must be at least 1")
    if max_attempts < 1:
        raise ValueError(f"max_attempts={max_attempts} must be at least 1")
    if workers > limit:
        logger.warning(f"workers={workers} > limit={limit}: most workers "
                       f"will finish idle")
    if gc_min == 0:
        logger.warning(f"gc_min=0%: 0-GC barcodes are allowed")
    elif gc_min < 1 / n:
        logger.warning(f"gc_min={gc_min:.0%} is below 1/{n}; the minimum GC "
                       f"count is 1, so every barcode ends up with GC "
                       f">= {1 / n:.0%}")


def parse_args():
    parser = argparse.ArgumentParser(
        description='Design DNA barcodes with guaranteed distance >= k',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False)
    parser.add_argument('--length', '-l', type=int, default=12,
                        help='barcode length (n)')
    parser.add_argument('--distance', '-d', type=int, default=3,
                        help='minimum pairwise distance (k)')
    parser.add_argument('--limit', type=int, default=1000,
                        help='target number of output barcodes')
    parser.add_argument('--gc_min', type=int, default=40,
                        help='minimum GC percentage (0-100)')
    parser.add_argument('--gc_max', type=int, default=60,
                        help='maximum GC percentage (0-100)')
    parser.add_argument('--homopolymer', type=int, default=3,
                        help='reject runs of >= this many identical bases')
    parser.add_argument('--exclude', type=str, default='ATG',
                        help='regex pattern to exclude from barcodes (enter '
                             'empty string to disable)')
    parser.add_argument('--seed', type=int, default=42,
                        help='random seed')
    parser.add_argument('--workers', type=int, default=None,
                        help='parallel workers (default: min(4, cpus))')
    parser.add_argument('--use_levenshtein', action='store_true',
                        help='use Levenshtein distance instead of Hamming')
    parser.add_argument('--hairpin_max', type=int, default=0,
                        help='reject barcodes with reverse-complement stem '
                             'longer than this (0 = disabled)')
    parser.add_argument('--max_to_check', type=int,
                        default=DEFAULT_MAX_TO_CHECK,
                        help='validation pair cap')
    parser.add_argument('--max_attempts', type=int,
                        default=DEFAULT_MAX_ATTEMPTS,
                        help='total sample attempts budget across all workers '
                             '(shared); raise it for dense/strict constraints')
    parser.add_argument('--verbosity', '-v', type=int, default=1,
                        help='0=warnings, 1=info, 2=debug')
    parser.add_argument('--out', type=str, default='',
                        help='output CSV path (default: auto-named)')
    return parser.parse_args()


def main():
    multiprocessing.freeze_support()
    args = parse_args()

    verbosity_map = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    logging.basicConfig(
        format='%(asctime)s -- %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=verbosity_map.get(args.verbosity, logging.INFO))

    args.gc_min /= 100
    args.gc_max /= 100
    if args.workers is None:
        args.workers = min(4, os.cpu_count() or 1)
    validate_params(args.length, args.distance,
                    args.gc_min, args.gc_max, args.homopolymer,
                    limit=args.limit, workers=args.workers,
                    hairpin_max=args.hairpin_max,
                    max_to_check=args.max_to_check,
                    max_attempts=args.max_attempts)
    args.exclude = args.exclude or None

    if args.workers > 1 and args.length > EXHAUSTIVE_MAX_LENGTH:
        logger.info(f"Using {args.workers} parallel workers")

    progress = _tqdm if (args.verbosity >= 1 and _tqdm is not None) else None

    logger.info("=" * 70)
    logger.info("DNA Barcode Generator v6")
    logger.info("=" * 70)
    logger.info(f"length={args.length} distance={args.distance} "
                f"limit={args.limit} metric={'Levenshtein' if args.use_levenshtein else 'Hamming'}")
    logger.info(f"GC {args.gc_min:.0%}-{args.gc_max:.0%} "
                f"homopolymer<{args.homopolymer} exclude={args.exclude!r}")
    logger.info("=" * 70)

    t_start = time.time()
    exclude_pat = re.compile(args.exclude) if args.exclude else None

    if args.length <= EXHAUSTIVE_MAX_LENGTH:
        barcodes = run_exhaustive(
            args.length, args.distance, args.gc_min, args.gc_max,
            args.homopolymer, exclude_pat, args.limit,
            args.workers, args.seed, args.use_levenshtein,
            progress=progress)
    else:
        barcodes = run_sampling(args, progress=progress)
        if barcodes is None:
            logger.error("Sampling failed; no output written.")
            raise SystemExit(1)

    gen_time = time.time() - t_start
    logger.info(f"Generation completed in {gen_time:.1f}s "
                f"-> {len(barcodes):,} barcodes")

    # Validation
    failures = validate_barcode_set(
        barcodes, args.distance, args.use_levenshtein,
        max_to_check=args.max_to_check, seed=args.seed)
    if failures:
        logger.error(f"!! VALIDATION FAILED: {len(failures)} errors !!")
        for a, b, d in failures[:5]:
            logger.error(f"  {a} vs {b} (d={d})")
    else:
        logger.info(f"Validation passed ({len(barcodes):,} barcodes, "
                    f"{len(barcodes) * (len(barcodes) - 1) // 2:,} pairs)")

    # Output
    min_hamming, min_lev, mean_hamming, mean_lev = \
        pairwise_distance_stats(barcodes)
    if len(barcodes) > 1:
        h_vals = sorted(min_hamming)
        l_vals = sorted(min_lev)
        logger.info(f"Pairwise min-distance (median): Hamming={h_vals[len(h_vals) // 2]} "
                    f"Levenshtein={l_vals[len(l_vals) // 2]}")
        mh_vals = sorted(mean_hamming)
        ml_vals = sorted(mean_lev)
        logger.info(f"Pairwise mean distance (median): "
                    f"Hamming={mh_vals[len(mh_vals) // 2]:.3f} "
                    f"Levenshtein={ml_vals[len(ml_vals) // 2]:.3f}")

    # rank barcodes by mean distance (descending): mean_levenshtein primary,
    # mean_hamming secondary
    order = sorted(range(len(barcodes)),
                   key=lambda i: (-mean_lev[i], -mean_hamming[i]))
    barcodes = [barcodes[i] for i in order]
    min_hamming = [min_hamming[i] for i in order]
    min_lev = [min_lev[i] for i in order]
    mean_hamming = [mean_hamming[i] for i in order]
    mean_lev = [mean_lev[i] for i in order]

    if not args.out:
        tag = 'lev' if args.use_levenshtein else 'ham'
        stamp = time.strftime('%Y%m%d%H%M')
        args.out = f"barcodes_n{args.length}_k{args.distance}_{tag}_{len(barcodes)}_{stamp}.csv"
    write_output(barcodes, min_hamming, min_lev, mean_hamming, mean_lev,
                 args.out)
    logger.info(f"Saved {len(barcodes):,} barcodes to {args.out}")
    logger.info(f"Total time: {time.time() - t_start:.1f}s")


if __name__ == '__main__':
    main()
