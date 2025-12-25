#!/usr/bin/env python
"""
run_ilastik.py
--------------
Run ilastik in headless mode **in parallel** over all images in every
``*_Measurement X`` folder.

Features
--------
* Discovers measurements automatically
* Splits images into batches → one ilastik process per batch
* Each batch logs to ``<root>/logs/ilastik/<measurement>/batch_XXX.log``
* Respects existing output (skip unless --overwrite)
* Optional subset pattern & test mode
* Full CLI with defaults (``--help``)
"""

from __future__ import annotations

import argparse
import math
import random
import logging
import re
import sys
import time
from datetime import datetime
from natsort import natsorted
from pathlib import Path
from subprocess import Popen
from typing import List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from image_helper import setup_logger, find_measurement_dirs



# ----------------------------------------------------------------------- #
# Single ilastik job (one batch of images)
# ----------------------------------------------------------------------- #
def _run_ilastik_batch(
    image_paths: List[Path],
    ilastik_proj: Path,
    ilastik_exec: Path,
    log_path: Path,
) -> int:
    """
    Execute ilastik on a **list** of image paths.

    Parameters
    ----------
    image_paths
        Full paths to input TIFFs.
    ilastik_proj
        Path to ``.ilp`` project.
    ilastik_exec
        Path to ``run_ilastik.sh`` (or executable).
    log_path
        Where to write stdout+stderr.

    Returns
    -------
    int
        Process exit code.
    """
    if not image_paths:
        return 0

    # Build command: "path/to/run_ilastik.sh" --headless --project="proj.ilp" "img1.tiff" "img2.tiff" ...
    quoted_imgs = " ".join(f'"{p}"' for p in image_paths)
    cmd = f'"{ilastik_exec}" --headless --project="{ilastik_proj}" {quoted_imgs}'

    start = datetime.now()
    print(f"[{start}] START batch → {len(image_paths)} images | log: {log_path.name}")

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", buffering=1, encoding="utf-8") as log_file:
        log_file.write(f"=== ilastik START: {start} ===\n")
        log_file.write(f"CMD: {cmd}\n\n")
        log_file.flush()

        try:
            proc = Popen(
                cmd,
                shell=True,
                text=True,
                stdout=log_file,
                stderr=log_file,
            )
            proc.communicate()  # blocks until done

            end = datetime.now()
            duration = end - start
            status = "SUCCESS" if proc.returncode == 0 else f"FAILED (code {proc.returncode})"
            log_file.write(f"\n=== ilastik END: {end} | Duration: {duration} | {status} ===\n")
            print(f"[{end}] END batch → {status}")
            return proc.returncode

        except Exception as e:
            error = f"EXCEPTION: {e}\n"
            log_file.write(error)
            print(f"[ERROR] {error}")
            return -1



# ----------------------------------------------------------------------- #
# Process one measurement folder
# ----------------------------------------------------------------------- #
def ilastik_process_measurement(
    measurement_dir: Path,
    image_subdir: str,
    ilastik_proj: Path,
    ilastik_exec: Path,
    output_suffix: str,
    subset_pattern: Optional[str],
    image_suffix: str,
    max_workers: int,
    overwrite: bool,
    test: bool,
    logger: logging.Logger,
) -> bool:
    """
    Run ilastik on **one** measurement.

    Returns
    -------
    bool
        ``True`` if all batches succeeded.
    """
    img_dir = measurement_dir / image_subdir
    if not img_dir.is_dir():
        logger.error(f"Images directory missing: {img_dir}")
        return False

    # ------------------------------------------------------------------- #
    # Gather input files
    # ------------------------------------------------------------------- #
    all_files = natsorted(img_dir.glob(f"*{image_suffix}"))
    # Exclude already-generated probability maps
    input_files = [
        p for p in all_files
        if output_suffix not in p.name
    ]

    if subset_pattern:
        pattern = re.compile(subset_pattern)
        input_files = [p for p in input_files if pattern.search(p.name)]

    if not input_files:
        logger.info(f"No images to process in {measurement_dir.name}")
        return True

    # ------------------------------------------------------------------- #
    # Skip existing outputs unless overwrite
    # ------------------------------------------------------------------- #
    if not overwrite:
        prob_paths = {p.parent / f"{p.stem}_{output_suffix}.tiff" for p in input_files}
        existing = [p for p in prob_paths if p.exists()]
        to_skip = {p.with_name(p.name.replace(f"_{output_suffix}.tiff", ".tiff")) for p in existing}
        input_files = [p for p in input_files if p not in to_skip]
        logger.info(f"Skipping {len(existing)} existing probability maps")

    if test:
        input_files = random.sample(input_files, min(4, len(input_files)))
        logger.info(f"TEST MODE: processing {len(input_files)} random images")

    if not input_files:
        logger.info("Nothing to do after filtering")
        return True

    # ------------------------------------------------------------------- #
    # Split into batches
    # ------------------------------------------------------------------- #
    batch_size = math.ceil(len(input_files) / max_workers)
    batches: List[List[Path]] = [
        input_files[i:i + batch_size]
        for i in range(0, len(input_files), batch_size)
    ]
    logger.info(f"Splitting {len(input_files)} images into {len(batches)} batches "
                f"(batch_size ≈ {batch_size})")

    # ------------------------------------------------------------------- #
    # Run in parallel
    # ------------------------------------------------------------------- #
    log_dir = measurement_dir / "logs" / "ilastik"
    failed = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {
            executor.submit(
                _run_ilastik_batch,
                image_paths=batch,
                ilastik_proj=ilastik_proj,
                ilastik_exec=ilastik_exec,
                log_path=log_dir / f"batch_{idx:03d}.log",
            ): idx
            for idx, batch in enumerate(batches)
        }

        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            code = future.result()
            if code != 0:
                failed.append(idx)

    # ------------------------------------------------------------------- #
    # Summary
    # ------------------------------------------------------------------- #
    total_time = sum((end - start).total_seconds() for start, end in [
        (datetime.now(), datetime.now())  # placeholder
    ])  # not tracked per batch, but we log per‑batch duration

    logger.info(
        f"ilastik complete → batches: {len(batches)} | "
        f"failed: {len(failed)} | images: {len(input_files)}"
    )
    if failed:
        logger.error("Failed batches: " + ", ".join(f"batch_{i:03d}" for i in failed))
        return False
    return True


# ----------------------------------------------------------------------- #
# CLI
# ----------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(description="Run ilastik (headless) in parallel over multiple Measurement folders",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("root_dir", type=str, help="Root folder containing *_Measurement X folders")
    parser.add_argument("--ilastik-proj", type=str, required=True, help="Path to ilastik project file (.ilp)")
    parser.add_argument("--ilastik-exec", type=str, default="/home/haohe/ilastik-1.4.2b5-Linux/run_ilastik.sh", help="Path to ilastik launcher script (e.g. run_ilastik.sh)")
    parser.add_argument("--image-subdir", type=str, default="Images", help="Subfolder containing raw TIFFs")
    parser.add_argument("--image-suffix", type=str, default=".tiff", help="File extension of input images")
    parser.add_argument("--output-suffix", type=str, default="Probabilities", help="Suffix added to output probability maps")
    parser.add_argument("--subset-pattern", type=str, default=None, help="Regex pattern to select subset of images (applied to filename)")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel ilastik processes")
    parser.add_argument("--overwrite", action="store_true", help="Reprocess images even if output exists")
    parser.add_argument("--test", action="store_true", help="Process only 4 random images per measurement (for debugging)")
    args = parser.parse_args()

    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    root_dir = Path(args.root_dir)
    if not root_dir.is_dir():
        print(f"Root directory not found: {root_dir}")
        sys.exit(1)

    ilastik_proj = Path(args.ilastik_proj)
    if not ilastik_proj.exists():
        print(f"ilastik project not found: {ilastik_proj}")
        sys.exit(1)

    ilastik_exec = Path(args.ilastik_exec)
    if not ilastik_exec.exists():
        print(f"ilastik executable not found: {ilastik_exec}")
        sys.exit(1)

    # Logging (shared with other modules)
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="ilastik")
    logger.info(f"ilastik pipeline started – root: {root_dir}")
    logger.info(f"Project: {ilastik_proj} | Exec: {ilastik_exec} | Workers: {args.workers}")

    # Discover measurements
    measurements = find_measurement_dirs(root_dir)
    if not measurements:
        logger.error("No *_Measurement X folders found")
        sys.exit(1)
    logger.info(f"Found {len(measurements)} measurement directories")

    # ------------------------------------------------------------------- #
    # Run per measurement
    # ------------------------------------------------------------------- #
    overall_success = True
    global_start = time.time()

    for mdir in measurements:
        logger.info(f"--- Processing {mdir.name} ---")
        try:
            ok = ilastik_process_measurement(
                measurement_dir=mdir,
                image_subdir=args.image_subdir,
                ilastik_proj=ilastik_proj,
                ilastik_exec=ilastik_exec,
                output_suffix=args.output_suffix,
                subset_pattern=args.subset_pattern,
                image_suffix=args.image_suffix,
                max_workers=args.workers,
                overwrite=args.overwrite,
                test=args.test,
                logger=logger,
            )
            if not ok:
                overall_success = False
        except Exception as e:
            logger.error(f"Unexpected error in {mdir.name}: {e}", exc_info=True)
            overall_success = False

    total_time = time.time() - global_start
    logger.info(f"ilastik pipeline finished in {total_time:.1f}s")
    sys.exit(0 if overall_success else 1)


if __name__ == "__main__":
    main()