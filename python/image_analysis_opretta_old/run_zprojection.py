
#!/usr/bin/env python
"""
MIP (Maximum Intensity Projection) Processor
Converts Z-stacks → single-plane MIP images.
Works recursively on multiple Measurement folders.
"""

import argparse
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from tifffile import imread, imwrite
from tqdm import tqdm

from image_helper import setup_logger, find_measurement_dirs


def mip_stacks_in_directory(
    directory: Path,
    delete_originals: bool = True,
    min_stack_size: int = 2,
    logger: Optional[logging.Logger] = None,
) -> Dict:
    """Perform MIP on Z-stacks grouped by (row, col, field, channel)."""
    if logger is None:
        logger = logging.getLogger("mip")

    directory = Path(directory)
    if not directory.is_dir():
        logger.error(f"Directory not found: {directory}")
        return {"processed": 0, "skipped_single": 0, "skipped_error": 0, "saved": [], "deleted": [], "errors": []}

    pattern = re.compile(r'r(\d+)c(\d+)f(\d+)p(\d+)-ch(\d+)sk\d+fk1fl1\.tiff')
    groups = defaultdict(list)

    # Group files
    for file_path in directory.glob("*.tiff"):
        match = pattern.match(file_path.name)
        if not match:
            continue
        r, c, f, p, ch = map(int, match.groups())
        groups[(r, c, f, ch)].append((p, file_path))

    if not groups:
        logger.info("No valid image stacks found.")
        return {"processed": 0, "skipped_single": 0, "skipped_error": 0, "saved": [], "deleted": [], "errors": []}

    summary = {
        "processed": 0,
        "skipped_single": 0,
        "skipped_error": 0,
        "saved": [],
        "deleted": [],
        "errors": []}

    progress_bar = tqdm(groups.items(), desc="MIP Processing", unit="group")

    for (r, c, f, ch), stack in progress_bar:
        progress_bar.set_postfix({"r": f"{r:02d}", "c": f"{c:02d}", "f": f"{f:02d}", "ch": ch, "n": len(stack)})

        if len(stack) < min_stack_size:
            summary["skipped_single"] += 1
            continue

        stack.sort(key=lambda x: x[0])  # sort by plane
        paths = [p for _, p in stack]

        images = []
        for path in paths:
            try:
                img = imread(path)
                if img.dtype != np.uint16:
                    img = img.astype(np.uint16)
                images.append(img)
            except Exception as e:
                msg = f"Read failed: {path.name} → {e}"
                logger.error(msg)
                summary["errors"].append(msg)
                images = None
                break

        if images is None or len(images) == 0:
            summary["skipped_error"] += 1
            continue

        try:
            mip = np.max(np.stack(images, axis=0), axis=0).astype(np.uint16)
        except Exception as e:
            msg = f"MIP computation failed r{r:02d}c{c:02d}f{f:02d}ch{ch}: {e}"
            logger.error(msg)
            summary["errors"].append(msg)
            summary["skipped_error"] += 1
            continue

        new_name = f"r{r:02d}c{c:02d}f{f:02d}p00-ch{ch}sk1fk1fl1.tiff"
        output_path = directory / new_name

        try:
            imwrite(output_path, mip, compression='zlib')
            summary["saved"].append(new_name)
        except Exception as e:
            msg = f"Save failed: {new_name} → {e}"
            logger.error(msg)
            summary["errors"].append(msg)
            summary["skipped_error"] += 1
            continue

        if delete_originals:
            for path in paths:
                try:
                    path.unlink()
                    summary["deleted"].append(path.name)
                except Exception as e:
                    msg = f"Delete failed: {path.name} → {e}"
                    logger.error(msg)
                    summary["errors"].append(msg)
        else:
            summary["deleted"].extend([f"(kept) {p.name}" for p in paths])

        summary["processed"] += 1

    logger.info(
        f"MIP Complete → Processed: {summary['processed']:,} | "
        f"Skipped (single plane): {summary['skipped_single']:,} | "
        f"Saved: {len(summary['saved']):,} | "
        f"{'Deleted' if delete_originals else 'Kept'}: {len(summary['deleted']):,}"
    )
    if summary["errors"]:
        logger.warning(f"{len(summary['errors']):,} errors occurred.")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Z-stack → MIP (Maximum Intensity Projection)", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("root_dir", type=str, help="Root directory containing Measurement folders")
    parser.add_argument("--min-stack-size", type=int, default=2, help="Minimum number of planes required to perform MIP")
    parser.add_argument("--delete-originals", action="store_true", help="Do not delete original Z-stack files after MIP")
    args = parser.parse_args()
    
    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    root_dir = Path(args.root_dir)
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="mip_projection")
    logger.info("Starting MIP projection")

    measurements = find_measurement_dirs(root_dir)
    if not measurements:
        logger.error("No Measurement directories found.")
        return

    logger.info(f"Found {len(measurements)} Measurement directories.")

    for mdir in measurements:
        logger.info(f"Processing Measurement: {mdir.name}")
        img_dir = mdir / "Images"
        if not img_dir.exists():
            logger.warning(f"Images directory not found: {img_dir}")
            continue

        mip_stacks_in_directory(
            img_dir,
            delete_originals=args.delete_originals,
            min_stack_size=args.min_stack_size,
            logger=logger)

    logger.info("MIP processing finished.")


if __name__ == "__main__":
    main()