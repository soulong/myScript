#!/usr/bin/env python
"""
Tile Splitter for Microscopy Images
Splits every TIFF image into non-overlapping tiles of given pixel dimensions.
Works recursively on multiple Measurement folders.
"""

import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Optional
import numpy as np
from tifffile import imread, imwrite
from tqdm import tqdm

from image_helper import setup_logger, find_measurement_dirs


def split_image_into_tiles(
    image_path: Path,
    tile_w_px: int,
    tile_h_px: int,
    logger: logging.Logger,
) -> Tuple[int, List[str]]:
    """
    Split a single image into tiles of exactly tile_w_px × tile_h_px pixels.
    The last tiles on right/bottom are cropped if image size is not divisible.
    Returns (total_tiles_created, list_of_saved_filenames)
    """
    try:
        img = imread(image_path)
        if img.ndim != 2:
            logger.warning(f"Skipping {image_path.name}: expected 2D image, got {img.ndim}D")
            return 0, []

        if img.dtype != np.uint16:
            img = img.astype(np.uint16)

        height, width = img.shape

        if tile_w_px <= 0 or tile_h_px <= 0:
            logger.error(f"Invalid tile size {tile_w_px}×{tile_h_px} for {image_path.name}")
            return 0, []

        saved_files = []
        tile_idx = 0

        for y_start in range(0, height, tile_h_px):
            for x_start in range(0, width, tile_w_px):
                y_end = min(y_start + tile_h_px, height)
                x_end = min(x_start + tile_w_px, width)

                tile = img[y_start:y_end, x_start:x_end]

                # Pad tile name with row/col index for easier sorting later if desired
                row = y_start // tile_h_px
                col = x_start // tile_w_px
                stem = image_path.stem
                new_name = f"{stem}--t{tile_idx:04d}_r{row:03d}c{col:03d}.tif"
                output_path = image_path.parent / new_name

                imwrite(output_path, tile, compression='zlib')
                saved_files.append(new_name)
                tile_idx += 1

        return len(saved_files), saved_files

    except Exception as e:
        logger.error(f"Failed to process {image_path.name}: {e}")
        return 0, []


def tile_images_in_directory(
    directory: Path,
    tile_w_px: int = 512,
    tile_h_px: int = 512,
    delete_originals: bool = False,
    logger: Optional[logging.Logger] = None,
) -> dict:
    """
    Process all .tiff/.tif files in one Images directory using fixed tile pixel size.
    """
    if logger is None:
        logger = logging.getLogger("tile")

    directory = Path(directory)
    if not directory.is_dir():
        logger.error(f"Directory not found: {directory}")
        return {"processed": 0, "tiles_created": 0, "skipped": 0, "saved": [], "deleted": [], "errors": []}

    tiff_files = list(directory.glob("*.tif")) + list(directory.glob("*.tiff"))
    if not tiff_files:
        logger.info(f"No TIFF images found in {directory}")
        return {"processed": 0, "tiles_created": 0, "skipped": 0, "saved": [], "deleted": [], "errors": []}

    summary = {
        "processed": 0,
        "tiles_created": 0,
        "skipped": 0,
        "saved": [],
        "deleted": [],
        "errors": []
    }

    progress = tqdm(tiff_files, desc="Tiling Images", unit="image")

    for img_path in progress:
        progress.set_postfix({"file": img_path.name[:30], "tile": f"{tile_w_px}×{tile_h_px}"})

        n_tiles, tile_names = split_image_into_tiles(img_path, tile_w_px, tile_h_px, logger)

        if n_tiles == 0:
            summary["skipped"] += 1
            continue

        summary["processed"] += 1
        summary["tiles_created"] += n_tiles
        summary["saved"].extend(tile_names)

        if delete_originals:
            try:
                img_path.unlink()
                summary["deleted"].append(img_path.name)
            except Exception as e:
                msg = f"Failed to delete original {img_path.name}: {e}"
                logger.error(msg)
                summary["errors"].append(msg)
        else:
            summary["deleted"].append(f"(kept) {img_path.name}")

    action = "Deleted" if delete_originals else "Kept"
    logger.info(
        f"Tiling complete in {directory.name} → "
        f"Images: {summary['processed']:,} | "
        f"Tiles ({tile_w_px}×{tile_h_px}px): {summary['tiles_created']:,} | "
        f"Skipped: {summary['skipped']:,} | "
        f"{action} originals"
    )
    if summary["errors"]:
        logger.warning(f"{len(summary['errors'])} errors occurred.")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Split microscopy TIFF images into fixed-size tiles (pixel dimensions)",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("root_dir", type=str, help="Root directory containing Measurement folders")
    parser.add_argument("-w", "--tile-width", type=int, default=1024, help="Tile width in pixels")
    parser.add_argument("-h", "--tile-height", type=int, default=1024, help="Tile height in pixels")
    parser.add_argument("--delete-originals", action="store_true", help="Delete original full images after successful tiling")
    args = parser.parse_args()

    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    if args.tile_width <= 0 or args.tile_height <= 0:
        print("Error: tile-width and tile-height must be positive integers")
        return
    
    root_dir = Path(args.root_dir)
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="tile_splitter")
    logger.info("Starting tile splitting")

    measurements = find_measurement_dirs(root_dir)
    if not measurements:
        logger.error("No Measurement directories found.")
        return

    logger.info(f"Found {len(measurements)} Measurement directories.")
    logger.info(f"Splitting each image into {args.tile_width} x {args.tile_height} tiles")

    global_summary = {
        "measurements_processed": 0,
        "total_images": 0,
        "total_tiles": 0}

    for mdir in measurements:
        logger.info(f"Processing Measurement: {mdir.name}")
        img_dir = mdir / "Images"
        if not img_dir.exists():
            logger.warning(f"Images directory not found: {img_dir}")
            continue

        result = tile_images_in_directory(
                    img_dir,
                    tile_w_px=args.tile_width,
                    tile_h_px=args.tile_height,
                    delete_originals=args.delete_originals,
                    logger=logger)

        global_summary["measurements_processed"] += 1
        global_summary["total_images"] += result["processed"]
        global_summary["total_tiles"] += result["tiles_created"]

    logger.info(
        f"All done! Processed {global_summary['measurements_processed']} Measurement folders | "
        f"{global_summary['total_images']:,} images → {global_summary['total_tiles']:,} tiles "
        f"({args.tile_width}×{args.tile_height}px) created")


if __name__ == "__main__":
    main()