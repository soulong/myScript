#!/usr/bin/env python
"""
run_dataloader.py
-----------------
Generate `cp_dataloader.csv` (CellProfiler input) for every ``*_Measurement X`` folder.

* Uses `images_to_dataset(..., cellprofiler_style=True)` from `image_helper.py`
* Writes CSV to `<measurement>/cp_dataloader.csv`
* Skips if file exists (unless `--overwrite`)
* Full CLI with defaults (``--help``)
"""

from __future__ import annotations

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional

# Import necessary helpers
from image_helper import setup_logger, find_measurement_dirs, images_to_dataset, parse_dataset_kwargs, load_config_and_args


# ----------------------------------------------------------------------- #
# Process one measurement
# ----------------------------------------------------------------------- #
def generate_dataloader_for_measurement(
    measurement_dir: Path,
    image_subdir: str,
    image_suffix: str,
    mask_suffix: Optional[str],
    subset_pattern: Optional[str],
    overwrite: bool,
    test: bool,
    logger: logging.Logger,
) -> bool:
    """
    Generate `cp_dataloader.csv` for a single measurement.

    Returns
    -------
    bool
        ``True`` on success.
    """
    # This function is not used in the original `main`, but kept for completeness
    # It would be called inside the main loop if the architecture was different.
    raise NotImplementedError("This function is not fully implemented or used in the provided main.")


# ----------------------------------------------------------------------- #
# CLI
# ----------------------------------------------------------------------- #
def _get_dataloader_parser() -> argparse.ArgumentParser:
    """Defines and returns the ArgumentParser for run_dataloader.py."""
    parser = argparse.ArgumentParser(description="Generate cp_dataloader.csv for CellProfiler in multiple Measurement folders",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", type=str, help="Path to a YAML configuration file to load parameters from.", default=None)
    parser.add_argument("--root_dir", type=str, required=True, help="Root folder containing *_Measurement X folders")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing cp_dataloader.csv")
    parser.add_argument("--test", action="store_true", help="Generate dataloader with ~100 random rows per measurement")
    parser.add_argument("--dataset_kwargs", type=str, default="{}", 
                        help="JSON string of extra kwargs to pass to images_to_dataset(), e.g. "
                        "'{\"image_suffix\": \".tif\", \"mask_suffix\": \".png\", \"image_extractor\": \"...\"}' "
                        "Use null for None, true/false for booleans.")
    return parser


def main() -> None:
    # ------------------------------------------------------------------- #
    # Setup parser and load args
    # ------------------------------------------------------------------- #
    parser = _get_dataloader_parser()
    args = load_config_and_args(parser)
    dataset_kwargs = parse_dataset_kwargs(args.dataset_kwargs)
    
    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    root_dir = Path(args.root_dir)
    if not root_dir.is_dir():
        print(f"Root directory not found: {root_dir}")
        sys.exit(1)

    # Logging
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="dataloader")
    logger.info(f"cp_dataloader generation started - root: {root_dir}")

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

    for mdir in measurements:
        logger.info(f"--- Processing {mdir.name} ---")
        try:
            output_path = mdir / "cp_dataloader.csv"
            if output_path.exists() and not args.overwrite:
                logger.info(f"cp_dataloader.csv already exists → skip {mdir.name}")
                continue

            logger.info(f"Generating dataloader → {output_path.name}")

            # Merge user overrides
            final_kwargs = {"remove_na_row":True, "image_subdir":"Images", "image_suffix":".tiff", "cellprofiler_style":True, **dataset_kwargs}
            dataset = images_to_dataset(mdir, final_kwargs)

            if not dataset or dataset["df"].empty:
                logger.warning(f"No valid images found in {mdir.name}")
                continue

            df = dataset["df"]
            total_rows = len(df)

            if args.test:
                df = df.sample(n=min(100, total_rows), random_state=42)
                logger.info(f"TEST MODE: using {len(df)} random rows (out of {total_rows})")

            # Write CSV
            df.to_csv(output_path, index=False)
            logger.info(f"Wrote {len(df)} rows to {output_path}")

        except Exception as e:
            logger.error(f"Failed to generate dataloader for {mdir.name}: {e}", exc_info=True)
            overall_success = False # Changed to set flag instead of continue

    logger.info("cp_dataloader generation finished")
    sys.exit(0 if overall_success else 1)


if __name__ == "__main__":
    main()