#!/usr/bin/env python
"""
BaSiC Shading Correction Pipeline
Supports:
  • fit   – train BaSiC model per channel (saved as .pkl + profiles)
  • transform – apply previously trained models
Works on multiple Measurement directories.
"""

import argparse
import logging
import pickle
import random
import sys
from pathlib import Path
from typing import List, Optional, Dict

import numpy as np
from basicpy import BaSiC
from tifffile import imread, imwrite
from tqdm import tqdm
from shutil import copy

# Import necessary helpers
from image_helper import setup_logger, find_measurement_dirs, images_to_dataset, parse_dataset_kwargs, load_config_and_args


# ==============================
# BaSiC Helpers
# ==============================
def _basic_fit(
    image_paths: List[Path],
    n_image: int = 50,
    enable_darkfield: bool = True,
    working_size: int = 128,
    logger: Optional[logging.Logger] = None,
) -> Optional[BaSiC]:
    if logger is None:
        logger = logging.getLogger("basic")

    if len(image_paths) > n_image:
        image_paths = random.sample(image_paths, k=n_image)

    try:
        imgs = np.stack([imread(p) for p in image_paths])
        # ... (rest of function unchanged)
        basic = BaSiC(
            get_darkfield=enable_darkfield,
            smoothness_flatfield=1,
            smoothness_darkfield=1,
            working_size=working_size,
            max_workers=8,
        )
        basic.fit(imgs)
        logger.info(f"BaSiC model fitted on {len(imgs)} images.")
        return basic
    except Exception as e:
        logger.error(f"BaSiC fitting failed: {e}")
        return None


def _basic_transform(
    image_paths: List[Path],
    model: BaSiC,
    target_dir: Optional[Path] = None,
    timelapse: bool = False,
    logger: Optional[logging.Logger] = None,
) -> bool:
    if logger is None:
        logger = logging.getLogger("basic")

    try:
        imgs = np.stack([imread(p) for p in image_paths])
        dtype_in = imgs.dtype
        corrected = model.transform(imgs, timelapse=timelapse)

        # Clip to original bit depth
        if dtype_in == np.uint16:
            corrected = np.clip(corrected, 0, 65535)
        elif dtype_in == np.uint8:
            corrected = np.clip(corrected, 0, 255)

        corrected = corrected.astype(dtype_in)

        for src_path, corr_img in zip(image_paths, corrected):
            if target_dir is not None:
                rel = src_path.relative_to(src_path.parents[2])  # assumes /Measurement/Images/...
                dst_path = target_dir / rel
                dst_path.parent.mkdir(parents=True, exist_ok=True)
                imwrite(dst_path, corr_img, compression="zlib")
            else:
                imwrite(src_path, corr_img, compression="zlib")
        return True
    except Exception as e:
        logger.error(f"BaSiC transform failed: {e}")
        return False


# ==============================
# Fit mode
# ==============================
def fit_basic_models(
    dataset: Dict,
    channels: List[str],
    n_image: int,
    working_size: int,
    darkfield: bool,
    output_root: Optional[Path],
    logger: logging.Logger,
) -> bool:
    df = dataset["df"]
    intensity_cols = dataset["intensity_colnames"]
    images_dir = Path(df['directory'].iloc[0])
    measurement_dir = images_dir.parent

    valid_channels = [c for c in channels if c in intensity_cols]
    if not valid_channels:
        logger.warning("No requested channels found in this measurement.")
        return False

    basic_dir = (output_root / measurement_dir.name / "BaSiC_model") if output_root else (measurement_dir / "BaSiC_model")
    basic_dir.mkdir(parents=True, exist_ok=True)

    success = True
    for chan in valid_channels:
        logger.info(f"Fitting BaSiC model for channel {chan} in {measurement_dir.name}")
        paths = [images_dir / f for f in df[chan].dropna()]
        model = _basic_fit(paths, n_image, darkfield, working_size, logger)
        if model is None:
            success = False
            continue

        # Save model + profiles
        with open(basic_dir / f"model_{chan}.pkl", "wb") as f:
            pickle.dump(model, f)
        imwrite(basic_dir / f"model_{chan}_flatfield.tiff", model.flatfield.astype(np.float32), compression="zlib")
        if darkfield:
            imwrite(basic_dir / f"model_{chan}_darkfield.tiff", model.darkfield.astype(np.float32), compression="zlib")

        # Test on 2 random images
        test_paths = random.sample(paths, min(2, len(paths)))
        corrected = model.transform(np.stack([imread(p) for p in test_paths]))
        for i, p in enumerate(test_paths):
            copy(p, basic_dir / p.name)
            # Check suffix for correct dtype, assuming .tiff is uint16
            dtype = np.uint16 if p.suffix == ".tiff" else np.uint8
            imwrite(basic_dir / f"{p.stem}_corrected.tiff", corrected[i].astype(dtype), compression="zlib")

    return success


# ==============================
# Transform mode
# ==============================
def apply_basic_correction(
    dataset: Dict,
    channels: List[str],
    output_root: Optional[Path],
    logger: logging.Logger,
) -> bool:
    df = dataset["df"]
    intensity_cols = dataset["intensity_colnames"]
    valid_channels = [c for c in channels if c in intensity_cols]
    images_dir = Path(df['directory'].iloc[0])
    measurement_dir = images_dir.parent

    model_base = (output_root / measurement_dir.name) if output_root else measurement_dir
    success = True

    for chan in valid_channels:
        model_path = model_base / "BaSiC_model" / f"model_{chan}.pkl"
        if not model_path.exists():
            logger.error(f"Model missing for {chan}: {model_path}")
            success = False
            continue

        with open(model_path, "rb") as f:
            model = pickle.load(f)

        paths = [images_dir / f for f in df[chan].dropna()]
        logger.info(f"Correcting {len(paths)} images – {chan} – {measurement_dir.name}")

        for batch_start in tqdm(range(0, len(paths), 50), desc=f"Correct {chan}"):
            batch = paths[batch_start:batch_start + 50]
            # output_root is passed as the target_dir
            if not _basic_transform(batch, model, output_root, logger=logger):
                success = False

    # Copy unprocessed channels
    if output_root:
        unprocessed = set(intensity_cols) - set(valid_channels)
        target_img_dir = output_root / measurement_dir.name / "Images"
        target_img_dir.mkdir(parents=True, exist_ok=True)
        for ch in unprocessed:
            logger.info(f"Copying untouched channel {ch}")
            for fname in df[ch].dropna():
                src = images_dir / fname
                dst = target_img_dir / fname
                dst.parent.mkdir(parents=True, exist_ok=True)
                if src.exists():
                    copy(src, dst)

    return success


# ==============================
# CLI
# ==============================
def _get_basic_parser() -> argparse.ArgumentParser:
    """Defines and returns the ArgumentParser for run_basic_correction.py."""
    parser = argparse.ArgumentParser(description="BaSiC Shading Correction (fit or transform)",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", type=str, help="Path to a YAML configuration file to load parameters from.", default=None)
    parser.add_argument("--root_dir", type=str, required=True, help="Root folder containing Measurement directories")
    parser.add_argument("--mode", choices=["fit", "transform"], required=True)
    parser.add_argument("--channels", nargs="+", default=["ch1"], help="Channels to process, e.g. ch1 ch2")
    parser.add_argument("--n_image", type=int, default=50, help="Number of images used for fitting")
    parser.add_argument("--working_size", type=int, default=128, help="BaSiC working size (resized during fit)")
    parser.add_argument("--darkfield", action="store_true", help="Enable darkfield estimation (default: True)")
    parser.add_argument("--output_dir", type=str, default=None, help="Optional output root (for transform/copy)")
    parser.add_argument("--dataset_kwargs", type=str, default="{}", 
                        help="JSON string of extra kwargs to pass to images_to_dataset(), e.g. "
                         "'{\"image_suffix\": \".tif\", \"mask_suffix\": \".png\", \"image_extractor\": \"...\"}' "
                         "Use null for None, true/false for booleans.")
    return parser


def main():
    parser = _get_basic_parser()
    args = load_config_and_args(parser)
    dataset_kwargs = parse_dataset_kwargs(args.dataset_kwargs)
    
    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    root_dir = Path(args.root_dir)
    output_dir = Path(args.output_dir) if args.output_dir else None
    
    # Initial checks for root_dir
    if not root_dir.is_dir():
        print(f"Root directory not found: {root_dir}")
        sys.exit(1)
        
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="basic_correction")
    logger.info(f"BaSiC {args.mode.upper()} mode started")

    measurements = find_measurement_dirs(root_dir)
    if not measurements:
        logger.error("No Measurement folders found.")
        sys.exit(1)

    logger.info(f"Found {len(measurements)} Measurement directories")
    overall_success = True

    for mdir in measurements:
        logger.info(f"{'Fitting' if args.mode == 'fit' else 'Transforming'} → {mdir.name}")
        try:
            # Merge user overrides
            final_kwargs = {"remove_na_row":True, "image_subdir":"Images", "image_suffix":".tiff", **dataset_kwargs}
            dataset = images_to_dataset(mdir, final_kwargs)
            if not dataset:
                logger.error(f"Failed to build dataset for {mdir.name}")
                overall_success = False
                continue
        
            if args.mode == "fit":
                ok = fit_basic_models(
                    dataset, args.channels, args.n_image, args.working_size, args.darkfield, output_dir, logger)
                if not ok:
                    overall_success = False
            else:  # transform
                ok = apply_basic_correction(dataset, args.channels, output_dir, logger)
                if not ok:
                    overall_success = False
                
        except Exception as e:
            logger.error(f"Failed on {mdir.name}: {e}", exc_info=True)
            overall_success = False

    logger.info("BaSiC processing completed.")
    sys.exit(0 if overall_success else 1)


if __name__ == "__main__":
    main()