#!/usr/bin/env python
"""
run_cellpose.py
---------------
Cellpose-SAM segmentation for multiple Measurement folders.

* Cellpose-SAM uses the **first three channels** of the input image.
* The script builds those three channels according to the user-defined groups.
* Supports:
    1. single channel → C1 = image, C2 = 0
    2. two channels  → C1 = chA, C2 = chB
    3. two merged channels → C1 = merge(chan1), C2 = merge(chan2)
* All arguments have sensible defaults and are displayed with ``--help``.
"""

from __future__ import annotations

# import gc
import re, sys, argparse, logging
logging.getLogger('cellpose').setLevel(logging.WARNING)
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np
import pandas as pd
import torch
from cellpose import models
from cellpose import io as cp_io
from skimage.morphology import closing
from skimage.transform import rescale
from tifffile import imread
from tqdm import tqdm

from image_helper import setup_logger, find_measurement_dirs, images_to_dataset, parse_dataset_kwargs, load_config_and_args



# ----------------------------------------------------------------------- #
# Image loading & channel construction
# ----------------------------------------------------------------------- #
def _merge_channels(
    paths: List[Path],
    method: str = "mean",
    resize_factor: float = 1.0,
) -> np.ndarray:
    """Read, (optionally) resize and merge a list of images."""
    imgs = [imread(p) for p in paths]
    stacked = np.stack(imgs, axis=0)

    if method == "mean":
        merged = np.mean(stacked, axis=0)
    elif method == "max":
        merged = np.max(stacked, axis=0)
    elif method == "min":
        merged = np.min(stacked, axis=0)
    else:
        raise ValueError(f"Unsupported merge method: {method}")

    if resize_factor != 1.0:
        merged = rescale(
            merged,
            resize_factor,
            anti_aliasing=True,
            preserve_range=True,
        ).astype(stacked.dtype)

    return merged


def _build_cellpose_image(
    df: pd.DataFrame,
    row_idx: int,
    chan1: List[str],
    chan2: Optional[List[str]],
    merge1: str,
    merge2: str,
    resize_factor: float,
) -> np.ndarray:
    """
    Build a **(2, H, W)** image for Cellpose‑SAM.

    * C1 = merge(chan1)  (mean / max / min)
    * C2 = merge(chan2)  (if chan2 supplied, else 0)
    """
    row = df.iloc[row_idx]
    dir_path = Path(row["directory"])

    # ---- group 1 (always required) ----
    ch1_paths = [dir_path / row[ch] for ch in chan1 if ch in row and pd.notna(row[ch])]
    if not ch1_paths:
        raise ValueError(f"Missing images for channel group 1: {chan1}")
    c1 = _merge_channels(ch1_paths, merge1, resize_factor)

    # ---- group 2 (optional) ----
    if chan2:
        ch2_paths = [dir_path / row[ch] for ch in chan2 if ch in row and pd.notna(row[ch])]
        if not ch2_paths:
            raise ValueError(f"Missing images for channel group 2: {chan2}")
        c2 = _merge_channels(ch2_paths, merge2, resize_factor)
    else:
        c2 = np.zeros_like(c1)

    # Stack → (2, H, W)
    return np.stack([c1, c2], axis=0)



# ----------------------------------------------------------------------- #
# Segmentation per measurement
# ----------------------------------------------------------------------- #
def cellpose_segment_measurement(
    dataset: Dict,
    chan1: List[str],
    chan2: Optional[List[str]],
    merge1: str,
    merge2: str,
    model_name: str,
    diameter: Optional[float],
    normalize: Optional[Dict],
    resize_factor: float,
    mask_name: str,
    reduce_mask_name: bool,
    overwrite_mask: bool,
    flow_threshold: float,
    cellprob_threshold: float,
    device: torch.device,
    gpu_batch_size: int,
    logger: logging.Logger,
) -> bool:
    """Run Cellpose-SAM on a single measurement folder."""

    df = dataset["df"]
    intensity_cols = dataset["intensity_colnames"]

    # sanity‑check requested channels
    missing = [ch for ch in (chan1 + (chan2 or [])) if ch not in intensity_cols]
    if missing:
        logger.error(f"Requested channels not present: {set(missing)}")
        return False

    # load model once per measurement
    try:
        model = models.CellposeModel(device=device, pretrained_model=model_name)
    except Exception as e:
        logger.error(f"Cannot load model '{model_name}': {e}")
        return False

    diameter_val = None if diameter is None or diameter <= 0 else int(diameter * resize_factor)

    logger.info(
        f"Config → model:{model_name} | diameter:{diameter_val} | "
        f"chan1:{chan1}({merge1}) | chan2:{chan2 or 'none'}({merge2 or '-'}) | "
        f"resize:{resize_factor:.2f} | mask:{mask_name}"
    )

    success = True
    for idx in tqdm(range(len(df)), desc=f"Cellpose", unit="img"):
        row = df.iloc[idx]

        # use first channel of chan1 to derive output stem
        stem_ch = chan1[0]
        src_path = Path(row["directory"]) / row[stem_ch]
        if not src_path.exists():
            logger.warning(f"Source image missing: {src_path}")
            success = False
            continue

        # mask filename
        save_stem = src_path.parent / f"{src_path.stem}_cp_masks"
        if reduce_mask_name:
            save_stem = Path(re.sub(r"_ch\d+_", "_ch0_", str(save_stem)))
        mask_path = save_stem.with_name(f"{save_stem.name}_{mask_name}.png")

        if mask_path.exists() and not overwrite_mask:
            continue

        try:
            img = _build_cellpose_image(
                df, idx,
                chan1, chan2, merge1, merge2,
                resize_factor,
                )

            masks, flows, styles = model.eval(
                img,
                batch_size=gpu_batch_size,
                # channels=[0, 1, 2],   # Cellpose‑SAM uses first three channels
                normalize=normalize,
                diameter=diameter_val,
                flow_threshold=flow_threshold,
                cellprob_threshold=cellprob_threshold,
                )

            # rescale mask back if we down‑sampled the input
            if resize_factor != 1.0:
                masks = rescale(masks, 1.0 / resize_factor, order=0).astype(np.uint16)

            if len(np.unique(masks)) <= 1:
                logger.debug(f"No objects detected in {src_path.name}")
                continue

            cp_io.save_masks(
                img[0], # original merged C1 (for overlay)
                closing(masks),
                flows,
                file_names=str(save_stem),
                suffix=f"_{mask_name}",
            )

        except torch.cuda.OutOfMemoryError:
            logger.error(f"GPU OOM on {src_path.name} – skipping")
            torch.cuda.empty_cache()
            success = False
        except Exception as e:
            logger.error(f"Segmentation error on {src_path.name}: {e}", exc_info=True)
            success = False

        # periodic cleanup
        if idx % 200 == 0:
            torch.cuda.empty_cache()
            # gc.collect()

    return success



# ----------------------------------------------------------------------- #
# CLI & argument handling
# ----------------------------------------------------------------------- #
def _parse_normalize(s: str) -> Optional[Dict]:
    """Parse ``--normalize percentile:0.5,99.5`` → dict (or True/None)."""
    s = s.strip().lower()
    if s in ("true", "1"):
        return True
    if s in ("false", "0", "none"):
        return None
    if ":" not in s:
        raise argparse.ArgumentTypeError("normalize must be 'percentile:low,high'")
    kind, vals = s.split(":", 1)
    low, high = map(float, vals.split(","))
    return {kind: [low, high]}


def _get_cellpose_parser() -> argparse.ArgumentParser:
    """Defines and returns the ArgumentParser for run_cellpose.py."""
    parser = argparse.ArgumentParser(description="Cellpose-SAM segmentation for multiple *-Measurement folders",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", type=str, help="Path to a YAML configuration file to load parameters from.", default=None)
    
    # Channel groups (group‑1 is mandatory, group‑2 optional)
    parser.add_argument("--root_dir", type=str, required=True, help="Root folder containing *_Measurement X sub‑folders")
    parser.add_argument("--chan1", nargs="+", type=str, required=True, help="First channel group (C1). At least one channel required.")
    parser.add_argument("--chan2", nargs="*", type=str, default=[], help="Optional second channel group (C2). If omitted → C2 = 0.")
    parser.add_argument("--merge1", choices=["mean", "max", "min"], default="mean", help="How to merge channels **inside** chan1 (default: mean)")
    parser.add_argument("--merge2", choices=["mean", "max", "min"], default="mean", help="How to merge channels **inside** chan2 (default: mean)")
    # Cellpose‑SAM parameters
    parser.add_argument("--model", default="cpsam", help="Cellpose-SAM model name or local path")
    parser.add_argument("--diameter", type=float, default=0, help="Approximate object diameter in pixels (0 → auto‑estimate)")
    parser.add_argument("--normalize", type=_parse_normalize, default="percentile:0.1,99.9", help="Normalization: 'percentile:low,high', True, or None")
    parser.add_argument("--resize", type=float, default=1.0, help="Resize factor applied **before** segmentation (1.0 = original size)")
    parser.add_argument("--mask_name", default="cell", help="Suffix added to mask filenames")
    parser.add_argument("--no_reduce_mask_name", action="store_true", help="Do not replace '_chN_' with '_ch0_' in mask filenames")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing mask files")
    parser.add_argument("--flow_thresh", type=float, default=0.4, help="Flow error threshold")
    parser.add_argument("--cellprob_thresh", type=float, default=0.0, help="Cell probability threshold")
    parser.add_argument("--device", choices=["auto", "gpu", "cpu"], default="auto", help="Computation device")
    parser.add_argument("--batch_size", type=int, default=16, help="GPU batch size (larger → faster but more VRAM)")
    # Note: Using 'dataset_kwargs' with underscores for command line compatibility.
    parser.add_argument("--dataset_kwargs", type=str, default="{}", 
                        help="JSON string of extra kwargs to pass to images_to_dataset(), e.g. "
                         "'{\"image_suffix\": \".tif\", \"mask_suffix\": \".png\", \"image_extractor\": \"...\"}' "
                         "Use null for None, true/false for booleans.")
    return parser


def main() -> None:
    # ------------------------------------------------------------------- #
    # Setup parser
    # ------------------------------------------------------------------- #
    parser = _get_cellpose_parser()
    args = load_config_and_args(parser)
    dataset_kwargs = parse_dataset_kwargs(args.dataset_kwargs)

    # ------------------------------------------------------------------- #
    # Setup measurement
    # ------------------------------------------------------------------- #
    root_dir = Path(args.root_dir)
    if not root_dir.is_dir():
        print(f"Root directory does not exist: {root_dir}")
        sys.exit(1)

    # device
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    elif args.device == "gpu":
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    # logging
    log_dir = root_dir / "logs"
    logger = setup_logger(log_dir, name="cellpose")
    logger.info(f"Cellpose-SAM started - root: {root_dir}")
    logger.info(
        f"Device: {device} | Model: {args.model} | Resize: {args.resize} | Normalization: {args.normalize}"
        f"Channel 1: {args.chan1}({args.merge1}) | Channel 2: {args.chan2 or 'none'}({args.merge2})")

    # discover measurements
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
            # Merge user overrides
            final_kwargs = {"remove_na_row":True, "image_subdir":"Images", "image_suffix":".tiff", **dataset_kwargs}
            dataset = images_to_dataset(mdir, final_kwargs)
            if not dataset:
                logger.error(f"Failed to build dataset for {mdir}")
                return False

            ok = cellpose_segment_measurement(
                dataset=dataset,
                chan1=args.chan1,
                chan2=args.chan2 or None,
                merge1=args.merge1,
                merge2=args.merge2,
                model_name=args.model,
                diameter=args.diameter,
                normalize=args.normalize,
                resize_factor=args.resize,
                mask_name=args.mask_name,
                reduce_mask_name=not args.no_reduce_mask_name,
                overwrite_mask=args.overwrite,
                flow_threshold=args.flow_thresh,
                cellprob_threshold=args.cellprob_thresh,
                device=device,
                gpu_batch_size=args.batch_size,
                logger=logger,
            )
            if not ok:
                overall_success = False
        except Exception as e:
            logger.error(f"Unexpected error in {mdir.name}: {e}", exc_info=True)
            overall_success = False

    logger.info("Cellpose-SAM pipeline finished")
    sys.exit(0 if overall_success else 1)


if __name__ == "__main__":
    main()