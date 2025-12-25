#!/usr/bin/env python
"""
PerkinElmer (Opera/Phenix/Harmony) Polynomial Flat-Field Correction
- Reads flat-field profiles from Measurement/Image/Index.xml
- Supports multi-channel profiles in a single XML
- Auto-detects all correctable channels if --channels not given
- Writes metadata summary per Measurement
"""

import argparse
import logging
import ast
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import xml.etree.ElementTree as ET

import numpy as np
from tifffile import imread, imwrite
from tqdm import tqdm
from shutil import copy

from image_helper import setup_logger, find_measurement_dirs, images_to_dataset


# --------------------------------------------------------------------------- #
# XML Parsing – one Index.xml contains all channels
# --------------------------------------------------------------------------- #
def parse_index_xml(index_xml_path: Path) -> Dict[str, Dict]:
    """
    Parse Measurement/Image/Index.xml and extract flat-field profile for every channel.
    Returns dict: channel_name (e.g. "ch1", "EGFP") → profile dict
    """
    if not index_xml_path.exists():
        raise FileNotFoundError(f"Index.xml not found: {index_xml_path}")

    tree = ET.parse(index_xml_path)
    root = tree.getroot()

    profiles = {}

    # Find all <Map><Entry ChannelID="...">
    ns = {"pe": "http://www.perkinelmer.com/PEHH/HarmonyV6"}
    for entry in root.findall(".//pe:Map/pe:Entry", ns):
        channel_id = entry.get("ChannelID")
        if not channel_id:
            continue

        ff_node = entry.find("pe:FlatfieldProfile", ns)
        if ff_node is None:
            continue

        # Extract text safely
        def get_text(path, eval=False):
            elem = ff_node.find(path, ns)
            if elem is None or elem.text is None:
                return None
            txt = elem.text.strip()
            return ast.literal_eval(txt) if eval and txt else txt

        try:
            fg_coeffs = ast.literal_eval(ff_node.find("./pe:Foreground/pe:Profile/pe:Coefficients", ns).text)
            bg_coeffs = ast.literal_eval(ff_node.find("./pe:Background/pe:Profile/pe:Coefficients", ns).text)

            dims = ast.literal_eval(ff_node.find("./pe:Foreground/pe:Profile/pe:Dims", ns).text)
            origin = ast.literal_eval(ff_node.find("./pe:Foreground/pe:Profile/pe:Origin", ns).text)
            scale = ast.literal_eval(ff_node.find("./pe:Foreground/pe:Profile/pe:Scale", ns).text)

            bg_mean = float(ff_node.findtext("./pe:Background/pe:Mean", ns) or 0)
            channel_name = ff_node.findtext("pe:ChannelName", ns) or f"ch{channel_id}"

            h, w = dims

            profiles[channel_name] = {
                "channel_id": channel_id,
                "channel_name": channel_name,
                "fg_coeffs": fg_coeffs,
                "bg_coeffs": bg_coeffs,
                "height": h,
                "width": w,
                "origin": origin,
                "scale": scale,
                "bg_mean": bg_mean,
            }
        except Exception as e:
            logging.getLogger().warning(f"Failed to parse flatfield for ChannelID {channel_id}: {e}")

    return profiles


def build_map(coeffs_list: list, h: int, w: int, origin: tuple, scale: tuple) -> np.ndarray:
    oy, ox = origin
    sy, sx = scale
    y, x = np.mgrid[0:h, 0:w]
    x = (x - ox) * sx
    y = (y - oy) * sy

    map_2d = np.zeros((h, w), dtype=np.float64)
    pow_y = 1.0
    for row in coeffs_list:
        pow_x = 1.0
        for c in row:
            map_2d += c * pow_x * pow_y
            pow_x *= x
        pow_y *= y
    return map_2d


def apply_correction_batch(image_paths: List[Path], fg_map: np.ndarray, bg_map: np.ndarray, out_dir: Optional[Path]):
    imgs = np.stack([imread(p) for p in image_paths])
    dtype = imgs.dtype
    corrected = (imgs.astype(np.float64) - bg_map) / fg_map

    info = np.iinfo(dtype)
    corrected = np.clip(corrected, info.min, info.max)
    corrected = np.rint(corrected).astype(dtype)

    for p, img in zip(image_paths, corrected):
        if out_dir:
            rel = p.relative_to(p.parents[2])  # Images → Measurement → root
            dst = out_dir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
        else:
            dst = p
        imwrite(dst, img, compression="zlib")


def write_metadata_txt(measurement_dir: Path, profiles: Dict[str, Dict]):
    txt_path = measurement_dir / "flatfield_metadata.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("PerkinElmer Flat-Field Correction Metadata\n")
        f.write("="*50 + "\n\n")
        for ch, prof in profiles.items():
            f.write(f"Channel        : {ch} (ID: {prof['channel_id']})\n")
            f.write(f"Image Size     : {prof['width']} × {prof['height']}\n")
            f.write(f"Origin         : {prof['origin']}\n")
            f.write(f"Scale          : {prof['scale']}\n")
            f.write(f"Background Mean: {prof['bg_mean']:.5f}\n")
            f.write(f"Foreground Coeffs: {prof['fg_coeffs']}\n")
            f.write(f"Background Coeffs: {prof['bg_coeffs']}\n")
            f.write("-" * 50 + "\n")
    return txt_path


# --------------------------------------------------------------------------- #
# Main correction routine per Measurement
# --------------------------------------------------------------------------- #
def run_correction(measurement_dir: Path, channels_to_correct: Optional[List[str]], output_root: Optional[Path], logger):
    index_xml = measurement_dir / "Image" / "Index.xml"
    if not index_xml.exists():
        logger.warning(f"No Index.xml found in {measurement_dir}")
        return False

    try:
        profiles = parse_index_xml(index_xml)
    except Exception as e:
        logger.error(f"Failed to parse {index_xml}: {e}")
        return False

    if not profiles:
        logger.info(f"No flat-field profiles found in {index_xml}")
        return False

    # Write metadata summary
    metadata_txt = write_metadata_txt(measurement_dir, profiles)
    logger.info(f"Wrote flat-field metadata → {metadata_txt.name}")

    # Determine which channels to process
    available_channels = list(profiles.keys())
    if channels_to_correct is None:
        channels_to_correct = available_channels
        logger.info(f"No --channels specified → correcting all available: {available_channels}")
    else:
        channels_to_correct = [ch for ch in channels_to_correct if ch in available_channels]
        if not channels_to_correct:
            logger.warning("None of the requested channels have flat-field profiles")
            return False

    # Load dataset to know image filenames per channel
    dataset = images_to_dataset(measurement_dir, remove_na_row=True)
    if not dataset:
        logger.error("Failed to build image dataset")
        return False

    df = dataset["df"]
    intensity_cols = dataset["intensity_colnames"]

    success = True

    for ch in channels_to_correct:
        if ch not in intensity_cols:
            logger.warning(f"Channel {ch} not found in image data (skipping)")
            continue

        prof = profiles[ch]
        fg_map = build_map(prof["fg_coeffs"], prof["height"], prof["width"], prof["origin"], prof["scale"])
        bg_map = build_map(prof["bg_coeffs"], prof["height"], prof["width"], prof["origin"], prof["scale"])

        paths = [measurement_dir / "Images" / f for f in df[ch].dropna()]
        logger.info(f"Correcting {len(paths)} images → {ch} ({prof['width']}×{prof['height']})")

        for i in tqdm(range(0, len(paths), 50), desc=ch, leave=False):
            apply_correction_batch(paths[i:i+50], fg_map, bg_map, output_root)

    # Copy uncorrected channels (if output_root used)
    if output_root:
        corrected_set = set(channels_to_correct)
        uncorrected = [c for c in intensity_cols if c not in corrected_set]
        target_img_dir = output_root / measurement_dir.name / "Images"
        target_img_dir.mkdir(parents=True, exist_ok=True)
        for ch in uncorrected:
            for fname in df[ch].dropna():
                src = measurement_dir / "Images" / fname
                dst = target_img_dir / fname
                dst.parent.mkdir(parents=True, exist_ok=True)
                if src.exists():
                    copy(src, dst)

    return success


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main():
    parser = argparse.ArgumentParser(description="PerkinElmer Polynomial Flat-Field Correction from Index.xml",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("root", type=Path, help="Root folder containing Measurement directories")
    parser.add_argument("--channels", nargs="*", default=None, help="Channels to correct (e.g. ch1 EGFP). If omitted → all available in Index.xml")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output root (optional). If omitted → overwrite in place")
    args = parser.parse_args()
    
    # ------------------------------------------------------------------- #
    # Setup
    # ------------------------------------------------------------------- #
    logger = setup_logger(args.root / "logs", "polynomial_correction")
    logger.info("PerkinElmer Polynomial Correction Started")

    measurements = find_measurement_dirs(args.root)
    if not measurements: 
        logger.error("No Measurement directories found")
        return

    logger.info(f"Found {len(measurements)} measurement(s)")

    for mdir in measurements:
        logger.info(f"Processing → {mdir.name}")
        try:
            run_correction(mdir, args.channels, args.output_dir, logger)
        except Exception as e:
            logger.error(f"Failed on {mdir.name}: {e}", exc_info=True)

    logger.info("All done.")


if __name__ == "__main__":
    main()