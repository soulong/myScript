#!/usr/bin/env python
"""
image_helper.py
---------------
Shared utilities for the preprocessing pipeline:

* ``setup_logger`` - console + file logging
* ``find_measurement_dirs`` - discover ``*_Measurement X`` folders
* ``images_to_dataset`` - build the per-position DataFrame (unchanged logic,
  only Path-based, type-annotated, and robust)

All functions use **Path** objects only.
"""

from __future__ import annotations

import re, string, sys, logging, argparse
import yaml, json
from typing import Union, Iterable, Dict, List, Optional, Any
from pathlib import Path
from glob import glob
import pandas as pd
import numpy as np
from natsort import natsorted
import imageio.v3 as iio
import scipy.ndimage as ndi
from skimage import transform
from skimage.measure import regionprops, label


# --------------------------------------------------------------------------- #
# Logger
# --------------------------------------------------------------------------- #
def setup_logger(log_dir: Path, name="preprocessing") -> logging.Logger:
    """
    Create (or reuse) a logger that writes to console **and** to
    ``<log_dir>/preprocessing.log``.

    Parameters
    ----------
    log_dir: Path
        Directory where the log file will be created.

    Returns
    -------
    logging.Logger
    """
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(name)
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s")

    # console
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # file
    fh = logging.FileHandler(log_dir / f"{name}.log", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# --------------------------------------------------------------------------- #
# Measurement discovery
# --------------------------------------------------------------------------- #
def find_measurement_dirs(
    root_dir: Path | str,
    measurement_pattern: str = r".*-Measurement.*$",
    image_sub_directory: str = "Images"
) -> List[Path]:
    """
    Find all measurement directories (i.e. folders like 'Exp_Measurement 1')
    that contain an 'Images' subfolder.

    Returns sorted list of Path objects pointing to the measurement folders.
    """
    root = Path(root_dir)
    
    if not root.is_dir():
        raise NotADirectoryError(f"{root} is not a valid directory")
    
    if not root.exists():
        raise FileNotFoundError(f"{root} is not existed")

    # Compile regex once
    pattern = re.compile(measurement_pattern) # re.IGNORECASE

    # Directly search for the 'Images' folder - this is efficient
    images_dirs = root.rglob(f"*/{image_sub_directory}")

    measurement_dirs = []
    for img_dir in images_dirs:
        parent = img_dir.parent
        # Check if parent name matches the measurement pattern
        if pattern.match(parent.name):
            measurement_dirs.append(parent)

    # Natural sort by full path (or just name if preferred)
    return natsorted(measurement_dirs, key=lambda p: p.name)



# --------------------------------------------------------------------------- #
# Build per-position DataFrame (unchanged core logic, Path-based)
# --------------------------------------------------------------------------- #
def images_to_dataset(
    measurement_dir: Path | str,
    subset_pattern: Optional[str] = None,
    image_subdir: str = "Images",
    image_suffix: str = ".tiff",
    mask_suffix: Optional[str] = ".png",
    remove_na_row: bool = True,
    cellprofiler_style: bool = False,
    position_metadata: Optional[Path] = None,
    pos_index_col: str = "P Index",
    image_extractor: str = (
        r"r(?P<row>.*)c(?P<column>.*)f(?P<field>.*)p(?P<stack>.*)-ch(?P<channel>[0-9]{1,})"
        r"sk(?P<timepoint>[0-9]{1,})fk1fl1(?P<self_generated>.*)"),
    mask_extractor: str = (
        r"r(?P<row>.*)c(?P<column>.*)f(?P<field>.*)p(?P<stack>.*)-ch(?P<channel>[0-9]{1,})"
        r"sk(?P<timepoint>[0-9]{1,})fk1fl1_(cp_masks_)?(?P<mask_name>.*)"),
) -> Optional[Dict]:
    """
    Convert a folder of images (and optional masks) into a tidy DataFrame.

    Parameters
    ----------
    measurement_dir: Path
        Root of the measurement (the ``*-Measurement X`` folder).
    subset_pattern, image_subdir, image_suffix, mask_suffix, remove_na_row,
    cellprofiler_style, position_metadata, pos_index_col, image_extractor,
    mask_extractor
        See original docstring - behaviour unchanged.

    Returns
    -------
    dict | None
        ``{'df': pd.DataFrame,
           'metadata_colnames': list,
           'intensity_colnames': list,
           'mask_colnames': list | None}``
        ``None`` if the directory does not exist or no images are found.
    """
    # ------------------------------------------------------------------- #
    # Basic validation
    # ------------------------------------------------------------------- #
    if isinstance(measurement_dir, str): measurement_dir=Path(measurement_dir)
    if not measurement_dir.is_dir():
        return None
    
    # print(f"\n>>>>> processing: {measurement_dir}")

    image_extract_cols = re.findall(r"\?P<([^>]+)>", image_extractor)
    if not ("channel" in image_extract_cols or "self_generated" in image_extract_cols):
        raise ValueError("image_extractor must contain [self_generated, channel]")

    # ------------------------------------------------------------------- #
    # Intensity images
    # ------------------------------------------------------------------- #
    image_glob = str(measurement_dir / image_subdir / f"**/*{image_suffix}")
    image_files = glob(image_glob, recursive=True)
    image_files = [Path(p) for p in image_files]

    if subset_pattern:
        image_files = [p for p in image_files if re.search(subset_pattern, str(p))]

    if not image_files:
        print("no valid images found, check directory or subset_pattern")
        return None

    # Build DataFrame
    image_df = pd.DataFrame({
        "directory": [p.parent for p in image_files],
        "filename": [p.name for p in image_files],
    })

    image_meta = image_df["filename"].str.extract(image_extractor + re.escape(image_suffix))

    if "channel" in image_meta.columns:
        image_meta["channel"] = image_meta["channel"].astype(str) + image_meta["self_generated"]
    else:
        image_meta["channel"] = "0"
    image_meta.pop("self_generated")

    image_df = pd.concat([image_df, image_meta], axis=1)

    # Index columns (everything except channel)
    image_index_cols = ["directory"] + [
        c for c in image_extract_cols if c not in ("channel", "self_generated")]
    # set index for merging with mask
    image_df = image_df.set_index(image_index_cols)
    image_df["channel"] = "ch" + image_df["channel"]
    image_colnames = natsorted(image_df["channel"].unique().tolist())

    # Pivot to wide format
    image_df = image_df.pivot(columns="channel", values="filename")
    print(f"{len(image_df)} grouped intensity images: {image_colnames}")
    # print(image_df)

    # ------------------------------------------------------------------- #
    # Masks (optional)
    # ------------------------------------------------------------------- #
    if mask_suffix is None:
        mask_colnames = None
        df = image_df
    else:
        mask_glob = str(measurement_dir / image_subdir / f"**/*{mask_suffix}")
        mask_files = glob(mask_glob, recursive=True)
        mask_files = [Path(p) for p in mask_files]

        if not mask_files:
            mask_colnames = None
            df = image_df
        else:
            mask_df = pd.DataFrame({
                "directory": [p.parent for p in mask_files],
                "filename": [p.name for p in mask_files],
            })
            mask_meta = mask_df["filename"].str.extract(mask_extractor + re.escape(mask_suffix))

            if "channel" not in mask_meta.columns:
                mask_meta["channel"] = "0"

            mask_df = pd.concat([mask_df, mask_meta], axis=1)
            mask_df = mask_df.set_index(image_index_cols)

            mask_colnames = natsorted(mask_df["mask_name"].unique().tolist())
            mask_df = mask_df.pivot(columns="mask_name", values="filename")
            print(f"{len(mask_df)} grouped object masks: {mask_colnames}")

            df = pd.merge(image_df, mask_df, how="left", left_index=True, right_index=True)
    
    # to int if possible
    # print(df)
    # # print(df["row"].dtype)
    df = df.applymap(lambda x : pd.to_numeric(x, errors='ignore'))
    # # print(df["row"].dtype)
    # # print(df["directory"].dtype)

    # ------------------------------------------------------------------- #
    # Clean-up
    # ------------------------------------------------------------------- #
    if remove_na_row:
        df = df.dropna()
    df = df.reset_index()
    print(f"{len(df)} groups after merging images and masks")

    metadata_colnames = image_index_cols.copy()
    # print(df)

    # ------------------------------------------------------------------- #
    # merge row and column to well
    # ------------------------------------------------------------------- #
    if "well" not in df.columns and {"row", "column"}.issubset(df.columns):
        row_map = {i: letter for i, letter in enumerate(string.ascii_uppercase, 1)}
        # Only apply mapping if row values look like numbers (1,2,3…) and NOT letters (A,B,C…)
        row_sample = df["row"].dropna().astype(str).str.strip().str.upper()
        if row_sample.str.match(r"^\d+$").all():
            df["row"]=pd.to_numeric(df["row"], errors="coerce").map(row_map).fillna("")
        else:
            df["row"]=df["row"].astype(str)
        df["column"]=df["column"].astype("Int64").astype(str)
        # Insert well right after directory
        df.insert(loc=df.columns.get_loc("directory") + 1, column='well', 
                  value=df["row"] + df["column"])
        # Drop helper & original columns
        df = df.drop(columns=["row", "column"])
        # Update metadata list
        metadata_colnames = [c for c in metadata_colnames if c not in ("row", "column")]
        metadata_colnames.insert(metadata_colnames.index("directory") + 1, "well")
        # print(df)
        # print(metadata_colnames)

    # ------------------------------------------------------------------- #
    # Optional position metadata (Nikon ND2 export)
    # ------------------------------------------------------------------- #
    if position_metadata and (measurement_dir / position_metadata).exists():
        pos_meta = pd.read_csv(measurement_dir / position_metadata)
        pos_meta = pos_meta.dropna(subset=[pos_index_col])
        pos_meta = pos_meta.drop_duplicates(
            subset=["prefix", "Time [s]", pos_index_col], keep="first"
        )

        if not pd.isna(pos_meta.iloc[0].get("Position Name")):
            print("add position well metadata")
            well_meta = pos_meta["Position Name"].str.extract(
                r"(?P<row>[A-Z])(?P<column>[0-9]+)#(?P<field>.+)"
            )
            well_meta["well"] = well_meta["row"] + well_meta["column"]

            if "T Index" in pos_meta.columns:
                pos_meta["T Index"] = pos_meta["T Index"].fillna(1)
                pos_meta = pd.concat(
                    [
                        pos_meta["prefix"],
                        pos_meta[pos_index_col].rename("position"),
                        pos_meta["T Index"].rename("timepoint"),
                        well_meta,
                    ],
                    axis=1,
                )
                pos_meta = pos_meta.drop(["row", "column"], axis=1)
                pos_meta["prefix"] = pos_meta["prefix"].astype(str)
                df = pd.merge(
                    pos_meta, df, how="right", on=["prefix", "position", "timepoint"]
                )
            else:
                pos_meta = pd.concat(
                    [pos_meta["prefix"], pos_meta[pos_index_col].rename("position"), well_meta],
                    axis=1,
                )
                pos_meta = pos_meta.drop(["row", "column"], axis=1)
                pos_meta["prefix"] = pos_meta["prefix"].astype(str)
                df = pd.merge(pos_meta, df, how="right", on=["prefix", "position"])

            metadata_colnames.extend(["well", "field"])

    # ------------------------------------------------------------------- #
    # CellProfiler style (optional)
    # ------------------------------------------------------------------- #
    if cellprofiler_style:
        for ch in image_colnames:
            df[f"Image_PathName_{ch}"] = df["directory"]
            df = df.rename(columns={ch: f"Image_FileName_{ch}"})
        if mask_colnames:
            for m in mask_colnames:
                df[f"Image_ObjectsPathName_mask_{m}"] = df["directory"]
                df = df.rename(columns={m: f"Image_ObjectsFileName_mask_{m}"})
        for meta in metadata_colnames:
            df = df.rename(columns={meta: f"Metadata_{meta}"})

    return {
        "df": df,
        "metadata_colnames": metadata_colnames,
        "intensity_colnames": image_colnames,
        "mask_colnames": mask_colnames,
    }




# --------------------------------------------------------------------------- #
# Command line helper
# --------------------------------------------------------------------------- #
def parse_dataset_kwargs(s: str) -> dict:
    s = s.strip()
    if not s or s == "{}":
        return {}
    try:
        raw = json.loads(s)
        # Convert "null" → None, "true" → True, etc.
        return {
            k: (None if v is None else v)
            for k, v in raw.items()
        }
    except json.JSONDecodeError as e:
        raise argparse.ArgumentTypeError(f"Invalid JSON in --dataset-kwargs: {e}")
    

def load_config_and_args(
    parser: argparse.ArgumentParser
) -> argparse.Namespace:
    """
    Loads configuration from a YAML file, merges with command-line arguments,
    and returns a fully populated argparse.Namespace object.

    Command-line arguments always override YAML file settings.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        An ArgumentParser object with all required arguments defined.

    Returns
    -------
    argparse.Namespace
        The final Namespace object with merged and typed parameters.
    """
    # 1. Initial Parse: Get the config path
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", type=str, help="Path to a YAML configuration file.")
    
    # parse_known_args is used to prevent the main parser from failing if 
    # required arguments (like 'root_dir') are missing when only '--config' is passed.
    config_args, _ = config_parser.parse_known_args()

    # 2. Load parameters from YAML file
    yaml_config: Dict[str, Any] = {}
    if config_args.config:
        config_path = Path(config_args.config)
        if not config_path.exists():
            print(f"Config file not found: {config_path}")
            sys.exit(1)
        try:
            with open(config_path, 'r') as f:
                # Use safe_load to load the dictionary
                yaml_config = yaml.safe_load(f) or {}
            # print(f"Loaded config from {config_path}")
        except Exception as e:
            print(f"Error loading YAML config from {config_path}: {e}")
            sys.exit(1)

    # 3. Final Parse: Get command-line arguments
    # We must parse all arguments again to get the command-line values.
    # We use a temporary parser for this without any defaults set yet.
    cmd_args = parser.parse_args()
    cmd_args_dict = vars(cmd_args)

    # 4a. Get the YAML dict (if present). Remove from config to avoid merge confusion.
    yaml_dataset_kwargs = yaml_config.pop('dataset_kwargs', None)
    
    # 4b. Get the CLI string (which is always present in cmd_args_dict).
    # Remove from CLI args to avoid merge confusion.
    cmd_dataset_kwargs = cmd_args_dict.pop('dataset_kwargs')

    # 4c. Merge all *other* parameters: YAML defaults <- Command line overrides
    final_params = {**yaml_config, **cmd_args_dict}

    # 4c. Merge all *other* parameters: YAML defaults <- Command line overrides
    final_params = {**yaml_config, **cmd_args_dict}

    # 4d. Resolve the final value for dataset_kwargs (CLI string overrides YAML dict, 
    # unless CLI string is default "{}")
    if isinstance(yaml_dataset_kwargs, dict) and cmd_dataset_kwargs == "{}":
        # YAML dict provided, and command line used the default string value. Use YAML dict, but dump back to string.
        final_params['dataset_kwargs'] = json.dumps(yaml_dataset_kwargs)
    else:
        # Command line provided a custom value (even if it's the default "{}") or 
        # the YAML value was not a dict. Use the command-line value (string).
        final_params['dataset_kwargs'] = cmd_dataset_kwargs


    # 5. Final Namespace Generation with Correct Typing
    # The merged dictionary 'final_params' contains the final desired values.
    # To ensure all values are correctly cast (e.g., string '30' -> float 30.0),
    # we use set_defaults and parse an empty list.
    final_parser = parser.__class__(
        description=parser.description,
        formatter_class=parser.formatter_class
    )
    
    # # We need to re-add all actions from the original parser to the final_parser
    # # to maintain the correct argument definitions (types, etc.).
    # for action in parser._actions:
    #     if isinstance(action, (argparse._StoreAction, argparse._StoreConstAction, 
    #                            argparse._StoreTrueAction, argparse._StoreFalseAction,
    #                            argparse._AppendAction, argparse._CountAction,
    #                            argparse._HelpAction, argparse._VersionAction,
    #                            argparse._SubParsersAction, argparse._StoreStringAction)):
    #         # This is a bit complex, but essential for reusability. 
    #         # We skip adding the original actions here because the primary 
    #         # parser definition will handle this for us.
    #         pass

    # Use set_defaults on the original parser instance to force the final values
    parser.set_defaults(**final_params)
    
    # Parse an empty list. This triggers all the type conversion functions (like
    # _parse_normalize) and applies the defaults set above.
    final_args = parser.parse_args([])
    
    # Ensure mandatory fields are present, as 'required=True' is bypassed by parse_args([]).
    if not final_args.root_dir:
        print("Error: The 'root_dir' argument is required and was not provided via command line or config file.")
        sys.exit(1)
    # if not final_args.chan1:
    #     print("Error: The 'chan1' argument is required and was not provided via command line or config file.")
    #     sys.exit(1)
        
    return final_args







def crop_cells(
    mask_path: Union[str, Path],
    img_paths: Union[Union[str, Path], Iterable[Union[str, Path]]],
    cell_ids: Union[int, Iterable[int]],
    scale_factor: float = 65535.0,
    target_size: int = None,
    clip_mask: bool = False,
    pad_square: bool = True,
    rotate_horizontal: bool = False
) -> np.ndarray:
    """
    Crop multiple cells from multiple channel images using their instance segmentation mask IDs.
    
    Args:
        mask_path: Path to the instance segmentation mask (integer-labeled)
        img_paths: Path to the input image (generally multiple channels), mask will apply to all img_paths
        cell_ids: Single int or iterable (list, tuple, set, etc.) of cell IDs to crop
        scale_factor: normalized by dividing to scale_factor
        target_size: If provided, resize each cropped cell to (target_size, target_size)
        clip_mask: Zero out pixels outside the cell mask
        pad_square: Pad each crop to a square
        rotate_horizontal: Rotate each cell so its major axis is horizontal
    
    Returns:
        cropped cell images (np.ndarray, B-H-W-C), in the same order as input cell_ids.
        Failed/extra-missing cells return None.
    """
    try:
        imgs = np.stack([iio.imread(img_path) for img_path in img_paths], axis=2)
        mask = iio.imread(mask_path)
        if mask.ndim == 3:
            print("convert GRB mask")
            mask = mask[..., 0]  # handle RGBA masks

        # Convert single int to list for uniform processing
        if isinstance(cell_ids, int):
            cell_ids = [cell_ids]
        elif not isinstance(cell_ids, (list, tuple, set, np.ndarray)):
            cell_ids = list(cell_ids)

        results = []
        for cell_id in cell_ids:
            try:
                # Skip if cell ID not in mask
                if cell_id not in np.unique(mask):
                    results.append(None)
                    continue

                cell_mask = (mask == cell_id)
                if not cell_mask.any():
                    results.append(None)
                    continue

                # Get bounding box and orientation
                props = regionprops(label(cell_mask))[0]
                bbox = props.bbox  # ymin, xmin, ymax, xmax
                orientation = props.orientation

                y0, x0, y1, x1 = bbox
                cropped_mask = cell_mask[y0:y1, x0:x1]
                cropped_imgs = imgs[y0:y1, x0:x1, :]

                # # Ensure image has channel dimension
                # if cropped_img.ndim == 2:
                #     cropped_img = cropped_img[..., None]

                # Optional: clip image to mask
                if clip_mask:
                    cropped_imgs = cropped_imgs * cropped_mask[..., None]

                # Rotate to make major axis horizontal
                if rotate_horizontal:
                    max_dim = max(cropped_imgs.shape[:2])
                    pad_width = max_dim
                    pad_img = ((pad_width, pad_width), (pad_width, pad_width))
                    pad_mask = ((pad_width, pad_width), (pad_width, pad_width))

                    if cropped_imgs.ndim == 3:
                        cropped_imgs = np.pad(cropped_imgs, pad_img + ((0,0),), mode='constant', constant_values=0)
                    else:
                        cropped_imgs = np.pad(cropped_imgs, pad_img, mode='constant', constant_values=0)
                    cropped_mask = np.pad(cropped_mask, pad_mask, mode='constant', constant_values=0)

                    angle_degrees = -np.degrees(orientation) + 90
                    cropped_imgs = ndi.rotate(cropped_imgs, angle_degrees, reshape=False, order=1, mode='constant', cval=0)
                    cropped_mask = ndi.rotate(cropped_mask, angle_degrees, reshape=False, order=0, mode='constant', cval=0)

                    # Re-crop to tight bounding box after rotation
                    coords = np.argwhere(cropped_mask)
                    if len(coords) == 0:
                        results.append(None)
                        continue
                    y_min, x_min = coords.min(axis=0)
                    y_max, x_max = coords.max(axis=0)
                    cropped_imgs = cropped_imgs[y_min:y_max+1, x_min:x_max+1, :]

                # Pad to square
                if pad_square:
                    h, w = cropped_imgs.shape[:2]
                    size = max(h, w)
                    pad_h = (size - h) // 2
                    pad_w = (size - w) // 2
                    pad_img = ((pad_h, size - h - pad_h), (pad_w, size - w - pad_w))
                    if cropped_imgs.ndim == 3:
                        pad_img += ((0, 0),)
                    cropped_imgs = np.pad(cropped_imgs, pad_img, mode='constant', constant_values=0)

                # Resize to target size
                if target_size:
                    output_shape = (target_size, target_size) + cropped_imgs.shape[2:]
                    cropped_imgs = transform.resize(
                        cropped_imgs, output_shape, anti_aliasing=True, preserve_range=True).astype(imgs.dtype)

                results.append(cropped_imgs)

            except Exception as e:
                print(f"[ERROR] Failed to crop cell ID {cell_ids}: {e}")
                results.append(None)

        # list to array
        results = np.stack(results, axis=0)

        # normalize
        results = results / scale_factor

        return results

    except Exception as e:
        print(f"[FATAL ERROR] Failed to load images/mask: {e}")
        # Return list of None with same length as input cell_ids
        n = 1 if isinstance(cell_ids, int) else len(cell_ids)
        return [None] * n
    




def _get_images_to_dataset_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Find measurements",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", type=str, help="Path to a YAML configuration file to load parameters from.", default=None)
    parser.add_argument("--root_dir", type=str, help="The single Measurement directory to process, or a root containing it.")
    parser.add_argument("--dataset_kwargs", type=str, default="{}", 
                        help="JSON string of extra kwargs to pass to images_to_dataset(), e.g. "
                         "'{\"image_suffix\": \".tif\", \"mask_suffix\": \".png\", \"image_extractor\": \"...\"}' "
                         "Use null for None, true/false for booleans.")
    return parser


def main():
    """Main execution entry point for image_helper.py demonstration."""
    parser = _get_images_to_dataset_parser()
    print(parser)
    args = load_config_and_args(parser)
    print(args)
    dataset_kwargs = parse_dataset_kwargs(args.dataset_kwargs)
    print(dataset_kwargs)
    
    # Find measurement directory
    root_dir = Path(args.root_dir)
    measurement_dirs = find_measurement_dirs(root_dir)
    
    # If no measurement dirs found, check if the root_dir itself is a measurement
    if not measurement_dirs:
        print(f"Could not find any Measurement directory in or at {root_dir}")
        sys.exit(1)

    for mdir in measurement_dirs:
        try:
            # Merge user overrides and default settings
            final_kwargs = {"remove_na_row":True, "image_subdir":"Images", "image_suffix":".tiff", **dataset_kwargs}
            
            # Call the core function
            dataset = images_to_dataset(mdir, final_kwargs)
            
            if dataset and dataset.get("df") is not None:
                df = dataset["df"]
                print(df)
            else:
                print("images_to_dataset returned an empty or invalid dataset.")

        except Exception as e:
            print(f"Error while running images_to_dataset demo: {e}", exc_info=True)
            continue

    sys.exit(0)


if __name__ == "__main__":
    main()