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

import logging, sys, re, string, argparse, json
from typing import Union, Iterable, Dict, List, Optional, Any
from pathlib import Path
from glob import glob
import pandas as pd
import numpy as np
from natsort import natsorted
try:
    import imageio.v3 as iio
except ImportError:
    import imageio as iio
# import cv2
import scipy.ndimage as ndi
from skimage.transform import resize as sk_resize
from skimage.measure import regionprops, regionprops_table, label
from skimage.morphology import disk, closing
from tqdm import tqdm
import sqlite3

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
    measurement_pattern: str = r".*Measurement(?:\s*\d+)?$",
    image_subdir: str = "Images"
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
    pattern = re.compile(measurement_pattern, re.IGNORECASE) # re.IGNORECASE

    # Directly search for the 'Images' folder - this is efficient
    images_dirs = root.rglob(f"{image_subdir}")

    measurement_dirs = []
    for img_dir in images_dirs:
        parent = img_dir.parent
        # Check if parent name matches the measurement pattern
        if pattern.match(parent.name):
            measurement_dirs.append(parent)

    # Natural sort by full path (or just name if preferred)
    return natsorted(measurement_dirs, key=lambda p: str(p))


def parse_string_to_dict(s: str) -> dict:
    # Remove outer braces and whitespace
    s = s.strip().strip('{} ')
    if not s:
        return {}
    # Split by top-level commas (safe here because no values contain commas)
    pairs = s.split(',')
    result = {}
    for pair in pairs:
        if ':' not in pair:
            continue
        key, value = pair.split(':', 1)          # split only on first colon
        key = key.strip()
        value = value.strip()
        
        # Convert True/False/None to proper Python types
        if value == 'True':
            value = True
        elif value == 'False':
            value = False
        elif value == 'None':
            value = None
        # else: leave as string (including complex regex)
        
        result[key] = value
    
    # return {k: None if v is None else v for k, v in raw.items()}
    return result
    

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
        r"sk(?P<timepoint>[0-9]{1,})fk1fl1(?P<self_generated>.*)"
    ),
    mask_extractor: str = None, # if none, re-phrase from image_extractor
    # mask_extractor: str = (
    #     r"r(?P<row>.*)c(?P<column>.*)f(?P<field>.*)p(?P<stack>.*)-ch(?P<channel>[0-9]{1,})"
    #     r"sk(?P<timepoint>[0-9]{1,})fk1fl1_(cp_masks_)?(?P<mask_name>.*)"
    # ),
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
    # if not ("channel" in image_extract_cols or "self_generated" in image_extract_cols):
    #     raise ValueError("image_extractor must contain [self_generated, channel]")

    # ------------------------------------------------------------------- #
    # Intensity images
    # ------------------------------------------------------------------- #
    image_glob = str(measurement_dir / image_subdir / f"**/*{image_suffix}")
    image_files = glob(image_glob, recursive=True)
    image_files = [Path(p) for p in image_files]
    # print(len(image_files))

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
    # print(image_df)
    image_meta = image_df["filename"].str.extract(image_extractor + re.escape(image_suffix))

    if "channel" in image_meta.columns and "self_generated" in image_meta.columns:
        image_meta["channel"] = image_meta["channel"].astype(str) + image_meta["self_generated"]
    elif "channel" in image_meta.columns:
        image_meta["channel"] = image_meta["channel"].astype(str)
    elif "self_generated" in image_meta.columns:
        image_meta["channel"] = image_meta["self_generated"].astype(str)
    else:
        image_meta["channel"] = "0"
    
    image_meta.loc[image_meta["channel"] == "", "channel"] = "0"
    if "self_generated" in image_meta.columns: image_meta = image_meta.drop(columns=["self_generated"])

    image_df = pd.concat([image_df, image_meta], axis=1)

    # Index columns (everything except channel)
    image_index_cols = ["directory"] + [c for c in image_extract_cols if c not in ("channel", "self_generated")]
    # set index for merging with mask
    image_df = image_df.set_index(image_index_cols)
    image_df["channel"] = "ch" + image_df["channel"]
    image_colnames = natsorted(image_df["channel"].unique().tolist())

    # Pivot to wide format
    # print(image_df)
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

            if mask_extractor is None:
                mask_extractor = image_extractor.replace("(?P<self_generated>.*)", "_cp_masks_(?P<mask_name>.*)")
                # print("set mask_extractor:", mask_extractor)
            
            mask_meta = mask_df["filename"].str.extract(mask_extractor + re.escape(mask_suffix))
            # print(mask_meta)
           
            if "channel" not in mask_meta.columns:
                mask_meta["channel"] = "0"
            
            mask_df = pd.concat([mask_df, mask_meta], axis=1)
            mask_df = mask_df.set_index(image_index_cols)

            mask_colnames = natsorted(mask_df["mask_name"].unique().tolist())
            # print(mask_df)
            mask_df = mask_df.pivot(columns="mask_name", values="filename")
            print(f"{len(mask_df)} grouped object masks: {mask_colnames}")
            df = pd.merge(image_df, mask_df, how="left", left_index=True, right_index=True)
            
    df = df.reset_index(drop=False)  # keep index columns as regular columns
    # to int if possible
    # # print(df["row"].dtype)
    
    # print(df)
    for col in df.select_dtypes(include='object').columns:
        try:
            df[col] = pd.to_numeric(df[col], errors='raise')
        except (ValueError, TypeError):
            pass
    # print(df["row"].dtype)
    # # print(df["directory"].dtype)
    # print(df)

    # ------------------------------------------------------------------- #
    # Clean-up
    # ------------------------------------------------------------------- #
    if remove_na_row:
        # print("removing rows with NA")
        df = df.dropna()
    # print(df)
    # print(df.columns)
    # print(df.index.name)

    print(f"{len(df)} groups after merging images and masks")

    metadata_colnames = image_index_cols.copy()
    # print(metadata_colnames)

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



def measure_image(
        mask: np.ndarray, # shape HW
        img: np.ndarray, # shape HWC
        properties: list = ['label','bbox','area_filled','perimeter',
                            'equivalent_diameter_area','eccentricity','solidity',
                            'intensity_mean','inertia_tensor_eigvals'], # list of regionprops properties
        extra_properties: list = None # list of function
        ) -> pd.DataFrame:
    """
    Measure image properties
    Args:
        mask: np.ndarray, shape HW
        img: np.ndarray, shape HWC
        properties: list, list of regionprops properties
    Returns:
        pd.DataFrame, dataframe with image properties
    """
    # check input
    assert mask.ndim == 2, 'mask must be 2D'
    img = img if img.ndim == 3 else img[:,:,None]
    assert img.shape[2] < img.shape[0], 'img must be HWC format'
    
    prop = regionprops_table(mask, img, 
                             properties=properties,
                             extra_properties=extra_properties)
    prop = pd.DataFrame(prop)

    # convert dict to df
    if extra_properties is not None:
        extra_properties = [x.__name__ for x in extra_properties]
        for func in extra_properties:
            if func in prop.columns:
                pd_new = pd.DataFrame()
                for obj_idx, row_data in prop.iterrows():
                    data = row_data[func]
                    for item in list(data[0].keys()):
                        feature_len = len(data[0][item])
                        for ch_idx in range(len(data)):
                            for feature_idx in range(feature_len):
                                pd_new.loc[obj_idx, f'{item}_{feature_idx+1}-{ch_idx}'] = data[ch_idx][item][feature_idx]
            prop = prop.drop(columns=[func])
            prop = pd.merge(prop, pd_new, left_index=True, right_index=True)
    
    return prop

def measure_dataset(
        res, # resulf from images_to_dataset
        save_path: str = './result.db',
        **kwargs # kwargs for measure_image
        ) -> dict:
    """
    Measure dataset
    Args:
        res: result from images_to_dataset
        save_path: str, path to save the result
        **kwargs: kwargs for measure_image
    """
    df = res['df']
    intensity_colnames = res['intensity_colnames']
    mask_colnames = res['mask_colnames']

    # config db
    # fname = join(dirname(df['directory'].iloc[0]), save_db_name)
    save_path = Path(save_path)
    if save_path.exists():
        print(f'{save_path} exist, exit!')
        return None
    else:
        conn = sqlite3.connect(save_path)
        print(f'save measurement: {save_path}')
    
    props_dict = {}
    # measure each mask
    for mask_name in mask_colnames:
        try:
            props = []
            for idx, data in tqdm(df.iterrows(), total=len(df)):
            
                img = [iio.imread(Path(data['directory']) / data[x]) for x in intensity_colnames]
                img = np.stack(img, axis=2) # skimage regionprops accept imgs as HWC
                mask = iio.imread(Path(data['directory']) / data[mask_name])
                prop = measure_image(mask, img, **kwargs)

                # rename channel
                # some columns need to be renamed first in case of naming error
                prop.rename(columns={'area_filled':'area', 'equivalent_diameter_area':'diameter'}, inplace=True)
                prop.rename(columns=lambda x: re.sub('bbox-','bbox_',x), inplace=True)
                prop.rename(columns=lambda x: re.sub('eigvals-','eigvals_',x), inplace=True)
                for ch_idx, ch in enumerate(intensity_colnames):
                    prop.rename(columns=lambda x: x.replace(f'-{ch_idx}', f'-{ch}'), inplace=True)
                
                # insert metadata
                meta = data.to_dict()
                [meta.pop(col, None) for col in intensity_colnames + mask_colnames]
                prop = prop.assign(**meta)

                # save to db
                if save_path:
                    table_name = mask_name.replace('cp_masks_', '')
                    prop.to_sql(table_name, conn, if_exists='append', index=False)
    
                props.append(prop)

            props = pd.concat(props)

            props_dict[mask_name] = props

        except Exception as e: 
            print(e)
        finally: 
            next

    # close db
    if save_path is not None: 
        conn.close()

    return props_dict



def granularity(mask, image, feature_size=16, step=4, target_size=128):
    """
    Calculate granularity spectrum for a region in an image.
    Compatible with regionprops extra_properties.
    
    Parameters
    ----------
    mask : ndarray
        Binary mask of the region
    image : ndarray 
        Intensity image
    feature_size : int
        Maximum size of features to analyze
    step : int
        Step size between element sizes
    target_size : int
        Resize image and mask to this size, this is used to keep the same resolution for all images
        
    Returns
    -------
    list
        Granularity spectrum measurements
    """

    if target_size:
        image = sk_resize(image, (target_size , target_size), anti_aliasing=True, preserve_range=True)
        mask = sk_resize(mask, (target_size, target_size), order=0, anti_aliasing=False, preserve_range=True)

    # Create mask for the region
    mask = mask.astype(bool)

    # # Scale down image and mask if requested
    # if size_scale < 1.0:
    #     new_shape = (np.array(image.shape) * size_scale).astype(int)
    #     i, j = (np.mgrid[0:new_shape[0], 0:new_shape[1]].astype(float) / size_scale)
    #     image = map_coordinates(image, (i, j), order=1)
    #     mask = map_coordinates(mask.astype(float), (i, j)) > 0.9
    #     feature_size = int(feature_size * size_scale)
    
    # Get image data within mask
    pixels = image * mask
    
    # Generate element sizes
    element_sizes = np.linspace(2, feature_size, step, dtype=int)
    spectrum = []
    
    # Calculate initial mean intensity
    start_mean = np.mean(pixels[mask])
    if start_mean <= np.finfo(float).eps:
        return [0.0] * len(element_sizes)
    
    # Calculate granularity spectrum
    prev_mean = start_mean
    for size in element_sizes:
        # Create structuring element
        selem = disk(size)
        
        # Perform morphological opening
        opened = ndi.grey_opening(pixels, structure=selem)
        
        # Calculate mean intensity after opening
        curr_mean = np.mean(opened[mask])
        
        # Calculate normalized difference
        diff = (prev_mean - curr_mean) * 100.0 / start_mean
        spectrum.append(max(0, diff))
        
        prev_mean = curr_mean
    
    # element_sizes = [int(x / size_scale) for x in element_sizes]

    return {"granularity": spectrum}
    # return spectrum



def radial_distribution(mask, image, num_bins=5, display=False):
    """
    Calculate radial intensity distribution for a region in an image.
    Compatible with regionprops extra_properties.
    The radial bins follow the shape of the mask, with equal distances from centroid to edge.
    
    Parameters
    ----------
    mask : ndarray
        Binary mask of the region
    image : ndarray
        Intensity image
    num_bins : int
        Number of radial bins to divide the region into
    display_radial : bool
        If True, displays the radial regions and their measured values on the image
        
    Returns
    -------
    dict
        Dictionary containing mean and total intensity values for each radial bin
    """
    
    # # Get centroid of the region
    # y_indices, x_indices = np.nonzero(mask)
    # center_y = np.mean(y_indices)
    # center_x = np.mean(x_indices)
    # # print(center_y, center_x)
    
    # A better way to get visual centroid 
    # https://stackoverflow.com/questions/1203135/what-is-the-fastest-way-to-find-the-visual-center-of-an-irregularly-shaped-pol
    # mask = draw_contour_on_mask((H,W), cnt)
    
    dist_img = ndi.distance_transform_edt(mask)
    center_y, center_x = np.where(dist_img == dist_img.max())
    center_y, center_x = center_y.mean(), center_x.mean()
    
    # # cv2.distanceTransform don't accept bool type
    # if mask.dtype == bool:
    #     mask = mask.astype(np.uint8)
    # dist_img = cv2.distanceTransform(mask, distanceType=cv2.DIST_L2, maskSize=5).astype(np.float32)
    # center_y, center_x = np.where(dist_img==dist_img.max())
    # center_y, center_x = center_y.mean(), center_x.mean() # there are sometimes cases where there are multiple values returned for the visual center
    # # print(center_y, center_x)
    
    # Calculate distances from center for all pixels in mask
    y, x = np.ogrid[:image.shape[0], :image.shape[1]]
    distances = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    
    # For each angle, find the distance to the mask edge
    angles = np.arctan2(y - center_y, x - center_x)
    edge_distances = np.zeros_like(angles)
    
    # For each pixel in the mask, update the max distance for its angle
    mask = mask.astype(bool)
    mask_coords = np.where(mask)
    pixel_angles = angles[mask_coords]
    pixel_distances = distances[mask_coords]
    
    # Bin angles into discrete values for edge distance calculation
    angle_bins = 180
    angle_edges = np.linspace(-np.pi, np.pi, angle_bins+1)
    binned_angles = np.clip(np.digitize(pixel_angles, angle_edges) - 1, 0, angle_bins-1)
    
    # Find max distance for each angle bin
    max_distances = np.zeros(angle_bins)
    for i in range(angle_bins):
        angle_mask = binned_angles == i
        if np.any(angle_mask):
            max_distances[i] = np.max(pixel_distances[angle_mask])
    
    # Fill in any gaps in max_distances using interpolation
    valid_indices = np.nonzero(max_distances)[0]
    if len(valid_indices) > 1:
        invalid_indices = np.where(max_distances == 0)[0]
        max_distances[invalid_indices] = np.interp(
            invalid_indices, 
            valid_indices,
            max_distances[valid_indices],
            period=angle_bins
        )
    
    # Interpolate edge distances for all angles
    all_angles = angles[mask]
    all_binned = np.clip(np.digitize(all_angles, angle_edges) - 1, 0, angle_bins-1)
    edge_distances[mask] = max_distances[all_binned]
    
    # Scale distances relative to edge distance
    scaled_distances = np.zeros_like(distances)
    nonzero_mask = mask  # Changed to use full mask instead of checking edge_distances
    scaled_distances[nonzero_mask] = distances[nonzero_mask] / edge_distances[nonzero_mask]
    
    # Create bins from 0 to 1
    bin_edges = np.linspace(0, 1, num_bins + 1)
    mean_values = np.zeros(num_bins)
    total_values = np.zeros(num_bins)
    
    if display:
        import matplotlib.pyplot as plt
        display_img = np.zeros_like(image)
    
    # Calculate mean and total intensity in each bin
    for i in range(num_bins):
        bin_mask = (scaled_distances >= bin_edges[i]) & (scaled_distances < bin_edges[i+1]) & mask
        if np.any(bin_mask):
            mean_values[i] = np.mean(image[bin_mask])
            total_values[i] = np.sum(image[bin_mask])
            
            if display:
                display_img[bin_mask] = mean_values[i]
       
    if display:
        plt.figure(figsize=(10, 5))
        plt.subplot(121)
        plt.imshow(image * mask)
        plt.title('Original Image (Masked)')
        plt.colorbar()
        
        plt.subplot(122)
        plt.imshow(display_img)
        plt.title('Radial Regions with Mean Values')
        plt.colorbar()
        plt.show()
            
    return { "mean_frac": mean_values.tolist(), "total_frac": total_values.tolist() }



def crop_cells(
    mask: Union[str, Path, np.ndarray],
    imgs: Union[None, str, Path, List[Union[str, Path]], np.ndarray] = None,
    cell_ids: Union[int, Iterable[int], None] = None,
    scale_factor: Union[float, None] = 65535.0,
    target_size: Optional[int] = None,
    clip_mask: bool = False,
    pad_square: bool = True,
    rotate_horizontal: bool = False
) -> List[Dict[str, Any]]:
    """
    Crop cells from intensity images using an instance segmentation mask.
    Always returns a list of dictionaries with 'cell_id', 'cell_img' and 'cell_mask' keys.

    Args:
        mask: Path to mask or already loaded integer mask np.ndarray (H, W)
        imgs: 
            - None → only return masks
            - Path / str → single channel
            - List of Path / str → multiple channels
            - np.ndarray (H, W, C) or (H, W) → pre-loaded image(s)
        cell_ids: Cell ID(s) to crop. If None, crop all non-background cells.
        scale_factor: Divide image by this value for normalization (None to skip)
        target_size: Resize crops to (target_size, target_size)
        clip_mask: Zero out pixels outside cell mask in intensity image
        pad_square: Pad crops to square shape
        rotate_horizontal: Rotate cell so major axis is horizontal

    Returns:
        List of dicts: [{'cell_id': int, 'cell_img': np.ndarray or None, 'cell_mask': np.ndarray (uint8 binary)}, ...]
        cell_img shape is (H, W, C), if only one channel, then (H, W, 1), float if scale_factor set
        cell_mask shape is (H, W), uint8 type
        Failed cells return {'cell_id': None, 'cell_img': None, 'cell_mask': None}
        Order matches input cell_ids (or order in mask if cell_ids=None)
    """
    # Load mask
    if isinstance(mask, (str, Path)):
        mask = iio.imread(mask)
    if mask.ndim != 2:
        raise ValueError('mask must be 2D')
    # mask = mask.astype(np.uint16, copy=False)
    mask_arr = mask

    # Load images if provided
    img_arr: Optional[np.ndarray] = None
    if imgs is not None:
        if isinstance(imgs, (str, Path)):
            imgs = [imgs]
        if isinstance(imgs, list):
            channels = [iio.imread(f) for f in imgs]
            img_arr = np.stack(channels, axis=-1)  # (H, W, C)
        elif isinstance(imgs, np.ndarray):
            if imgs.ndim == 2:
                img_arr = imgs[..., np.newaxis]
            elif imgs.ndim == 3:
                img_arr = imgs
            else:
                raise ValueError("Image array must be 2D or 3D, channel is last (H-W-C)")
        else:
            raise TypeError("imgs must be path, list of paths, or np.ndarray")

        if img_arr.shape[:2] != mask_arr.shape:
            raise ValueError(f"Image shape {img_arr.shape[:2]} != mask shape {mask_arr.shape}")

    # Determine cell IDs
    if cell_ids is None:
        cell_ids = np.unique(mask_arr[mask_arr != 0])
    elif isinstance(cell_ids, int):
        cell_ids = [cell_ids]
    else:
        cell_ids = list(cell_ids)

    results = []
    for cell_id in cell_ids:
        try:
            if cell_id not in mask_arr:
                results.append({'img': None, 'mask': None})
                continue

            cell_mask_bool = (mask_arr == cell_id)
            if not cell_mask_bool.any():
                results.append({'img': None, 'mask': None})
                continue

            # Initial bounding box
            props = regionprops(label(cell_mask_bool))[0]
            y0, x0, y1, x1 = props.bbox
            orientation = props.orientation

            # Initial crops
            cropped_mask_bool = cell_mask_bool[y0:y1, x0:x1]
            cropped_mask_uint8 = cropped_mask_bool.astype(np.uint8)
            cropped_img = img_arr[y0:y1, x0:x1, :] if img_arr is not None else None

            # Clip image early
            if cropped_img is not None and clip_mask:
                cropped_img = cropped_img * cropped_mask_bool[..., None]

            # === Rotate to horizontal ===
            if rotate_horizontal:
                angle_deg = -np.degrees(orientation) + 90

                h, w = cropped_mask_bool.shape
                max_dim = max(h, w)
                pad_total = max_dim * 2
                pad_h = (pad_total - h) // 2
                pad_w = (pad_total - w) // 2
                pad_kwargs = ((pad_h, pad_total - h - pad_h), (pad_w, pad_total - w - pad_w))

                cropped_mask_bool = np.pad(cropped_mask_bool, pad_kwargs, mode='constant', constant_values=False)
                cropped_mask_uint8 = np.pad(cropped_mask_uint8, pad_kwargs, mode='constant', constant_values=0)
                if cropped_img is not None:
                    cropped_img = np.pad(cropped_img, pad_kwargs + ((0, 0),), mode='constant', constant_values=0)

                # Rotate
                if cropped_img is not None:
                    cropped_img = ndi.rotate(cropped_img, angle_deg, reshape=False, order=1,
                                             mode='constant', cval=0.0)
                cropped_mask_bool = ndi.rotate(cropped_mask_bool, angle_deg, reshape=False, order=0,
                                               mode='constant', cval=False)
                cropped_mask_uint8 = ndi.rotate(cropped_mask_uint8, angle_deg, reshape=False, order=0,
                                                mode='constant', cval=0)

                # Tight recrop
                coords = np.column_stack(np.where(cropped_mask_bool))
                if len(coords) == 0:
                    results.append({'img': None, 'mask': None})
                    continue
                ymin, xmin = coords.min(axis=0)
                ymax, xmax = coords.max(axis=0)
                cropped_mask_bool = cropped_mask_bool[ymin:ymax+1, xmin:xmax+1]
                cropped_mask_uint8 = cropped_mask_uint8[ymin:ymax+1, xmin:xmax+1]
                if cropped_img is not None:
                    cropped_img = cropped_img[ymin:ymax+1, xmin:xmax+1]

            # === Pad to square ===
            if pad_square:
                h, w = cropped_mask_bool.shape
                size = max(h, w)
                pad_h = (size - h) // 2
                pad_w = (size - w) // 2
                pad_kwargs = ((pad_h, size - h - pad_h), (pad_w, size - w - pad_w))

                cropped_mask_bool = np.pad(cropped_mask_bool, pad_kwargs, mode='constant', constant_values=False)
                cropped_mask_uint8 = np.pad(cropped_mask_uint8, pad_kwargs, mode='constant', constant_values=0)
                if cropped_img is not None:
                    cropped_img = np.pad(cropped_img, pad_kwargs + ((0, 0),), mode='constant', constant_values=0)

            # === Resize ===
            if target_size is not None:
                new_shape = (target_size, target_size)
                if cropped_img is not None:
                    cropped_img = sk_resize(cropped_img, new_shape + cropped_img.shape[2:],
                                            anti_aliasing=True, preserve_range=True).astype(cropped_img.dtype)
                cropped_mask_bool = sk_resize(cropped_mask_bool, new_shape, order=0, anti_aliasing=False, preserve_range=True) > 0.5
                cropped_mask_uint8 = sk_resize(cropped_mask_uint8, new_shape, order=0, anti_aliasing=False, preserve_range=True)
                cropped_mask_uint8 = (cropped_mask_uint8 > 0.5).astype(np.uint8)

            # === Normalize image ===
            if cropped_img is not None and scale_factor is not None and scale_factor > 0:
                cropped_img = cropped_img / scale_factor

            results.append({
                'cell_id': cell_id,
                'cell_img': cropped_img,
                'cell_mask': cropped_mask_uint8
            })

        except Exception as e:
            print(f"[ERROR] Failed to crop cell ID {cell_id}: {e}")
            results.append({"cell_id": None, 'cell_img': None, 'cell_mask': None})

    return results



def normalize_img(img: np.array,
                  method: str = "percentile",
                  percentile_value: tuple = (0.1, 99.9),
                  channel_index: int = None):
    """
    Normalizes an image, optionally by channel.

    Args:
        img: The input NumPy array representing the image.
        method: The normalization method ("scale", "percentile", "minmax", "mean", or None).
        percentile_value: The percentile values for the "percentile" method.
        channel_index: The index of the channel to normalize, or None to normalize the entire image.
    """
    if method not in ["scale", "percentile", "minmax", "mean", None]:
        raise ValueError(f'method should be one of ["scale","percentile","minmax","mean",None]')

    if channel_index is not None:
        if channel_index >= img.shape[-1] or channel_index < 0:
            raise ValueError(f"channel_index {channel_index} is out of range for image with shape {img.shape}")
        
        channel_img = img[..., channel_index]
        mask = channel_img > 0
        nonzero_values = channel_img[mask]
        
        if len(nonzero_values) == 0:
            return img  # If all zeros in channel, return original
        
        if method == 'mean':
            mean = np.mean(nonzero_values)
            if mean == 0: mean = 1.0
            img_norm_channel = channel_img / mean
            img_norm = img.copy()
            img_norm[..., channel_index] = img_norm_channel
            
        elif method == "scale":
            mean = np.mean(nonzero_values)
            std = np.std(nonzero_values)
            if std == 0: std = 1.0
            img_norm_channel = (channel_img - mean) / std
            img_norm = img.copy()
            img_norm[..., channel_index] = img_norm_channel

        elif method == "percentile":
            p_low = np.percentile(nonzero_values, percentile_value[0])
            p_high = np.percentile(nonzero_values, percentile_value[1])
            if p_high - p_low > 0:
                img_norm_channel = np.clip((channel_img - p_low) / (p_high - p_low), 0, 1)
            else:
                img_norm_channel = channel_img
            img_norm = img.copy()
            img_norm[..., channel_index] = img_norm_channel

        elif method == "minmax":
            min_val = np.min(nonzero_values)
            max_val = np.max(nonzero_values)
            if max_val - min_val > 0:
                img_norm_channel = (channel_img - min_val) / (max_val - min_val)
            else:
                img_norm_channel = channel_img
            img_norm = img.copy()
            img_norm[..., channel_index] = img_norm_channel

        else:
            img_norm = img

    else:  # Normalize entire image
        mask = img > 0
        nonzero_values = img[mask]

        if len(nonzero_values) == 0:
            return img  # If all zeros, return original
        
        if method == 'mean':
            mean = np.mean(nonzero_values)
            if mean == 0: mean = 1.0
            img_norm = img / mean
            
        elif method == "scale":
            mean = np.mean(nonzero_values)
            std = np.std(nonzero_values)
            if std == 0: std = 1.0
            img_norm = (img - mean) / std

        elif method == "percentile":
            p_low = np.percentile(nonzero_values, percentile_value[0])
            p_high = np.percentile(nonzero_values, percentile_value[1])
            if p_high - p_low > 0:
                img_norm = np.clip((img - p_low) / (p_high - p_low), 0, 1)
            else:
                img_norm = img

        elif method == "minmax":
            min_val = np.min(nonzero_values)
            max_val = np.max(nonzero_values)
            if max_val - min_val > 0:
                img_norm = (img - min_val) / (max_val - min_val)
            else:
                img_norm = img
                
        else:
            img_norm = img

    return img_norm


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate cp_dataloader.csv for CellProfiler in multiple Measurement folders",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("root_dir", type=str, help="Root folder containing *_Measurement X folders")
    parser.add_argument("--dry_run", action="store_true", help="Just display but not write anything")
    parser.add_argument("--cp_dataloder", action="store_true", help="write cellprofiler dataloader style")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing cp_dataloader.csv")
    parser.add_argument("--test", action="store_true", help="Generate dataloader with ~100 random rows per measurement")
    parser.add_argument("--dataset_kwargs", type=str, default="{}",
                        help="JSON string with overrides for images_to_dataset(). "
                            "Supports Python literals! Examples:\n"
                            "  --dataset_kwargs '{\"image_suffix\": \".tif\"}'\n"
                            "  --dataset_kwargs '{\"mask_suffix\": None, \"remove_na_row\": False}'\n"
                            "  --dataset_kwargs '{\"cellprofiler_style\": True}'")
    args = parser.parse_args()
    user_dataset_kwargs = parse_string_to_dict(args.dataset_kwargs)
    # print(user_dataset_kwargs)

    root_dir = Path(args.root_dir)
    if not root_dir.is_dir():
        print(f"Root directory not found: {root_dir}")
        sys.exit(1)

    # Discover measurements
    measurements = find_measurement_dirs(root_dir)
    if not measurements:
        print("No *_Measurement X folders found")
        sys.exit(1)


    for mdir in measurements:
        print("")
        print(f"--- Processing {mdir.name} ---")
        # try:
        final_kwargs = {"cellprofiler_style": args.cp_dataloder, **user_dataset_kwargs}
        # print(final_kwargs)
        dataset = images_to_dataset(mdir, **final_kwargs)

        if not dataset or dataset["df"].empty:
            print(f"No valid images found in {mdir.name}")
            continue

        df = dataset["df"]
        total_rows = len(df)

        if args.test:
            df = df.sample(n=min(100, total_rows), random_state=42)
            print(f"TEST MODE: using {len(df)} random rows (out of {total_rows})")

        # Write dataloader CSV
        if args.cp_dataloder and not args.dry_run:
            output_path = mdir / "cp_dataloader.csv"
            if output_path.exists() and not args.overwrite:
                print(f"cp_dataloader.csv already exists → skip {mdir.name}")
                continue
            df.to_csv(output_path, index=False)
            print(f"Wrote {len(df)} rows to {output_path}")

        # except Exception as e:
        #     print(f"Failed to generate dataset for {mdir.name}: {e}")
        #     continue


if __name__ == "__main__":
    main()