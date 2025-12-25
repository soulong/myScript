import sys
sys.path.append('/mnt/c/Users/haohe/Documents/GitHub/myScript/python/image_analysis')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import imageio.v3 as iio
from scipy import stats # Needed for robust transformation
from skimage.util import img_as_float # Good practice for normalization

import sqlite3
import json
from datetime import datetime
from pathlib import Path
from typing import Union, List, Tuple, Dict, Any, Optional

from image_helper import images_to_dataset, crop_cells


# =========================================================================
# 1. CORE ANALYSIS FUNCTION: intensity_bins_analysis
# =========================================================================

def intensity_bins_analysis(
    cell_data_list: List[Dict[str, Any]], 
    n_bins: int = 10,
    normalize: Optional[str] = None,  # None, 'mean', 'median', or 'quantile'
    bin_method: str = 'linear',  
    manual_edges: Optional[Union[List[float], np.ndarray]] = None,
    intensity_range: Tuple[float, float] = None, # (min_val, max_val)
    log_transform: bool = False # apply log transformation before normalization
) -> List[Dict[str, Any]]:
    """
    Analyze individual cell pixel intensity distribution across defined bins.

    - Now supports 'quantile' normalization using q0.01 and q0.99 clipping/scaling.
    """
    if normalize not in {None, 'mean', 'median', 'quantile'}:
        raise ValueError("normalize must be None, 'mean', 'median', or 'quantile'")

    results = []
    for data in cell_data_list: 
        cell_id = data.get('cell_id') 
        img = data.get('cell_img') 
        mask = data.get('cell_mask') 
        
        if img is None or mask is None:
            print(f"Skipping cell {cell_id or 'unknown'}: missing 'img' or 'mask' key.")
            continue

        gray = img_as_float(img[:,:,0]) 
        h, w = gray.shape[:2]

        cell_mask = (mask > 0).ravel()
        flat = gray.ravel()

        valid_pixels_orig = flat[cell_mask] 
        
        # --- Normalization ---
        norm_factor_mean_median = 1.0       
        norm_min_q = None               # Quantile min (q0.01)
        norm_max_q = None             # Quantile max (q0.99)
        normalized_pixels = valid_pixels_orig

        if log_transform:
            valid_pixels_orig = np.log(valid_pixels_orig + 1e-8)
            
        if len(valid_pixels_orig) > 0:
            if normalize == 'mean':
                mean_val = valid_pixels_orig.mean()
                norm_factor_mean_median = mean_val if mean_val > 1e-8 else 1.0
                normalized_pixels = valid_pixels_orig / norm_factor_mean_median
            
            elif normalize == 'median':
                median_val = np.median(valid_pixels_orig)
                norm_factor_mean_median = median_val if median_val > 1e-8 else 1.0
                normalized_pixels = valid_pixels_orig / norm_factor_mean_median
            
            elif normalize == 'quantile':
                # Calculate the 1st and 99th percentiles
                q01, q99 = np.percentile(valid_pixels_orig, (0.1, 99.9))
                norm_min_q = float(q01)
                norm_max_q = float(q99)

                norm_range_q = q99 - q01

                if norm_range_q > 1e-8:
                    # 1. Clip pixels to the [q01, q99] range
                    clipped_pixels = np.clip(valid_pixels_orig, q01, q99)
                    # 2. Normalize by scaling the clipped range to [0, 1]
                    normalized_pixels = (clipped_pixels - q01) / norm_range_q
                else:
                    # All pixels are essentially the same value
                    normalized_pixels = np.full_like(valid_pixels_orig, 0.5) 
                    norm_max_q = norm_min_q + 1e-8 # Ensure q99 > q01 in metadata
                    
            # Apply user-defined intensity range clipping after normalization
            # Clip normalized pixels to user-defined intensity range
            if intensity_range is not None:
                min_val, max_val = intensity_range
                normalized_pixels = np.clip(normalized_pixels, min_val, max_val)
            else:
                min_val = np.min(normalized_pixels)
                max_val = np.max(normalized_pixels)
            
            # --- Bin edge generation (moved to after normalization) ---
            if bin_method == 'manual':
                if manual_edges is None or len(manual_edges) < 2:
                    raise ValueError("manual_edges must be provided for 'manual' bin_method")
                bin_edges = np.array(manual_edges)
            else:
                if bin_method not in {'linear', 'sqrt', 'exp', 'logspace', 'geomspace'}:
                    raise ValueError("bin_method must be 'linear', 'sqrt', 'exp', 'logspace', 'geomspace', or 'manual'")
                    
                # Use the clipped range for bin generation
                actual_min = max(min_val, normalized_pixels.min())
                actual_max = min(max_val, normalized_pixels.max())
                
                inner_min = max(1e-8, actual_min) if bin_method in {'exp', 'logspace', 'geomspace'} else actual_min
                inner_max = actual_max
                
                if bin_method == 'linear':
                    bin_edges = np.linspace(inner_min, inner_max, n_bins + 1)
                elif bin_method == 'sqrt':
                    bin_edges = np.linspace(np.sqrt(inner_min), np.sqrt(inner_max), n_bins + 1) ** 2
                elif bin_method == 'exp':
                    bin_edges = np.exp(np.linspace(np.log(inner_min), np.log(inner_max), n_bins + 1))
                elif bin_method == 'logspace':
                    bin_edges = np.logspace(np.log10(inner_min), np.log10(inner_max), n_bins + 1)
                elif bin_method == 'geomspace':
                    bin_edges = np.geomspace(inner_min, inner_max, n_bins + 1)

                bin_edges[0] = actual_min
                bin_edges[-1] = actual_max
            
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            # --- Binning Logic ---
            bin_indices_flat_valid = np.digitize(normalized_pixels, bin_edges) - 1
            bin_indices_flat_valid = np.clip(bin_indices_flat_valid, 0, n_bins - 1)
            
            bin_map_flat = np.full(flat.size, -1, dtype=np.int8)
            bin_map_flat[cell_mask] = bin_indices_flat_valid
            bin_map = bin_map_flat.reshape(h, w)
        else:
            bin_map = np.full((h, w), -1, dtype=np.int8)
            bin_edges = np.array([min_val, max_val])  # Minimal bin edges for empty case
            bin_centers = np.array([(min_val + max_val) / 2])
        
        # if cell_id == 1:
        #     # print(n_bins)
        #     # print(bin_edges)
        #     formatted_list_str = [f"{num:.6f}" for num in bin_edges]
        #     print(f"[{', '.join(formatted_list_str)}]")
        
        # --- Statistics Calculation ---
        pixel_counts = np.zeros(n_bins, dtype=np.int64)
        mean_intensities = np.full(n_bins, np.nan)
        
        if len(valid_pixels_orig) > 0:
            for i in range(n_bins):
                mask_in_bin = (bin_map == i) 
                orig_in_bin = gray[mask_in_bin] 
                if len(orig_in_bin) > 0:
                    pixel_counts[i] = len(orig_in_bin)
                    mean_intensities[i] = orig_in_bin.mean()

        mean_int = float(valid_pixels_orig.mean()) if len(valid_pixels_orig) > 0 else 0.0
        median_int = float(np.median(valid_pixels_orig)) if len(valid_pixels_orig) > 0 else 0.0
        
        # Collect all results and metadata
        results.append({
            'cell_id': cell_id,
            'gray': gray,
            'bin_map': bin_map,
            'bin_edges': bin_edges.copy(),
            'bin_centers': bin_centers.copy(),
            'pixel_counts': pixel_counts,
            'mean_intensities': mean_intensities,
            # Metadata for DB export
            'valid_pixel_count': len(valid_pixels_orig),
            'normalize': normalize,
            'norm_factor_mean_median': norm_factor_mean_median, 
            'norm_min_q': norm_min_q,    # NEW: Quantile min (q0.01)
            'norm_max_q': norm_max_q,    # NEW: Quantile max (q0.99)
            'bin_method': bin_method,
            'intensity_range_min': min_val,
            'intensity_range_max': max_val,
            'n_bins': n_bins,
            'mean_intensity_cell': mean_int,
            'median_intensity_cell': median_int
        })
        
    return results


# =========================================================================
# 2. PLOTTING FUNCTIONS
# =========================================================================

def plot_intensity_bins(results: List[Dict[str, Any]], title_prefix: str = "Cell"):
    """Plot the pixel count and mean intensity distribution for all cells in results."""
    
    n_cells = len(results)
    if n_cells == 0:
        print("No results to plot.")
        return

    # Use a dynamic grid layout for multiple cells
    n_rows = int(np.ceil(n_cells / 3))
    n_cols = min(n_cells, 3)

    fig, axes = plt.subplots(n_rows, n_cols * 2, figsize=(4 * n_cols, 4 * n_rows), 
                             squeeze=False, gridspec_kw={'width_ratios': [1] * n_cols + [1] * n_cols})
    axes = axes.flatten()

    for i, res in enumerate(results):
        cell_id = res['cell_id']
        n_bins = res['n_bins']
        bin_centers = res['bin_centers']
        pixel_counts = res['pixel_counts']
        mean_intensities = res['mean_intensities']
        
        # Pixel Count Plot (Left column)
        ax_count = axes[i]
        ax_count.bar(bin_centers, pixel_counts, width=(bin_centers[1]-bin_centers[0])*0.8, color='skyblue', edgecolor='black')
        ax_count.set_title(f"{title_prefix} {cell_id or i+1}: Pixel Count")
        ax_count.set_xlabel(f"Normalized Intensity (Bin Center)")
        ax_count.set_ylabel("Pixel Count")
        ax_count.set_xticks(bin_centers)
        ax_count.tick_params(axis='x', rotation=45)
        
        # Mean Intensity Plot (Right column)
        ax_mean = axes[i + n_cols]
        # Filter out NaN mean intensities for plotting
        valid_centers = bin_centers[~np.isnan(mean_intensities)]
        valid_means = mean_intensities[~np.isnan(mean_intensities)]
        
        ax_mean.bar(valid_centers, valid_means, width=(valid_centers[1]-valid_centers[0])*0.8, color='salmon', edgecolor='black')
        ax_mean.set_title(f"{title_prefix} {cell_id or i+1}: Mean Intensity")
        ax_mean.set_xlabel("Normalized Intensity (Bin Center)")
        ax_mean.set_ylabel("Mean Raw Intensity in Bin")
        ax_mean.set_xticks(bin_centers)
        ax_mean.tick_params(axis='x', rotation=45)

    # Hide unused subplots
    for j in range(n_cells * 2, len(axes)):
        fig.delaxes(axes[j])
        
    fig.tight_layout()
    plt.show()


def plot_bin_colored_image(
    results: List[Dict[str, Any]],
    cell_idx: int = 0,
    raw_cmap: str = 'Greys_r', # Default set to reverse greys for dark background
    bin_cmap: str = 'viridis',
    raw_transform: str = 'none', # 'none', 'log', 'minmax', 'quantile'
    bg_color: str = 'black',
    title_prefix: str = "Cell"
):
    """
    Displays the raw image (transformed) and the bin-colored image side-by-side
    for a single cell.
    """
    if cell_idx >= len(results):
        print(f"Cell index {cell_idx} out of range (0-{len(results)-1}).")
        return
        
    res = results[cell_idx]
    gray = res['gray']
    bin_map = res['bin_map']
    n_bins = res['n_bins']
    cell_id = res['cell_id']
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    
    # --- 1. Raw Image Transformation ---
    raw_img = gray.copy()
    
    if raw_transform == 'log':
        raw_img = np.log1p(raw_img)
        raw_img /= np.max(raw_img) # Rescale to [0, 1]
    elif raw_transform == 'minmax':
        raw_img = (raw_img - raw_img.min()) / (raw_img.max() - raw_img.min() + 1e-8)
    elif raw_transform == 'quantile':
        # Use a robust quantile-based scaling (similar to histogram equalization effects)
        p2, p98 = np.percentile(raw_img, (2, 98))
        raw_img = np.clip(raw_img, p2, p98)
        raw_img = (raw_img - raw_img.min()) / (raw_img.max() - raw_img.min() + 1e-8)
        
    # Apply background color to the raw image subplot
    if bg_color == 'black':
        plt.rcParams['figure.facecolor'] = 'black'
        plt.rcParams['axes.facecolor'] = 'black'
        raw_cmap = raw_cmap if 's_r' not in raw_cmap else raw_cmap.replace('_r', '')
    else:
        plt.rcParams['figure.facecolor'] = 'white'
        plt.rcParams['axes.facecolor'] = 'white'
        
    # --- Plot Raw Image ---
    ax1 = axes[0]
    # Set background color for regions outside the image data area if needed
    ax1.imshow(raw_img, cmap=raw_cmap)
    ax1.set_title(f"{title_prefix} {cell_id or cell_idx}: Raw ({raw_transform})", 
                  color='white' if bg_color == 'black' else 'black')
    ax1.axis('off')

    # --- 2. Bin Colored Image ---
    ax2 = axes[1]
    
    # Create a custom colormap for the bins, plus one color for the background (-1)
    cmap = plt.get_cmap(bin_cmap, n_bins)
    
    # Ensure -1 (background) is transparent or a distinct background color
    # We use a ListedColormap with an extra entry for the background
    new_colors = cmap(np.linspace(0, 1, n_bins))
    
    # Use the bg_color for the excluded background (-1) in the bin map
    if bg_color == 'black':
        bg_rgb = (0., 0., 0., 1.) # Black background
    else:
        bg_rgb = (1., 1., 1., 1.) # White background

    # Combine background color (for index -1) with bin colors (for indices 0 to n_bins-1)
    # The new colormap will map index 0 to new_colors[0], ..., index n_bins-1 to new_colors[n_bins-1]
    # We need to manually set the color for the background index (-1).
    
    # Using pcolormesh is often better for discrete data, but imshow is simpler for fixed data
    
    # We can use boundary normalization to map the indices exactly
    bounds = np.arange(n_bins + 1) - 0.5
    norm = mcolors.BoundaryNorm(bounds, n_bins)
    
    # We create a masked array for plotting. Bin map already has -1 outside the cell.
    masked_bin_map = np.ma.masked_where(bin_map == -1, bin_map)
    
    # Plot the masked bin map. We rely on the axes background color for the -1 regions.
    ax2.imshow(masked_bin_map, cmap=cmap, norm=norm, interpolation='nearest')

    # Add colorbar, showing bin index
    cbar = fig.colorbar(ax2.images[0], ax=ax2, ticks=np.arange(n_bins))
    cbar.set_label('Bin Index')

    ax2.set_title(f"{title_prefix} {cell_id or cell_idx}: Bin Assignment",
                  color='white' if bg_color == 'black' else 'black')
    ax2.axis('off')
    
    fig.tight_layout()
    plt.show()
    # 


# =========================================================================
# 3. SQL EXPORT FUNCTION
# =========================================================================
def export_to_sqlite_simple(
    results: List[Dict[str, Any]],
    image_path: Union[str, Path],  
    db_path: str = "cell_analysis.db",
    table_name: str = "cells"
):
    """
    Exports cell metadata and bin distribution data to a SQLite database.
    Uses a **COMPOSITE PRIMARY KEY** (source_image_path, cell_id) to allow 
    safe iteration over multiple source images without data loss.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    image_path_str = str(Path(image_path).resolve())
    timestamp = datetime.now().isoformat()
    
    data_rows = []
    max_bins = 0

    # 1. Prepare data rows
    for res in results:
        max_bins = max(max_bins, res['n_bins'])
        
        row = {
            'timestamp': timestamp,
            'source_image_path': image_path_str,
            'cell_id': res['cell_id'],
            'total_pixels_in_cell': res['valid_pixel_count'],
            'n_bins': res['n_bins'],
            'normalize_method': res['normalize'],
            'normalization_factor': res['norm_factor_mean_median'] if res['normalize'] in ('mean', 'median') else None,
            'norm_min_q': res['norm_min_q'],    
            'norm_max_q': res['norm_max_q'],    
            'bin_method': res['bin_method'],
            'intensity_range_min': res['intensity_range_min'],
            'intensity_range_max': res['intensity_range_max'],
            'mean_intensity_cell': res['mean_intensity_cell'],
            'median_intensity_cell': res['median_intensity_cell'],
            'bin_edges_json': json.dumps(res['bin_edges'].tolist())
        }
        
        for i in range(res['n_bins']):
            row[f'bin{i}_count'] = int(res['pixel_counts'][i])
            row[f'bin{i}_mean_int'] = float(res['mean_intensities'][i])
            
        data_rows.append(row)

    if not data_rows:
        conn.close()
        return

    # 2. Define table schema
    standard_columns = {
        'timestamp': 'TEXT',
        'source_image_path': 'TEXT',
        'cell_id': 'TEXT', # No longer PRIMARY KEY alone
        'total_pixels_in_cell': 'INTEGER',
        'n_bins': 'INTEGER',
        'normalize_method': 'TEXT',
        'normalization_factor': 'REAL', 
        'norm_min_q': 'REAL',       
        'norm_max_q': 'REAL',       
        'bin_method': 'TEXT',
        'intensity_range_min': 'REAL',
        'intensity_range_max': 'REAL',
        'mean_intensity_cell': 'REAL',
        'median_intensity_cell': 'REAL',
        'bin_edges_json': 'TEXT'
    }

    bin_columns = {}
    for i in range(max_bins):
        bin_columns[f'bin{i}_count'] = 'INTEGER'
        bin_columns[f'bin{i}_mean_int'] = 'REAL'

    all_columns = {**standard_columns, **bin_columns}
    column_defs = ", ".join([f"{name} {dtype}" for name, dtype in all_columns.items()])
    
    # Define the composite primary key
    column_defs += f", PRIMARY KEY (source_image_path, cell_id)"
    
    create_table_sql = f"CREATE TABLE IF NOT EXISTS {table_name} ({column_defs})"
    cursor.execute(create_table_sql)
    
    # 3. Insert data
    insert_cols = list(all_columns.keys())
    placeholders = ", ".join(["?"] * len(insert_cols))
    
    # The ON CONFLICT clause must now reference the composite primary key
    update_fields = ", ".join([f"{col}=excluded.{col}" for col in insert_cols if col not in ['source_image_path', 'cell_id']])
    
    insert_sql = f"INSERT INTO {table_name} ({', '.join(insert_cols)}) VALUES ({placeholders}) ON CONFLICT(source_image_path, cell_id) DO UPDATE SET {update_fields}"
    
    insert_values = []
    for row in data_rows:
        values = []
        for col in insert_cols:
            if col in bin_columns and col not in row:
                values.append(None)
            else:
                values.append(row[col])
        insert_values.append(values)
        
    cursor.executemany(insert_sql, insert_values)
    
    conn.commit()
    conn.close()
    print(f"Successfully exported {len(data_rows)} cells to {db_path} table '{table_name}'.")

# =========================================================================
# 4. EXAMPLE USAGE: Integration Workflow
# =========================================================================


if __name__ == "__main__":
    
    # image_root_dir = Path("/mnt/d/LGLab/Project_NMN/DAPI_analysis_2025-12-15/IMR90/IMR90 NMN SIR DNA_SIM-Measurement")
    # image_root_dir = Path("/mnt/d/LGLab/Project_NMN/DAPI_analysis_2025-12-15/MEF/MEF Scramble KO Control_3D_1_SIM-Measurement")
    # image_root_dir = Path("/mnt/d/LGLab/Project_NMN/DAPI_analysis_2025-12-15/MEF/MEF Scramble KO Control_3D_2_SIM-Measurement")
    image_root_dir = Path("/mnt/d/LGLab/Project_NMN/DAPI_analysis_2025-12-15/MEF/MEF Scramble KO NMN_3D_1_SIM-Measurement")

    dataset = images_to_dataset(
        image_root_dir,
        subset_pattern=".*_z[123][13579]_ORG--t",
        image_subdir="Images", image_suffix=".tif",
        image_extractor="(?P<celltype>.*) (?P<treat>.*) SIR DNA_SIM²_z(?P<stack>.*)_ORG--t(?P<field>\\d+)_r\\d+c\\d+_?(?P<self_generated>.*)"
        # image_extractor="(?P<celltype>.*) Scramble KO (?P<treat>.*)_3D.*_SIM²_z(?P<stack>.*)_ORG--t(?P<field>\\d+)_r\\d+c\\d+_?(?P<self_generated>.*)"
        )
    df = dataset['df']
    
    for idx, row_data in df.iterrows():
    # for idx, row_data in df.iloc[0:2].iterrows():
        # print(row_data['ch0'])
        cell_data_list = crop_cells(
            mask=row_data['directory'] / row_data['cell'], # Mask file path
            imgs=row_data['directory'] / row_data['ch0'], # List of image channel paths
            cell_ids=None, # Extract all cells
            clip_mask=False, 
            pad_square=True, 
            rotate_horizontal=False)
        
        result = intensity_bins_analysis(
            cell_data_list, 
            n_bins=7,  
            normalize='quantile', 
            # bin_method='linear', 
            intensity_range=None,
            log_transform=True,
            bin_method='manual',
            manual_edges=[0.00, 0.2, 0.5, 0.7, 0.85, 0.92, 0.95, 1])
        export_to_sqlite_simple(result, row_data['ch0'], 'cell_analysis_7bins_quantile_manual.db', table_name='cells')
        # plot_intensity_bins(result)
        
        # plt.imshow(cell_data_list[1]['cell_img'])
        # plt.show()
        
    for _, row_data in df.iloc[21:31].iterrows():
        cell_data_list = crop_cells(
            mask=row_data['directory'] / row_data['cell'], # Mask file path
            imgs=row_data['directory'] / row_data['ch0'], # List of image channel paths
            cell_ids=None, # Extract all cells
            clip_mask=True, 
            pad_square=True, 
            rotate_horizontal=True)
        result = intensity_bins_analysis(
            cell_data_list, 
            n_bins=7,  
            normalize='quantile', 
            # bin_method='linear', 
            intensity_range=None,
            log_transform=True,
            bin_method='manual',
            manual_edges=[0.00, 0.2, 0.5, 0.7, 0.85, 0.92, 0.95, 1])
        for idx, _ in enumerate(cell_data_list):
            plot_bin_colored_image(
                result, 
                cell_idx=idx, # Plot the first cell
                raw_cmap='Blues', 
                bin_cmap='plasma', 
                raw_transform='log', # Transform raw image for visualization
                bg_color='white')
        
        
        
        # result = intensity_bins_analysis(
        #     cell_data_list, 
        #     n_bins=7,  
        #     normalize=None, 
        #     bin_method='logspace', 
        #     intensity_range=None)
        # export_to_sqlite_simple(result, row_data['ch0'], 'cell_analysis_none_7bins_log.db', table_name='cells')
        
        # result = intensity_bins_analysis(
        #     cell_data_list, 
        #     n_bins=7,  
        #     normalize='quantile', 
        #     bin_method='linear', 
        #     intensity_range=None)
        # export_to_sqlite_simple(result, row_data['ch0'], 'cell_analysis_quantile_7bins_linear.db', table_name='cells')
        
        # result = intensity_bins_analysis(
        #     cell_data_list, 
        #     n_bins=7,  
        #     normalize='quantile', 
        #     bin_method='logspace', 
        #     intensity_range=None)
        # export_to_sqlite_simple(result, row_data['ch0'], 'cell_analysis_quantile_7bins_log.db', table_name='cells')
        


# row_data = df.iloc[0] # Select the first image/mask pair
# # 2. Extract and crop individual cells from a single image/mask pair
# print(f"Loading and cropping cells from: {row_data['ch0']}")
# cell_data_list = crop_cells(
#     mask=row_data['directory'] / row_data['cell'], # Mask file path
#     imgs=[row_data['directory'] / row_data['ch0']], # List of image channel paths
#     cell_ids=None, # Extract all cells
#     clip_mask=True, 
#     pad_square=False, 
#     rotate_horizontal=False)
# result_linear = intensity_bins_analysis(
#     cell_data_list, 
#     n_bins=7,  
#     normalize='quantile', 
#     bin_method='linear', 
#     intensity_range=(0.0, 1.0))
# plot_intensity_bins(result_linear, title_prefix="Cell (Linear/Mean Norm)")
# plot_bin_colored_image(
#     result_linear, 
#     cell_idx=2, # Plot the first cell
#     raw_cmap='gray', 
#     bin_cmap='tab10', 
#     raw_transform=None, # Transform raw image for visualization
#     bg_color='white') 
# export_to_sqlite_simple(result_linear, row_data['ch0'], 'cell_analysis_quantile_bins.db', table_name='cells_linear')
