
import streamlit as st
import pandas as pd
import numpy as np
import os
from pathlib import Path
from skimage import exposure

from image_helper import find_measurement_dirs, images_to_dataset, crop_cells




# --------------------------------------------------------------------------- #
# Functions
# --------------------------------------------------------------------------- #
def normalize_image(img: np.ndarray, method: str):
    """
    img: numpy array
    method: one of []
    """
    if method == "None":
        return img
    elif method == "Percentile":
        p1, p2 = np.percentile(img, (0.1, 99.9))
        return exposure.rescale_intensity(img, in_range=(p1, p2), out_range=(0, 0.9999))
    elif method == "Equalize Histogram":
        return exposure.equalize_hist(img)
    elif method == "CLAHE":
        return exposure.equalize_adapthist(img, clip_limit=clip_limit)
    elif method == "Gamma":
        return exposure.adjust_gamma(img, gamma=gamma, gain=1)
    return img

def normalize_image_by_group(img: np.ndarray,
                             method: str,
                             by: str = "None") -> np.ndarray:
    """
    Normalize 4D image (B, H, W, C) along specified axes.
    
    Parameters
    ----------
    img : np.ndarray
        Shape (B, H, W, C)
    method : str
        Normalization method
    by : str
        "None"  -> no grouping, normalize whole volume at once
        "B"     -> normalize each batch item independently (HWC for each B)
        "C"     -> normalize each channel independently (across B,H,W)
        "BC"    -> normalize each (H,W) slice per batch and channel
    """
    if img.ndim != 4:
        raise ValueError("Input must be 4D array (B, H, W, C)")

    if method == "None" or by == "None":
        return normalize_image(img, method)

    img = img.astype(np.float32)  # important for most methods

    if by == "B":
        # Apply normalization independently to each batch element (axis 0)
        normalized_batches = []
        for i in range(img.shape[0]):
            normalized_batches.append(normalize_image(img[i], method))
        return np.stack(normalized_batches, axis=0)
    
        # Or more efficiently with vectorized functions (recommended):
        # return np.array([normalize_image(img[i], method) for i in range(img.shape[0])])

    elif by == "C":
        # Normalize each channel independently
        normalized_channels = []
        for c in range(img.shape[3]):
            channel = img[..., c]
            normalized_channels.append(normalize_image(channel, method)[..., np.newaxis])
        return np.concatenate(normalized_channels, axis=-1)

    elif by == "BC":
        # Normalize each (H,W) plane independently across B and C
        # Reshape to (B*C, H, W), normalize, then reshape back
        B, H, W, C = img.shape
        img_reshaped = img.transpose(0, 3, 1, 2).reshape(B * C, H, W)  # (B*C, H, W)
        normalized_reshaped = normalize_image(img_reshaped, method)
        return normalized_reshaped.reshape(B, C, H, W).transpose(0, 2, 3, 1)

    else:
        raise ValueError(f"Invalid 'by' parameter: {by}. Choose from 'B', 'C', 'BC', 'None'")




# --------------------------------------------------------------------------- #
# Sidebar: Inputs
# --------------------------------------------------------------------------- #
# Set page config
st.title("Single Cell Image Viewer")
st.set_page_config(page_title="SCIV", page_icon=":rocket:", layout="wide")

st.sidebar.header("Input Data")
root_dir = st.sidebar.text_input("Root Directory (containing Measurement folders)",
    value="/mnt/d/TP53/HCS", 
    help="Path to folder containing folders like '*__2025-11-17T15_50_02-Measurement 1'")
result_file = st.sidebar.file_uploader("Upload result_df.xlsx or .csv", type=["xlsx", "csv"])

st.sidebar.header("Cell Cropping")
clip_mask = st.sidebar.checkbox("Clip mask", value=False)
pad_to_square = st.sidebar.checkbox("Pad to square", value=False)
rotate_to_horizontal = st.sidebar.checkbox("Rotate horizontally", value=False)
target_crop_size = st.sidebar.number_input("Resize (px)", min_value=64, max_value=512, value=128, step=64)

st.sidebar.header("Cell Normalization")
norm_method = st.sidebar.selectbox("Normalization", ["Percentile", "Equalize Histogram", "CLAHE", "Gamma","None"])
if norm_method == "CLAHE":
    clip_limit = st.sidebar.number_input("Clip limit", 0.0, 1.0, 0.02, step=0.005)
if norm_method == "Gamma":
    gamma = st.sidebar.number_input("Gamma", 0.0, 5.0, 0.4, step=0.1)
norm_axis = st.sidebar.radio("Norm axis", ["B", "C", "BC"], horizontal=True)

st.sidebar.header("Display Option")
show_filename = st.sidebar.checkbox("Show filename", value=False)
n_fields_per_row = st.sidebar.number_input("N fields per row", min_value=5, max_value=50, value=10, step=5)





# --------------------------------------------------------------------------- #
# Load and Cache Data
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner="Building dataloader_all from all measurements...")
def build_dataloader_all(root_path: str):
    root = Path(root_path)
    if not root.exists():
        st.error(f"Root directory not found: {root}")
        return []
    measurement_dirs = find_measurement_dirs(root)

    if not measurement_dirs:
        st.error("No measurement directories found.")
        return None
    
    dfs = []
    for meas_dir in measurement_dirs:
        print(f"[DEBUG] Processing measurement: {meas_dir}")

        cp_csv = meas_dir / "cp_dataloader.csv"
        if cp_csv.exists():
            print(f"Loading from CSV: {cp_csv}")
            df = pd.read_csv(cp_csv)
        else:
            result = images_to_dataset(
                meas_dir,
                cellprofiler_style=True,
                mask_suffix=".png",
                image_suffix=".tiff"
            )
            if result is not None:
                df = result["df"]

        if not df.empty:
            dfs.append(df)
            print(f"[DEBUG] Added {len(df)} sites from {meas_dir.name}")
    
    if not dfs:
        st.error("No image data found.")
        return None
    
    dataloader_all = pd.concat(dfs, ignore_index=True)
    print(f"[DEBUG] Total sites loaded: {len(dataloader_all)}")
    return dataloader_all

@st.cache_data(show_spinner="Loading result table...")
def load_result_df(file):
    df = pd.read_csv(file) if file.name.endswith(".csv") else pd.read_excel(file)
    #     df = pd.read_csv(file)
    # else:
    #     df = pd.read_excel(file)
    print(f"[DEBUG] Result table loaded: {df.shape}")
    return df




# --------------------------------------------------------------------------- #
# Pre-process Data
# --------------------------------------------------------------------------- #
if not root_dir or not result_file:
    st.info("Please provide root directory and upload result file.")
    st.stop()

dataloader_all = build_dataloader_all(root_dir)
result_df = load_result_df(result_file)

if dataloader_all is None or result_df is None:
    st.stop()

# Normalize directory names for matching
def normalize_directory_name(path_str: str) -> str:
    """Extract short name like '2025-11-08_p53_mutants_dox_22' from full path"""
    path = Path(path_str)
    name = path.parent.name
    # Handle both "__2025-..." and direct short names
    if "__" in name:
        return name.split("__")[0]
    return name

dataloader_all["Metadata_directory"] = dataloader_all["Metadata_directory"].apply(normalize_directory_name)
# result_df["Metadata_directory"] = result_df["Metadata_directory"].apply(normalize_directory_name)
result_df["Metadata_directory"] = result_df["Metadata_directory"].astype(str)

# Find common metadata columns
metadata_cols = [c for c in result_df.columns if c.startswith("Metadata_")]
numeric_cols = [c for c in result_df.columns if not c.startswith("Metadata_")]




# --------------------------------------------------------------------------- #
# Dynamic filters/sorting/grouping
# --------------------------------------------------------------------------- #
filtered_result = result_df.copy()

with st.expander("Apply Filters", expanded=True):
    slider_n_per_row = 4
    cols = st.columns(slider_n_per_row, gap="small")
    filter_widgets = {}
    for i, col in enumerate(metadata_cols):
        with cols[i % slider_n_per_row]:
            options = ["All"] + sorted(result_df[col].dropna().unique().tolist())
            selected = st.multiselect(f"{col}", options, default=["All"])
            if "All" not in selected and selected:
                filtered_result = filtered_result[filtered_result[col].isin(selected)]

    for i, col in enumerate(numeric_cols):
        with cols[i % slider_n_per_row]:
            mn, mx = float(result_df[col].min()), float(result_df[col].max())
            vmin, vmax = st.slider(f"{col}", mn, mx, (mn, mx), step=(mx-mn)/20)
            filtered_result = filtered_result[(filtered_result[col] >= vmin) & (filtered_result[col] <= vmax)]

if len(filtered_result) == 0:
    st.warning("No cells match the filters.")
    st.stop()

st.success(f"Filtered to {len(filtered_result)} cells from result table.")

# get common keys
common_keys = ["Metadata_directory", "Metadata_well", "Metadata_field"]
for c in ["Metadata_stack", "Metadata_timepoint"]:
    if c in filtered_result.columns: 
        common_keys += [c]

# Merge on common keys
filtered_sites = pd.merge(filtered_result, dataloader_all, on=common_keys, how="inner")

if len(filtered_sites) == 0:
    st.error("No matching image sites found. Check directory naming.")
    st.stop()
st.info(f"Found {len(filtered_sites)} site-cell combinations with images.")


# Apply sorting
display_df = filtered_sites.copy()

# set input control
max_cells_per_image = st.sidebar.number_input("Max cells per image", min_value=5, max_value=100, value=20, step=5)
sort_col = st.sidebar.selectbox("Sort cells by", numeric_cols)
sort_order = st.sidebar.radio("Sort order", ["Increase", "Random", "Decrease", "As is"], horizontal=True)

if sort_order == "Increase":
    display_df = display_df.sort_values(sort_col, ascending=True)
elif sort_order == "Decrease":
    display_df = display_df.sort_values(sort_col, ascending=False)
elif sort_order == "Random":
    display_df = display_df.sample(frac=1, random_state=42)
else:
    pass
print("sorted", display_df)

# Group by directory -> well
group_cols = st.sidebar.multiselect(f"Grouping by", metadata_cols, default=["Metadata_directory"])
all_channels = [c for c in display_df.columns if c.startswith("Image_FileName_ch")]
show_channels = st.sidebar.multiselect(f"Show Channels", all_channels, default=all_channels[0])


st.sidebar.header("\n---- Set Column ID ----")
mask_file_column = st.sidebar.selectbox("Mask Column", [c for c in dataloader_all.columns if c.startswith("Image_ObjectsFileName_")])
cell_id_column = st.sidebar.selectbox(" Cell ID Column", numeric_cols, 2) # cell_id_column = 'Parent_mask_cell'





# --------------------------------------------------------------------------- #
# Display
# --------------------------------------------------------------------------- #
for group_key, group_df in display_df.groupby(group_cols, dropna=False):
    # print(group_key, group_df)
    group_key = [group_key] if isinstance(group_key, str) else group_key
    group_name = " → ".join([str(k) for k in group_key])

    with st.expander(f"{group_name} ({len(group_df)} cells)", expanded=True):

        cropped_all = []
        count = 0
        cols = st.columns(n_fields_per_row, gap="small")

        for per_image_key, per_image_df in group_df.groupby(show_channels, dropna=False):
            
            try:
                parent_dir = Path(per_image_df[show_channels[0].replace("Image_FileName_", "Image_PathName_")].values.tolist()[0])
                mask_path = parent_dir / per_image_df[mask_file_column].values.tolist()[0]
                image_channel_paths = [parent_dir / c for c in per_image_df[show_channels].iloc[0].tolist()]
                cell_ids = per_image_df[cell_id_column].values.tolist()

                if len(cell_ids) > max_cells_per_image: 
                    cell_ids = cell_ids[:max_cells_per_image]

                if not os.path.exists(mask_path):
                    print(f"[WARN] Mask missing: {mask_path}")
                    continue
                if any([not os.path.exists(img_path) for img_path in image_channel_paths]):
                    print(f"[WARN] img_paths missing: {image_channel_paths}")
                    continue
                
                cropped = crop_cells(
                    mask_path, image_channel_paths, cell_ids,
                    scale_factor=65535.0,
                    target_size=target_crop_size if target_crop_size > 0 else None,
                    clip_mask=clip_mask,
                    pad_square=pad_to_square,
                    rotate_horizontal=rotate_to_horizontal)
                cropped = [x['img'] for x in cropped] # extract img
                cropped = np.stack(cropped, axis=0)  # list to array
                
                if cropped is None: 
                    continue

                # normalization
                cropped_norm = normalize_image_by_group(cropped, norm_method, by=norm_axis)

                # concatenate cells in one image
                B, H, W, C = cropped_norm.shape
                cropped_norm = cropped_norm.reshape(B * H, W * C)

                # show image
                caption = f"{image_channel_paths[0].name.replace('fk1fl1.tiff','')}" if show_filename else None
                with cols[count % n_fields_per_row]:
                    st.image(cropped_norm, caption=caption, clamp=True, width='stretch')
                count += 1
            except Exception as e:
                print("mask_path: ", mask_path)
                print("image_channel_paths: ", image_channel_paths)
                print("cell_ids: ", cell_ids)
                print("cropped: ", cropped.shape)
                print("cropped_norm: ", cropped_norm.shape)
            finally:
                continue

st.success("Done! Use filters to explore.")


