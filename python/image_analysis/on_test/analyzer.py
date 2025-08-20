#%%
import numpy as np
import pandas as pd
from scipy.stats import skew
from scipy import ndimage
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries
from skimage import morphology, transform, measure
from skimage.morphology import disk
from typing import Dict, List, Tuple, Optional, Union, Any
import matplotlib.pyplot as plt
import cv2


class BaseImageAnalyzer:
    """
    Base class for multi-channel intensity analysis with optional cell-level analysis.
    
    Provides common functionality for image validation, preprocessing, and data management.
    """
    
    def __init__(self, channel_names: Optional[List[str]] = None):
        """
        Initialize BaseImageAnalyzer.
        
        Parameters:
        -----------
        channel_names : list, optional
            Default channel names to use. Can be overridden in measure()
        """
        self.results = None
        self.intensity_data = None
        self.mask_data = None
        self.channel_names = channel_names
    
    def _validate_inputs(self, intensity_array: np.ndarray, mask_array: Optional[np.ndarray] = None) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Validate and standardize input arrays.
        
        Parameters:
        -----------
        intensity_array : ndarray
            A 2D (H, W) or 3D (C, H, W) numpy array
        mask_array : ndarray, optional
            Shape (H, W) with integer labels for different cells/objects
            
        Returns:
        --------
        Tuple[np.ndarray, Optional[np.ndarray]]
            Validated intensity and mask arrays
        """
        if intensity_array.ndim not in [2, 3]:
            raise ValueError("intensity_array must be a 2D (H, W) or 3D (C, H, W) array.")

        # If image is 2D, expand to 3D for consistent processing
        if intensity_array.ndim == 2:
            intensity_array = intensity_array[np.newaxis, :, :]

        n_channels, height, width = intensity_array.shape
        
        if mask_array is not None:
            if mask_array.ndim != 2:
                raise ValueError("mask_array must be a 2D (H, W) array.")
            if mask_array.shape != (height, width):
                raise ValueError("Mask array must have shape (H, W) matching intensity array")
        
        return intensity_array, mask_array
    
    def _setup_channel_names(self, intensity_array: np.ndarray, channel_names: Optional[List[str]] = None) -> List[str]:
        """
        Setup and validate channel names.
        
        Parameters:
        -----------
        intensity_array : ndarray
            Shape (C, H, W)
        channel_names : list, optional
            Names for each channel
            
        Returns:
        --------
        List[str]
            Validated channel names
        """
        n_channels = intensity_array.shape[0]
        
        # Priority: parameter > instance attribute > auto-generate
        if channel_names is not None:
            if len(channel_names) != n_channels:
                raise ValueError("Number of channel names must match number of channels")
            return channel_names
        elif self.channel_names is not None:
            if len(self.channel_names) != n_channels:
                raise ValueError("Instance channel names don't match number of channels")
            return self.channel_names
        else:
            return [f"ch{i}" for i in range(n_channels)]
    
    def _store_data(self, intensity_array: np.ndarray, mask_array: Optional[np.ndarray], channel_names: List[str]):
        """Store data for later access (e.g., plotting)."""
        self.intensity_data = intensity_array.copy()
        self.mask_data = mask_array.copy() if mask_array is not None else None
        self.channel_names = channel_names
    
    def _reorder_columns(self, df: pd.DataFrame, channel_names: List[str]) -> pd.DataFrame:
        """
        Reorder dataframe columns so channel-specific columns come last.
        
        Parameters:
        -----------
        df : pd.DataFrame
            Input dataframe
        channel_names : List[str]
            List of channel names
            
        Returns:
        --------
        pd.DataFrame
            Dataframe with reordered columns
        """
        if df.empty:
            return df
        
        # Identify channel-specific columns
        channel_columns = []
        non_channel_columns = []
        
        for col in df.columns:
            is_channel_col = any(ch_name in col for ch_name in channel_names)
            if is_channel_col:
                channel_columns.append(col)
            else:
                non_channel_columns.append(col)
        
        # Sort channel columns by channel name order, then by metric name
        channel_columns_sorted = []
        for ch_name in channel_names:
            ch_cols = [col for col in channel_columns if ch_name in col]
            ch_cols.sort()  # Sort by metric name within channel
            channel_columns_sorted.extend(ch_cols)
        
        # Add any remaining channel columns not caught above
        remaining_ch_cols = [col for col in channel_columns if col not in channel_columns_sorted]
        channel_columns_sorted.extend(sorted(remaining_ch_cols))
        
        # Reorder: non-channel columns first, then channel columns
        new_column_order = non_channel_columns + channel_columns_sorted
        
        return df[new_column_order]
    
    def measure(self, intensity_array: np.ndarray, mask_array: Optional[np.ndarray] = None, 
                channel_names: Optional[List[str]] = None, **kwargs) -> Dict[str, Any]:
        """
        Abstract method to be implemented by subclasses.
        
        Returns:
        --------
        Dict[str, Any]
            Results dictionary with measurements
        """
        raise NotImplementedError("Subclasses must implement the measure method")


class IntensityAnalyzer(BaseImageAnalyzer):
    """
    Analyzer for basic intensity statistics with optional cell-level analysis.
    """
    
    def __init__(self, channel_names: Optional[List[str]] = None, percentiles: List[float] = [0.1, 1, 99, 99.9]):
        """
        Initialize IntensityAnalyzer.
        
        Parameters:
        -----------
        channel_names : list, optional
            Default channel names to use
        percentiles : list, optional
            Default percentiles to calculate
        """
        super().__init__(channel_names)
        self.percentiles = percentiles
    
    def measure(self, intensity_array: np.ndarray, mask_array: Optional[np.ndarray] = None, 
                channel_names: Optional[List[str]] = None, 
                percentiles: Optional[List[float]] = None,
                measure_global: bool = True
                ) -> Dict[str, pd.DataFrame]:
        """
        Measure intensity statistics from multi-channel data.
        
        Parameters:
        -----------
        intensity_array : ndarray
            A 2D (H, W) or 3D (C, H, W) numpy array
        mask_array : ndarray, optional
            Shape (H, W) with integer labels for different cells/objects
        channel_names : list, optional
            Names for each channel
        percentiles : list, optional
            List of percentiles to calculate
        measure_global: bool, 
            measure at per-image level
            
        Returns:
        --------
        Dict[str, pd.DataFrame]
            Dictionary with 'global' and 'objects' keys
        """
        # Validate inputs
        intensity_array, mask_array = self._validate_inputs(intensity_array, mask_array)
        channel_names = self._setup_channel_names(intensity_array, channel_names)
        if channel_names is None:
            return None

        percentiles_to_use = percentiles if percentiles is not None else self.percentiles
        if not measure_global and mask_array is None:
            raise ValueError("measure_global and mask_array musted provied at least one")

        # Store data for plotting
        self._store_data(intensity_array, mask_array, channel_names)
        
        results = {}
        # Measure global intensity
        if measure_global:
            global_stats = self._measure_global_intensity(intensity_array, percentiles_to_use, channel_names)
            # global_df = self._reorder_columns(global_df, channel_names)
            results['global'] = global_stats
        
        # Measure object-level intensities if mask provided
        if mask_array is not None:
            object_stats = self._measure_object_intensity(intensity_array, mask_array, channel_names)
            object_stats = self._reorder_columns(object_stats, channel_names)
            results['objects'] = object_stats
        
        self.results = results
        return self.results
    
    def _measure_global_intensity(self, intensity_array: np.ndarray, percentiles: List[float], channel_names: List[str]) -> Dict[str, float]:
        """Measure global intensity statistics for all channels."""
        stats = {}
        
        for ch_idx, ch_name in enumerate(channel_names):
            ch_data = intensity_array[ch_idx]
            
            # Basic statistics
            stats.update({
                f"mean_{ch_name}": np.mean(ch_data),
                f"std_{ch_name}": np.std(ch_data),
                f"min_{ch_name}": np.min(ch_data),
                f"max_{ch_name}": np.max(ch_data),
                f"sum_{ch_name}": np.sum(ch_data),
                f"median_{ch_name}": np.median(ch_data)
            })
            
            # Percentile statistics
            if percentiles:
                percentile_values = np.percentile(ch_data, percentiles)
                for p, p_val in zip(percentiles, percentile_values):
                    stats[f"p{int(p*10)/10}_{ch_name}"] = p_val
        
        return pd.DataFrame([stats])
    
    def _measure_object_intensity(self, intensity_array: np.ndarray, mask_array: np.ndarray, channel_names: List[str]) -> pd.DataFrame:
        """Measure object-level intensity statistics using regionprops."""
        results_list = []
        
        # Measure all objects at once for each channel using regionprops
        for ch_idx, ch_name in enumerate(channel_names):
            ch_data = intensity_array[ch_idx]
            
            # Get regionprops for this channel
            props = regionprops(mask_array, intensity_image=ch_data, extra_properties=[_skew_intensity])
            
            # Store channel-specific measurements
            if ch_idx == 0:
                # First channel: initialize results with morphological properties
                for prop in props:
                    result_row = {
                        'object_id': prop.label,
                        'centroid_y': prop.centroid[0],
                        'centroid_x': prop.centroid[1],
                        'bbox_min_y': prop.bbox[0],
                        'bbox_min_x': prop.bbox[1],
                        'bbox_max_y': prop.bbox[2],
                        'bbox_max_x': prop.bbox[3],
                        'area': prop.area,
                        'diameter': prop.equivalent_diameter,
                        'perimeter': prop.perimeter,
                        'eccentricity': prop.eccentricity,
                        'solidity': prop.solidity,
                        'major_axis_length': prop.major_axis_length,
                        'minor_axis_length': prop.minor_axis_length
                    }
                    
                    # Add intensity measurements for this channel
                    result_row.update({
                        f"mean_{ch_name}": prop.intensity_mean,
                        f"min_{ch_name}": prop.intensity_min,
                        f"max_{ch_name}": prop.intensity_max,
                        f"std_{ch_name}": np.std(ch_data[mask_array == prop.label]),
                        f"sum_{ch_name}": np.sum(ch_data[mask_array == prop.label]),
                        f"median_{ch_name}": np.median(ch_data[mask_array == prop.label]),
                        f"skew_{ch_name}": prop._skew_intensity,
                    })
                    
                    results_list.append(result_row)
            else:
                # Subsequent channels: add intensity measurements to existing results
                for i, prop in enumerate(props):
                    if i < len(results_list):  # Safety check
                        results_list[i].update({
                            f"mean_{ch_name}": prop.intensity_mean,
                            f"min_{ch_name}": prop.intensity_min,
                            f"max_{ch_name}": prop.intensity_max,
                            f"std_{ch_name}": np.std(ch_data[mask_array == prop.label]),
                            f"sum_{ch_name}": np.sum(ch_data[mask_array == prop.label]),
                            f"median_{ch_name}": np.median(ch_data[mask_array == prop.label]),
                            f"skew_{ch_name}": prop._skew_intensity,
                        })
        
        return pd.DataFrame(results_list)
    
    def plot_overlay(self, channel_index: int = 0, figsize: Tuple[int, int] = (15, 6), 
                    show_labels: bool = True, show_centroids: bool = True, 
                    colormap: str = 'viridis', alpha_mask: float = 0.3):
        """
        Plot intensity image with overlay of measurement results.
        
        Parameters:
        -----------
        channel_index : int
            Which channel to display
        figsize : tuple
            Figure size
        show_labels : bool
            Whether to show object ID labels
        show_centroids : bool
            Whether to show centroids
        colormap : str
            Colormap for intensity image
        alpha_mask : float
            Transparency of mask overlay
        """
        if self.intensity_data is None:
            raise ValueError("No data to plot. Run measure() first.")
        
        if channel_index >= len(self.channel_names):
            raise ValueError(f"Channel index {channel_index} out of range")
        
        # Determine number of subplots
        n_plots = 3 if self.mask_data is not None else 2
        fig, axes = plt.subplots(1, n_plots, figsize=figsize)
        if n_plots == 2:
            axes = [axes[0], axes[1]]
        
        # Get the channel data
        channel_data = self.intensity_data[channel_index]
        channel_name = self.channel_names[channel_index]
        
        # Plot 1: Original intensity image
        im1 = axes[0].imshow(channel_data, cmap=colormap)
        axes[0].set_title(f'{channel_name} - Original')
        axes[0].axis('off')
        plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
        
        # Plot 2: With measurements overlay
        im2 = axes[1].imshow(channel_data, cmap=colormap)
        title = f'{channel_name} - With Measurements'
        if self.mask_data is not None and 'objects' in self.results:
            title += f' ({len(self.results["objects"])} objects)'
        axes[1].set_title(title)
        axes[1].axis('off')
        
        if self.mask_data is not None:
            # Create colored mask overlay
            unique_labels = np.unique(self.mask_data)
            if 0 in unique_labels:
                unique_labels = unique_labels[unique_labels != 0]
            
            mask_colored = np.zeros_like(self.mask_data, dtype=float)
            for i, label_id in enumerate(unique_labels):
                mask_colored[self.mask_data == label_id] = i + 1
            
            # Overlay mask with transparency
            mask_colored_ma = np.ma.masked_where(mask_colored == 0, mask_colored)
            axes[1].imshow(mask_colored_ma, alpha=alpha_mask, cmap='tab20')
            
            # Add boundaries
            boundaries = find_boundaries(self.mask_data, mode='outer')
            axes[1].contour(boundaries, colors='white', linewidths=1)
            
            # Add labels and centroids
            if 'objects' in self.results:
                for _, row in self.results['objects'].iterrows():
                    if 'centroid_x' in row:
                        x, y = row['centroid_x'], row['centroid_y']
                        
                        if show_centroids:
                            axes[1].plot(x, y, 'r+', markersize=10, markeredgewidth=2)
                        
                        if show_labels:
                            axes[1].text(x, y-8, f"{int(row['object_id'])}", 
                                       ha='center', va='bottom', color='white', 
                                       fontsize=9, fontweight='bold',
                                       bbox=dict(boxstyle='round,pad=0.2', 
                                               facecolor='black', alpha=0.8))
            
            # Plot 3: Mask only (if mask exists)
            if n_plots == 3:
                im3 = axes[2].imshow(self.mask_data, cmap='tab20')
                axes[2].set_title('Object Mask')
                axes[2].axis('off')
                plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)
        
        plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.show()
        
        # Print summary statistics
        self._print_summary(channel_index)
    
    def _print_summary(self, channel_index: int = 0):
        """Print summary statistics."""
        print(f"\n{'='*60}")
        print(f"INTENSITY MEASUREMENT SUMMARY")
        print(f"{'='*60}")
        
        # Global statistics
        if 'global' in self.results:
            print(f"\nGLOBAL IMAGE STATISTICS:")
            print(f"Image dimensions: {self.intensity_data.shape}")
            print(f"Number of channels: {len(self.channel_names)}")
            print(f"Channel names: {', '.join(self.channel_names)}")
            
            global_stats = self.results['global'].iloc[0]
            for ch_name in self.channel_names:
                mean_key = f"mean"
                std_key = f"std"
                if mean_key in global_stats:
                    print(f"{ch_name}: {global_stats[mean_key]:.2f} ± {global_stats[std_key]:.2f}")
                    
                    # Print percentiles if available
                    percentile_keys = [k for k in global_stats.keys() if k.startswith(f"p")]
                    if percentile_keys:
                        percentile_str = ", ".join([f"P{k.split('_p')[1]}: {global_stats[k]:.1f}" 
                                                  for k in sorted(percentile_keys)])
                        print(f"  Percentiles - {percentile_str}")
        
        # Object-level statistics
        if 'objects' in self.results:
            object_df = self.results['objects']
            print(f"\nOBJECT-LEVEL STATISTICS:")
            print(f"Number of objects: {len(object_df)}")
            if 'area' in object_df.columns:
                print(f"Mean area: {object_df['area'].mean():.2f} ± {object_df['area'].std():.2f} pixels")
            
            # Channel-specific object statistics
            channel_name = self.channel_names[channel_index]
            mean_col = f"mean_{channel_name}"
            if mean_col in object_df.columns:
                mean_intensity = object_df[mean_col].mean()
                std_intensity = object_df[mean_col].std()
                print(f"{channel_name} object mean intensity: {mean_intensity:.2f} ± {std_intensity:.2f}")
        
        print(f"{'='*60}")


class RadialDistributionAnalyzer(BaseImageAnalyzer):
    """
    Analyzer for radial distribution of intensities within objects.
    """
    
    def __init__(self, channel_names: Optional[List[str]] = None, num_bins: int = 5, 
                 center_method: str = "object_center", angle_bins: int = 180):
        """
        Initialize RadialDistributionAnalyzer.
        
        Parameters:
        -----------
        channel_names : list, optional
            Default channel names to use
        num_bins : int
            Number of concentric bins
        center_method : str
            Method to determine center: "object_center" or "centroid"
        angle_bins : int
            Number of angular bins for edge distance calculation
        """
        super().__init__(channel_names)
        self.num_bins = num_bins
        self.center_method = center_method
        self.angle_bins = angle_bins
        
        if self.center_method not in ["object_center", "centroid"]:
            raise ValueError("center_method must be 'object_center' or 'centroid'")
    
    def measure(self, intensity_array: np.ndarray, mask_array: np.ndarray, 
                channel_names: Optional[List[str]] = None) -> Dict[str, Union[pd.DataFrame, Dict]]:
        """
        Measure radial distribution of intensities within objects.
        
        Parameters:
        -----------
        intensity_array : ndarray
            A 2D (H, W) or 3D (C, H, W) numpy array
        mask_array : ndarray
            Shape (H, W) with integer labels
        channel_names : list, optional
            Names for each channel
            
        Returns:
        --------
        Dict[str, Union[pd.DataFrame, Dict]]
            Dictionary with 'radialdistribution' and 'metadata' keys
        """
        # Validate inputs
        intensity_array, mask_array = self._validate_inputs(intensity_array, mask_array)
        channel_names = self._setup_channel_names(intensity_array, channel_names)
        if channel_names is None:
            return None

        if mask_array is None:
            raise ValueError("RadialDistributionAnalyzer requires a mask_array")
        
        # Store data for plotting
        self._store_data(intensity_array, mask_array, channel_names)
        
        object_labels = np.unique(mask_array)
        object_labels = object_labels[object_labels != 0]

        wide_results_list = []
        object_metadata = {}

        for obj_id in object_labels:
            obj_mask = (mask_array == obj_id)

            if self.center_method == "centroid":
                props = measure.regionprops(obj_mask.astype(np.uint8))[0]
                center_y, center_x = props.centroid
            else:
                center_y, center_x = self._get_object_center(obj_mask)

            # Get normalized radial positions
            radial_positions, radial_bins_image = self._calculate_shape_aware_radial_positions(
                obj_mask, center_y, center_x
            )
            
            if radial_positions is None or len(radial_positions) == 0:
                continue

            obj_coords = np.where(obj_mask)
            
            object_metadata[int(obj_id)] = {
                'center': (center_y, center_x),
                'radial_bins_image': radial_bins_image
            }

            obj_results = {'object_id': int(obj_id)}
            for channel_idx, channel_name in enumerate(channel_names):
                intensities = intensity_array[channel_idx][obj_coords]
                total_intensity = np.sum(intensities)
                total_pixels = len(intensities)

                pixel_bin_indices = radial_bins_image[obj_mask]

                for bin_num in range(1, self.num_bins + 1):
                    bin_mask_indices = (pixel_bin_indices == bin_num)
                    if not np.any(bin_mask_indices):
                        frac_at_d, mean_frac, radial_cv = 0.0, 0.0, 0.0
                    else:
                        bin_intensities = intensities[bin_mask_indices]
                        bin_intensity_sum = np.sum(bin_intensities)
                        pixel_fraction = len(bin_intensities) / total_pixels
                        frac_at_d = bin_intensity_sum / total_intensity if total_intensity > 0 else 0
                        mean_frac = (frac_at_d / pixel_fraction) if pixel_fraction > 0 else 0
                        mean_intensity = np.mean(bin_intensities)
                        radial_cv = (np.std(bin_intensities) / mean_intensity) if mean_intensity > 0 else 0.0
                    
                    obj_results[f"radialFracAtD_{bin_num}_{channel_name}"] = frac_at_d
                    obj_results[f"radialMeanFrac_{bin_num}_{channel_name}"] = mean_frac
                    obj_results[f"radialCV_{bin_num}_{channel_name}"] = radial_cv
            
            wide_results_list.append(obj_results)
        
        radialdistribution_df = pd.DataFrame(wide_results_list)
        radialdistribution_df = self._reorder_columns(radialdistribution_df, channel_names)
        
        results = {
            'radialdistribution': radialdistribution_df,
            'metadata': object_metadata
        }
        
        self.results = results
        return self.results
    
    def _get_object_center(self, binary_mask: np.ndarray) -> Tuple[float, float]:
        """Find center as point farthest from any edge using distance transform."""
        if binary_mask.dtype == bool:
            binary_mask = binary_mask.astype(np.uint8)
        
        distance_transform = cv2.distanceTransform(binary_mask, cv2.DIST_L2, 5).astype(np.float32)
        center_coords = np.where(distance_transform == distance_transform.max())
        center_y, center_x = center_coords[0].mean(), center_coords[1].mean()
        
        return float(center_y), float(center_x)

    def _calculate_shape_aware_radial_positions(self, obj_mask: np.ndarray, center_y: float, center_x: float) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Calculate radial positions using shape-aware approach."""
        if not np.any(obj_mask):
            return None, None
            
        y, x = np.ogrid[:obj_mask.shape[0], :obj_mask.shape[1]]
        distances_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
        angles = np.arctan2(y - center_y, x - center_x)
        
        mask_coords = np.where(obj_mask)
        if len(mask_coords[0]) == 0:
            return None, None
            
        pixel_distances = distances_from_center[mask_coords]
        pixel_angles = angles[mask_coords]
        
        # Bin angles for edge distance calculation
        angle_edges = np.linspace(-np.pi, np.pi, self.angle_bins + 1)
        binned_angles = np.clip(np.digitize(pixel_angles, angle_edges) - 1, 0, self.angle_bins - 1)
        
        # Find max distance for each angle bin
        max_distances = np.zeros(self.angle_bins)
        for i in range(self.angle_bins):
            angle_mask = binned_angles == i
            if np.any(angle_mask):
                max_distances[i] = np.max(pixel_distances[angle_mask])
        
        # Fill gaps using circular interpolation
        valid_indices = np.nonzero(max_distances)[0]
        if len(valid_indices) > 1:
            invalid_indices = np.where(max_distances == 0)[0]
            if len(invalid_indices) > 0:
                max_distances[invalid_indices] = np.interp(
                    invalid_indices, 
                    valid_indices,
                    max_distances[valid_indices],
                    period=self.angle_bins
                )
        
        edge_distances = max_distances[binned_angles]
        edge_distances = np.maximum(edge_distances, 1e-8)
        
        # Calculate normalized radial positions
        radial_positions = pixel_distances / edge_distances
        radial_positions = np.clip(radial_positions, 0, 1)
        
        # Create bins
        bin_edges = np.linspace(0, 1, self.num_bins + 1)
        pixel_bin_indices = np.digitize(radial_positions, bin_edges, right=True)
        pixel_bin_indices = np.clip(pixel_bin_indices, 1, self.num_bins)
        
        # Store radial bin assignments
        radial_bins_image = np.zeros(obj_mask.shape, dtype=int)
        radial_bins_image[mask_coords] = pixel_bin_indices
        
        return radial_positions, radial_bins_image

    def plot_radial_distribution(self, image: np.ndarray, mask: np.ndarray,
                                channel_name: str, object_id: Optional[int] = None,
                                metric_name: str = 'radialMeanFrac',
                                cmap_name: str = 'viridis', alpha: float = 0.5):
        """
        Plot radial distribution overlay on an image.
        
        Parameters:
        -----------
        image : np.ndarray
            The original 2D image to draw on
        mask : np.ndarray
            The mask array with labeled objects
        channel_name : str
            The name of the channel to plot
        object_id : Optional[int]
            The ID of the object to visualize (None for all objects)
        metric_name : str
            The metric to plot ('radialFracAtD', 'radialMeanFrac', 'radialCV')
        cmap_name : str
            Matplotlib colormap name
        alpha : float
            Transparency of overlay
        """
        if self.results is None:
            raise ValueError("No results to plot. Run measure() first.")
            
        if metric_name not in ['radialFracAtD', 'radialMeanFrac', 'radialCV']:
            raise ValueError(f"metric_name must be one of 'radialFracAtD', 'radialMeanFrac', or 'radialCV'")

        radial_df = self.results['radialdistribution']
        metadata = self.results['metadata']
        
        if object_id is not None:
            # Single object plotting
            if object_id not in metadata:
                raise ValueError(f"Object {object_id} not found in results.")
            
            radial_bins_image = metadata[object_id]['radial_bins_image']
            obj_row = radial_df[radial_df['object_id'] == object_id]
            
            bin_values = []
            for i in range(1, self.num_bins + 1):
                column_name = f"{metric_name}_{i}_{channel_name}"
                if column_name in obj_row.columns:
                    bin_values.append(obj_row[column_name].values[0])
                else:
                    raise ValueError(f"Column '{column_name}' not found in results.")
            
            bin_values = np.array(bin_values)
            norm = plt.Normalize(vmin=bin_values.min(), vmax=bin_values.max())
            
            bin_value_overlay = np.full(image.shape, np.nan, dtype=float)
            for i, value in enumerate(bin_values, 1):
                bin_value_overlay[radial_bins_image == i] = value

            fig, ax = plt.subplots(1, 2, figsize=(16, 8))
            
            ax[0].imshow(image, cmap='gray')
            ax[0].set_title(f"Original Image with Object {object_id}")
            ax[0].axis('off')
            
            ax[1].imshow(image, cmap='gray')
            im = ax[1].imshow(bin_value_overlay, cmap=cmap_name, norm=norm, alpha=alpha)
            ax[1].set_title(f"{metric_name} for {channel_name} (Object {object_id})")
            ax[1].axis('off')

            cbar = fig.colorbar(im, ax=ax[1], orientation='vertical', fraction=0.046, pad=0.04)
            cbar.set_label(f'{metric_name}')

            plt.tight_layout()
            plt.show()
            
        else:
            # All objects plotting
            if radial_df.empty:
                print("No objects to plot.")
                return

            bin_value_overlay = np.full(image.shape, np.nan, dtype=float)
            all_radial_values = []
            
            for obj_id in metadata:
                obj_row = radial_df[radial_df['object_id'] == obj_id]
                radial_bins_image = metadata[obj_id]['radial_bins_image']
                
                bin_values = []
                for i in range(1, self.num_bins + 1):
                    column_name = f"{metric_name}_{i}_{channel_name}"
                    if column_name in obj_row.columns:
                        value = obj_row[column_name].values[0]
                        bin_values.append(value)
                        all_radial_values.append(value)
                    else:
                        raise ValueError(f"Column '{column_name}' not found for object {obj_id}.")
                
                for i, value in enumerate(bin_values, 1):
                    bin_value_overlay[radial_bins_image == i] = value
            
            if not all_radial_values:
                print("No valid radial distribution data to plot.")
                return

            norm = plt.Normalize(vmin=min(all_radial_values), vmax=max(all_radial_values))
            
            fig, ax = plt.subplots(1, 2, figsize=(16, 8))
            ax[0].imshow(image, cmap='gray')
            ax[0].set_title("Original Image (All Objects)")
            ax[0].axis('off')
            
            ax[1].imshow(image, cmap='gray')
            im = ax[1].imshow(bin_value_overlay, cmap=cmap_name, norm=norm, alpha=alpha)
            ax[1].set_title(f"{metric_name} for {channel_name} (All Objects)")
            ax[1].axis('off')

            cbar = fig.colorbar(im, ax=ax[1], orientation='vertical', fraction=0.046, pad=0.04)
            cbar.set_label(f'{metric_name}')

            plt.tight_layout()
            plt.show()


class GranularityAnalyzer(BaseImageAnalyzer):
    """
    Analyzer for object-based granularity in multi-channel images.
    """
    
    def __init__(self, channel_names: Optional[List[str]] = None, 
                 granularity_spectrum_range: List[int] = [1,3,5,8,12,16],
                 object_size: int = 30, subsampling_factor: float = 0.25):
        """
        Initialize GranularityAnalyzer.
        
        Parameters:
        -----------
        channel_names : list, optional
            Default channel names to use
        granularity_spectrum_range : list
            List of granularity sizes to measure
        object_size : int
            Radius for background subtraction
        subsampling_factor : float
            Proportional factor for downsampling
        """
        super().__init__(channel_names)
        self.granularity_spectrum_range = sorted(list(set(granularity_spectrum_range)))
        self.object_size = object_size
        self.subsampling_factor = subsampling_factor
        
        if not all(isinstance(x, int) and x > 0 for x in granularity_spectrum_range):
            raise ValueError("granularity_spectrum_range must be a list of positive integers.")
        
        if not 0.0 < subsampling_factor <= 1.0:
            raise ValueError("subsampling_factor must be between 0 and 1.")
    
    def measure(self, intensity_array: np.ndarray, mask_array: np.ndarray, 
                channel_names: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
        """
        Measure granularity spectrum for all objects across all channels.
        
        Parameters:
        -----------
        intensity_array : ndarray
            A 2D (H, W) or 3D (C, H, W) numpy array
        mask_array : ndarray
            Shape (H, W) with integer labels
        channel_names : list, optional
            Names for each channel
            
        Returns:
        --------
        Dict[str, pd.DataFrame]
            Dictionary with 'granularity' key
        """
        # Validate inputs
        intensity_array, mask_array = self._validate_inputs(intensity_array, mask_array)
        channel_names = self._setup_channel_names(intensity_array, channel_names)
        if channel_names is None:
            return None
        
        if mask_array is None:
            raise ValueError("GranularityAnalyzer requires a mask_array")
        
        # Store data for plotting
        self._store_data(intensity_array, mask_array, channel_names)
        
        # Downsample if requested
        if self.subsampling_factor < 1.0:
            new_shape = (
                int(mask_array.shape[0] * self.subsampling_factor),
                int(mask_array.shape[1] * self.subsampling_factor),
            )
            processed_mask = self._downsample_mask(mask_array, new_shape)
            processed_images = np.array([
                self._downsample_image(intensity_array[c], new_shape)
                for c in range(intensity_array.shape[0])
            ])
        else:
            processed_mask = mask_array
            processed_images = intensity_array

        channel_dfs = []
        for i, channel_name in enumerate(channel_names):
            df = self._process_single_channel(processed_images[i], processed_mask)
            
            # Rename columns to be channel-specific
            rename_dict = {
                col: f'{col}_{channel_name}'
                for col in df.columns if 'granularity' in col
            }
            df = df.rename(columns=rename_dict)
            channel_dfs.append(df)

        if not channel_dfs:
            granularity_df = pd.DataFrame()
        else:
            # Merge all dataframes on 'object_id'
            final_df = channel_dfs[0]
            for i in range(1, len(channel_dfs)):
                if 'object_id' in channel_dfs[i].columns and not channel_dfs[i].empty:
                    final_df = pd.merge(final_df, channel_dfs[i], on='object_id', how='outer')
            granularity_df = self._reorder_columns(final_df, channel_names)
        
        results = {'granularity': granularity_df}
        self.results = results
        return self.results
    
    def _downsample_image(self, image: np.ndarray, new_shape: Tuple[int, int]) -> np.ndarray:
        """Downsample image using interpolation."""
        return transform.resize(image, new_shape, anti_aliasing=True, preserve_range=True)

    def _downsample_mask(self, mask: np.ndarray, new_shape: Tuple[int, int]) -> np.ndarray:
        """Downsample mask using nearest neighbor."""
        return transform.resize(mask, new_shape, order=0, anti_aliasing=False, preserve_range=True).astype(mask.dtype)

    def _calculate_granularity_spectrum(self, obj_image: np.ndarray, obj_mask: np.ndarray, sizes: List[int]) -> List[float]:
        """Calculate granularity spectrum for a single object."""
        granularity_values = []
        initial_signal = np.sum(obj_image[obj_mask])

        if initial_signal == 0:
            return [0.0] * len(sizes)

        opened_image = obj_image.copy()
        for size in sizes:
            selem = disk(size)
            opened_image = morphology.opening(opened_image, selem)
            opened_image *= obj_mask
            
            remaining_signal = np.sum(opened_image)
            granularity = (initial_signal - remaining_signal) / initial_signal
            granularity_values.append(granularity)
            
            if remaining_signal == 0:
                granularity_values.extend([1.0] * (len(sizes) - len(granularity_values)))
                break

        return granularity_values

    def _process_single_channel(self, image_channel: np.ndarray, mask: np.ndarray) -> pd.DataFrame:
        """Process a single image channel to measure granularity."""
        # Background subtraction
        if self.object_size > 0:
            selem = disk(self.object_size)
            background = morphology.opening(image_channel, selem)
            image_channel = np.maximum(image_channel - background, 0)

        # Get unique object labels
        object_labels = np.unique(mask)
        object_labels = object_labels[object_labels > 0]

        results = []
        
        # Calculate granularity for each object
        for obj_label in object_labels:
            obj_mask = (mask == obj_label)
            obj_image = image_channel * obj_mask
            
            granularity_values = self._calculate_granularity_spectrum(
                obj_image, obj_mask, self.granularity_spectrum_range
            )
            
            result_dict = {'object_id': int(obj_label)}
            for i, size in enumerate(self.granularity_spectrum_range):
                result_dict[f'granularity_{size}'] = granularity_values[i]
            
            results.append(result_dict)
            
        if not results:
            return pd.DataFrame({'object_id': pd.Series(dtype=int)})

        return pd.DataFrame(results)

    def plot_granularity(self, channel_to_display: int,
                        granularity_data: Optional[pd.DataFrame] = None,
                        granularity_size_to_display: Optional[int] = None,
                        colormap: str = 'viridis', show: bool = True):
        """
        Plot granularity results with optional overlay.
        
        Parameters:
        -----------
        channel_to_display : int
            The index of the image channel to display
        granularity_data : Optional[pd.DataFrame]
            The DataFrame from measure() for overlay
        granularity_size_to_display : Optional[int]
            The granularity size to visualize
        colormap : str
            The colormap for the granularity overlay
        show : bool
            If True, displays the plot immediately
        """
        if self.intensity_data is None:
            raise ValueError("No data to plot. Run measure() first.")
            
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        
        if not 0 <= channel_to_display < len(self.channel_names):
            raise IndexError(f"channel_to_display {channel_to_display} out of bounds")

        base_image = self.intensity_data[channel_to_display]

        # Plot with Granularity Overlay
        if granularity_data is not None and granularity_size_to_display is not None:
            channel_name = self.channel_names[channel_to_display]
            column_name = f'granularity_{granularity_size_to_display}_{channel_name}'
            if column_name not in granularity_data.columns:
                raise ValueError(f"Column '{column_name}' not found in DataFrame.")
            
            granularity_map = np.zeros_like(self.mask_data, dtype=float)
            for _, row in granularity_data.iterrows():
                if pd.notna(row['object_id']) and pd.notna(row[column_name]):
                    obj_id = int(row['object_id'])
                    granularity_value = row[column_name]
                    granularity_map[self.mask_data == obj_id] = granularity_value

            ax.imshow(base_image, cmap='gray')
            masked_map = np.ma.masked_where(self.mask_data == 0, granularity_map)
            im = ax.imshow(masked_map, cmap=colormap, alpha=0.6)
            
            fig.colorbar(im, ax=ax, label=f"Granularity Value (Size {granularity_size_to_display})")
            ax.set_title(f"Granularity Overlay on Channel {channel_to_display}")
        else:
            # Plot only the intensity image
            ax.imshow(base_image, cmap='gray')
            ax.set_title(f"Intensity Image: Channel {channel_to_display}")

        ax.axis('off')

        if show:
            plt.tight_layout()
            plt.show()


# Helper function for skew calculation
def _skew_intensity(regionmask, intensity_image):
    """Calculate skewness of intensity within a region."""
    return skew(intensity_image[regionmask])



#%% Example usage
if __name__ == "__main__":
    # Example with synthetic data
    np.random.seed(42)
    
    # Create synthetic data
    intensity_data = np.random.randint(0, 255, (2, 120, 120)).astype(np.float64)  # 2 channels
    mask_data = np.zeros((120, 120), dtype=int)
    mask_data[20:40, 20:40] = 1  # Object 1
    mask_data[60:100, 60:100] = 2  # Object 2
    
    channel_names = ['DAPI', 'GFP']
    
    print("Testing IntensityAnalyzer...")
    # Test IntensityAnalyzer
    intensity_analyzer = IntensityAnalyzer(channel_names=channel_names)
    intensity_results = intensity_analyzer.measure(intensity_data, mask_data)
    print("Global results shape:", intensity_results['global'].shape)
    print("Global columns:", list(intensity_results['global'].columns))
    print("Object results shape:", intensity_results['objects'].shape)
    print("Object columns (last 10):", list(intensity_results['objects'].columns)[-10:])
    
    print("\nTesting RadialDistributionAnalyzer...")
    # Test RadialDistributionAnalyzer
    radial_analyzer = RadialDistributionAnalyzer(channel_names=channel_names, num_bins=3)
    radial_results = radial_analyzer.measure(intensity_data, mask_data)
    print("Radial distribution results shape:", radial_results['radialdistribution'].shape)
    print("Radial columns (last 10):", list(radial_results['radialdistribution'].columns)[-10:])
    print("Number of objects in metadata:", len(radial_results['metadata']))
    
    print("\nTesting GranularityAnalyzer...")
    # Test GranularityAnalyzer
    granularity_analyzer = GranularityAnalyzer(
        channel_names=channel_names,
        granularity_spectrum_range=[1, 3, 5],
        subsampling_factor=1.0
    )
    granularity_results = granularity_analyzer.measure(intensity_data, mask_data)
    print("Granularity results shape:", granularity_results['granularity'].shape)
    print("Granularity columns:", list(granularity_results['granularity'].columns))
    
    print("\nTesting plotting functions...")
    # Test plotting (uncomment to see plots)
    # intensity_analyzer.plot_overlay()
    # radial_analyzer.plot_radial_distribution(intensity_data[0], mask_data, 'DAPI')
    # granularity_analyzer.plot_granularity(0, granularity_results['granularity'], 3)
    
    print("All tests completed successfully!")
