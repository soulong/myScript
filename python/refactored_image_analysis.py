import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
from PIL import Image
import tifffile
import pandas as pd
from collections import defaultdict


class ImageDataset:

    def __init__(self, 
                directory: str,
                intensity_regex: str = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_ch?(?P<channel>[0-9]{1,})(?P<self_generated>.*)',
                mask_regex: str = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_ch?(?P<channel>[0-9]{1,})_(?P<mask_name>.*)'


                intensity_pattern: str = "{prefix}__t{timepoint}_p{position}_z{zstack}_ch{channel}.tiff",
                mask_pattern: str = "{prefix}__t{timepoint}_p{position}_z{zstack}_ch{channel}_cp_masks_{mask}.png",
                index_cols: List[str] = ['prefix','timepoint','position','stack','self_generated'], 
                external_metadata_csv: Optional[str] = None,
                column_to_use = ['prefix', 'timepoint', 'position', 'well']
                ):
        self.directory = Path(directory)

        self.intensity_files = self.get_files(self.directory, Path(intensity_pattern).suffix)
        self.intensity_regex, self.intensity_fields = self.pattern_to_regex(intensity_pattern)
        self.intensity = self.pharse_filename(self.intensity_files, self.intensity_regex)
        self.intensity_df = pd.DataFrame(self.intensity)

        self.mask_files = self.get_files(self.directory,  Path(mask_pattern).suffix)
        self.mask_regex, self.mask_fields = self.pattern_to_regex(mask_pattern)
        self.mask = self.pharse_filename(self.mask_files, self.mask_regex)
        self.mask_df = pd.DataFrame(self.mask)

        self.external_metadata_csv = external_metadata_csv
        self.column_to_use = column_to_use

    # def __getitem__(self):
        

    # def __len__(self):
    #     pass

    # def __repr__(self):
    #     # print(f'diretory {str(self.intensity)}')
    #     # print(f'total files {len([self.mask])}')
    #     pass
    
    @staticmethod
    def get_files(path, suffix: str = '.tiff') -> List[str]:
        files = []
        for filepath in path.iterdir():
            if filepath.is_file() and filepath.suffix == suffix:
                files.append(filepath)
        return files

    @staticmethod
    def pattern_to_regex(pattern: str) -> Tuple[str, List[str]]:
        """
        Convert a pattern template to regex and extract field names.
        
        Args:
            pattern (str): Pattern with {field} placeholders
            
        Returns:
            Tuple[str, List[str]]: (regex_pattern, field_names)
        """
        # Find all field placeholders
        field_matches = re.findall(r'\{(\w+)\}', pattern)
        # Escape special regex characters in the pattern
        escaped_pattern = re.escape(pattern)
        # Replace escaped placeholders with capture groups
        regex_pattern = escaped_pattern
        for field in field_matches:
            escaped_placeholder = re.escape('{' + field + '}')
            # Use named capture groups for better organization
            regex_pattern = regex_pattern.replace(escaped_placeholder, f'(?P<{field}>[^_./\\\\]+)')
        # Make it match the entire filename
        regex_pattern = '^' + regex_pattern + '$'
        return regex_pattern, field_matches
    
    @staticmethod
    def pharse_filename(paths: List[Path], regex: str) -> List[Dict]:
        metadata = []
        for f_path in paths:
            filename = f_path.name
            filepath = f_path.parent
            filename_match = re.match(regex, filename)
            if filename_match:
                meta = filename_match.groupdict()
            else:
                meta = {}
            meta['filepath'] = str(filepath)
            meta['filename'] = str(filename)
            metadata.append(meta)
        return metadata

    @staticmethod
    def pivot_df_wide(df: pd.DataFrame, index: List[str]) -> pd.DataFrame:
        # df = images_df.reset_index()
        index_cols = ['directory'] + [x for x in image_extract_cols if x not in ['channel', 'self_generated']]
        image_df = image_df.set_index(image_index_cols)
        # add prefix to channel
        image_df["channel"] = "ch" + image_df["channel"]
        # get unique channel names
        image_colnames = natsorted(image_df['channel'].unique().tolist())
        # reshape channel to wide
        image_df = image_df.pivot(columns="channel", values="filename")
        print(f"{len(image_df)} grouped intensity images: {list(image_df)}")

        return df_new
    
    def merge_group(x1, x2):
        pass

x = ImageDataset('/media/hao/Data/test')
x.intensity_fields
x.mask_fields
x.intensity_df.iloc[13]
x.intensity_df.columns


pd.DataFrame(x.intensity)


x._pharse_filename(x.intensity_files[0:5], x.intensity_regex)





class MicroscopyDataset:
    """
    A flexible class to handle microscopy image datasets with customizable naming conventions.
    
    Supports M:N relationships between masks and intensity channels, where:
    - M different mask types can exist for each image position
    - N different intensity channels can exist for each image position
    - External metadata can be loaded from CSV files
    
    The combined data structure organizes as: ch1 ch2 chX mask1 mask2 maskY
    """
    
    def __init__(self, 
                 directory: str, 
                 intensity_pattern: str = "{prefix}__t{timepoint}_p{position}_z{zstack}_ch{channel}.tiff",
                 mask_pattern: str = "{prefix}__t{timepoint}_p{position}_z{zstack}_ch{channel}_cp_masks_{mask}.png",
                 external_metadata_csv: Optional[str] = None,
                 column_to_use = ['prefix', 'timepoint', 'position', 'well']
                 ):
        """
        Initialize the microscopy dataset with flexible naming patterns.
        
        Args:
            directory (str): Path to directory containing image files
            intensity_pattern (str): Template pattern for intensity files using {field} placeholders
            mask_pattern (str, optional): Template pattern for mask files. If None, will auto-detect masks
            external_metadata_csv (str, optional): Path to CSV file with additional metadata
            column_to_use: columns to be used in external_metadata_csv
        
        Examples of patterns:
            "{prefix}_t{timepoint}_z{stack}_ch{channel}_u1{userdefine}.tif"
            "{experiment}_{condition}_T{time}_Z{z}_C{channel}.tiff"
            "{sample}_d{day}_t{timepoint}_p{position}_ch{channel}_masks_{mask_type}.png"
        """
        self.directory = Path(directory)
        self.intensity_pattern = intensity_pattern
        self.mask_pattern = mask_pattern
        self.external_metadata_csv = external_metadata_csv
        self.column_to_use = column_to_use
                
        # Parse patterns to create regex and extract field names
        self.intensity_regex, self.intensity_fields = self._pattern_to_regex(intensity_pattern)
        
        if mask_pattern:
            self.mask_regex, self.mask_fields = self._pattern_to_regex(mask_pattern)
        else:
            # Auto-detect mask pattern by looking for common mask indicators
            self.mask_regex = None
            self.mask_fields = []
        
        # Storage for parsed file information
        self.files: Dict[str, Dict] = {}  # key: filepath, value: metadata
        
        # Organized indices for all unique field values
        self.unique_values: Dict[str, set] = {}
        self.files_by_metadata: Dict[str, List[str]] = {}
        
        # M:N relationship storage
        self.image_groups: Dict[str, Dict] = {}  # key: group_id, value: {intensities: [], masks: [], metadata: {}}
        
        # External metadata
        self.external_metadata: Dict[str, Dict] = {}
        
        self._scan_directory()
        self._organize_files()
        self._load_external_metadata()
        self._build_image_groups()
    
    def _pattern_to_regex(self, pattern: str) -> Tuple[str, List[str]]:
        """
        Convert a pattern template to regex and extract field names.
        
        Args:
            pattern (str): Pattern with {field} placeholders
            
        Returns:
            Tuple[str, List[str]]: (regex_pattern, field_names)
        """
        # Find all field placeholders
        field_matches = re.findall(r'\{(\w+)\}', pattern)
        
        # Escape special regex characters in the pattern
        escaped_pattern = re.escape(pattern)
        
        # Replace escaped placeholders with capture groups
        regex_pattern = escaped_pattern
        for field in field_matches:
            escaped_placeholder = re.escape('{' + field + '}')
            # Use named capture groups for better organization
            regex_pattern = regex_pattern.replace(escaped_placeholder, f'(?P<{field}>[^_./\\\\]+)')
        
        # Make it match the entire filename
        regex_pattern = '^' + regex_pattern + '$'
        
        return regex_pattern, field_matches
    
    def _scan_directory(self):
        """Scan directory and parse file names."""
        if not self.directory.exists():
            raise FileNotFoundError(f"Directory {self.directory} does not exist")
        
        for filepath in self.directory.iterdir():
            if filepath.is_file():
                self._parse_filename(filepath)
    
    def _parse_filename(self, filepath: Path):
        """Parse individual filename and extract metadata."""
        filename = filepath.name
        
        # Try intensity pattern first
        intensity_match = re.match(self.intensity_regex, filename)
        if intensity_match:
            metadata = intensity_match.groupdict()
            
            # Convert numeric fields to integers where possible
            metadata = self._convert_numeric_fields(metadata)
            
            # Add file type and path information
            metadata.update({
                'type': 'intensity',
                'filepath': str(filepath),
                'filename': filename,
                'extension': filepath.suffix.lower()
            })
            
            self.files[str(filepath)] = metadata
            self._update_unique_values(metadata)
            return
        
        # Try mask pattern if provided
        if self.mask_regex:
            mask_match = re.match(self.mask_regex, filename)
            if mask_match:
                metadata = mask_match.groupdict()
                metadata = self._convert_numeric_fields(metadata)
                
                metadata.update({
                    'type': 'mask',
                    'filepath': str(filepath),
                    'filename': filename,
                    'extension': filepath.suffix.lower()
                })
                
                self.files[str(filepath)] = metadata
                self._update_unique_values(metadata)
                return

    #     # Auto-detect masks by common patterns
    #     if not self.mask_regex:
    #         self._try_auto_detect_mask(filepath, filename)
    
    # def _try_auto_detect_mask(self, filepath: Path, filename: str):
    #     """Try to auto-detect mask files based on common patterns."""
    #     mask_indicators = ['cp_masks', 'mask', 'seg', 'segmentation']
        
    #     for indicator in mask_indicators:
    #         if indicator in filename.lower():
    #             # Create a simple metadata structure for auto-detected masks
    #             metadata = {
    #                 'type': 'mask',
    #                 'filepath': str(filepath),
    #                 'filename': filename,
    #                 'extension': filepath.suffix.lower(),
    #                 'mask_type': indicator,
    #                 'auto_detected': True
    #             }
                
    #             # Try to extract some common fields from filename using heuristics
    #             self._extract_common_fields(metadata, filename)
                
    #             self.files[str(filepath)] = metadata
    #             self._update_unique_values(metadata)
    #             return
    
    # def _extract_common_fields(self, metadata: Dict, filename: str):
    #     """Extract common fields from filename using heuristics."""
    #     # Common patterns to look for
    #     patterns = [
    #         (r't(\d+)', 'timepoint'),
    #         (r'p(\d+)', 'position'),
    #         (r'z(\d+)', 'zstack'),
    #         (r'ch(\d+)', 'channel'),
    #         (r'c(\d+)', 'channel'),
    #         (r'T(\d+)', 'timepoint'),
    #         (r'Z(\d+)', 'zstack'),
    #         (r'C(\d+)', 'channel'),
    #     ]
        
    #     for pattern, field_name in patterns:
    #         match = re.search(pattern, filename, re.IGNORECASE)
    #         if match:
    #             metadata[field_name] = int(match.group(1))
    
    def _convert_numeric_fields(self, metadata: Dict) -> Dict:
        """Convert string values to integers where appropriate."""
        converted = {}
        for key, value in metadata.items():
            if isinstance(value, str) and value.isdigit():
                converted[key] = int(value)
            else:
                converted[key] = value
        return converted
    
    def _update_unique_values(self, metadata: Dict):
        """Update sets of unique values for each metadata field."""
        for key, value in metadata.items():
            if key not in ['filepath', 'filename', 'type']:  # Skip non-metadata fields
                if key not in self.unique_values:
                    self.unique_values[key] = set()
                self.unique_values[key].add(value)
    
    def _organize_files(self):
        """Organize files by metadata combinations for quick lookup."""
        for filepath, metadata in self.files.items():
            # Create lookup keys for each metadata field
            for key, value in metadata.items():
                if key not in ['filepath', 'filename']:  # Skip path information
                    lookup_key = f"{key}_{value}"
                    if lookup_key not in self.files_by_metadata:
                        self.files_by_metadata[lookup_key] = []
                    self.files_by_metadata[lookup_key].append(filepath)
    
    def _load_external_metadata(self):
        """Load external metadata from CSV file if provided."""
        if not self.external_metadata_csv:
            return
        
        csv_path = Path(self.external_metadata_csv)
        if not csv_path.exists():
            print(f"Warning: External metadata CSV file {csv_path} not found")
            return
        
        try:
            df = pd.read_csv(csv_path)
            
            # Convert DataFrame to dictionary keyed by combination of identifying fields
            for _, row in df.iterrows():
                # Create a key from common identifying fields
                key_fields = []
                for field in self.column_to_use:
                    if field in row and pd.notna(row[field]):
                        key_fields.append(f"{field}_{row[field]}")
                
                if key_fields:
                    key = "_".join(key_fields)
                    self.external_metadata[key] = row.to_dict()
            
            print(f"Loaded external metadata for {len(self.external_metadata)} entries")
        except Exception as e:
            print(f"Warning: Could not load external metadata: {e}")
    
    def _build_image_groups(self):
        """Build M:N relationships between masks and intensity channels."""
        # Group files by their identifying metadata (excluding channel and mask-specific fields)
        groups = defaultdict(lambda: {'intensities': [], 'masks': [], 'metadata': {}})
        
        for filepath, metadata in self.files.items():
            # Create group key from non-channel/mask specific fields
            key_fields = []
            for field in self.column_to_use:
                if field in metadata:
                    key_fields.append(f"{field}_{metadata[field]}")
            
            # # Add other custom fields (excluding channel, mask, type, filepath, filename, extension)
            # exclude_fields = {'channel', 'mask', 'type', 'filepath', 'filename', 'extension'}
            # for field, value in metadata.items():
            #     if field not in exclude_fields and field not in ['prefix', 'timepoint', 'position', 'zstack']:
            #         key_fields.append(f"{field}_{value}")
            
            group_key = "_".join(key_fields) if key_fields else "default"
            
            # Add to appropriate list
            if metadata['type'] == 'intensity':
                groups[group_key]['intensities'].append(metadata)
            elif metadata['type'] == 'mask':
                groups[group_key]['masks'].append(metadata)
            
            # Store common metadata (first occurrence)
            if not groups[group_key]['metadata']:
                common_metadata = {k: v for k, v in metadata.items()}
                                #  if k not in exclude_fields}
                groups[group_key]['metadata'] = common_metadata
        
        # Add external metadata if available
        for group_key, group_data in groups.items():
            if group_key in self.external_metadata:
                group_data['metadata'].update(self.external_metadata[group_key])
        
        self.image_groups = dict(groups)
    
    def get_image_group(self, group_key: str) -> Optional[Dict]:
        """
        Get a specific image group by key.
        
        Args:
            group_key (str): Group identifier
            
        Returns:
            Dict: Group data with intensities, masks, and metadata
        """
        return self.image_groups.get(group_key)
    
    def get_image_groups_by_criteria(self, **criteria) -> Dict[str, Dict]:
        """
        Get image groups matching specified criteria.
        
        Args:
            **criteria: Any metadata field as keyword argument
            
        Returns:
            Dict[str, Dict]: Dictionary of matching groups
        """
        matching_groups = {}
        
        for group_key, group_data in self.image_groups.items():
            metadata = group_data['metadata']
            match = True
            
            for criterion_key, criterion_value in criteria.items():
                if criterion_key not in metadata or metadata[criterion_key] != criterion_value:
                    match = False
                    break
            
            if match:
                matching_groups[group_key] = group_data
        
        return matching_groups
    
    def get_combined_image_stack(self, group_key: str, 
                                intensity_channels: Optional[List[Any]] = None,
                                mask_types: Optional[List[Any]] = None) -> Tuple[np.ndarray, Dict]:
        """
        Get combined image stack for a group with format: [ch1, ch2, ..., chN, mask1, mask2, ..., maskM]
        
        Args:
            group_key (str): Group identifier
            intensity_channels (List, optional): List of specific channels to include
            mask_types (List, optional): List of specific mask types to include
            
        Returns:
            Tuple[np.ndarray, Dict]: Combined image stack and metadata
        """
        if group_key not in self.image_groups:
            raise KeyError(f"Group {group_key} not found")
        
        group_data = self.image_groups[group_key]
        
        # Filter and sort intensity channels
        intensities = group_data['intensities']
        if intensity_channels is not None:
            intensities = [img for img in intensities if img.get('channel') in intensity_channels]
        
        # Sort by channel number
        intensities.sort(key=lambda x: x.get('channel', 0))
        
        # Filter and sort masks
        masks = group_data['masks']
        if mask_types is not None:
            masks = [img for img in masks if img.get('mask') in mask_types or img.get('mask_type') in mask_types]
        
        # Sort by mask type/name
        masks.sort(key=lambda x: str(x.get('mask', x.get('mask_type', ''))))
        
        # Load all images
        image_stack = []
        image_info = []
        
        # Add intensity channels
        for intensity_meta in intensities:
            img = self.load_image(intensity_meta['filepath'])
            if img.ndim == 2:
                img = img[np.newaxis, ...]  # Add channel dimension
            elif img.ndim == 3 and img.shape[0] > img.shape[2]:  # Likely HWC format
                img = np.transpose(img, (2, 0, 1))  # Convert to CHW
            
            for i in range(img.shape[0]):
                image_stack.append(img[i])
                info = intensity_meta.copy()
                info['stack_type'] = 'intensity'
                info['stack_index'] = len(image_stack) - 1
                if img.shape[0] > 1:
                    info['sub_channel'] = i
                image_info.append(info)
        
        # Add masks
        for mask_meta in masks:
            img = self.load_image(mask_meta['filepath'])
            if img.ndim == 2:
                img = img[np.newaxis, ...]  # Add channel dimension
            elif img.ndim == 3 and img.shape[0] > img.shape[2]:  # Likely HWC format
                img = np.transpose(img, (2, 0, 1))  # Convert to CHW
            
            for i in range(img.shape[0]):
                image_stack.append(img[i])
                info = mask_meta.copy()
                info['stack_type'] = 'mask'
                info['stack_index'] = len(image_stack) - 1
                if img.shape[0] > 1:
                    info['sub_mask'] = i
                image_info.append(info)
        
        if not image_stack:
            raise ValueError(f"No images found for group {group_key}")
        
        # Stack all images
        combined_stack = np.stack(image_stack, axis=0)
        
        # Create comprehensive metadata
        combined_metadata = {
            'group_key': group_key,
            'group_metadata': group_data['metadata'],
            'image_info': image_info,
            'n_intensity_channels': len([info for info in image_info if info['stack_type'] == 'intensity']),
            'n_mask_channels': len([info for info in image_info if info['stack_type'] == 'mask']),
            'total_channels': len(image_info),
            'shape': combined_stack.shape
        }
        
        return combined_stack, combined_metadata
    
    def get_all_combined_stacks(self, **criteria) -> List[Tuple[np.ndarray, Dict]]:
        """
        Get combined image stacks for all groups matching criteria.
        
        Args:
            **criteria: Any metadata field as keyword argument
            
        Returns:
            List[Tuple[np.ndarray, Dict]]: List of (image_stack, metadata) tuples
        """
        matching_groups = self.get_image_groups_by_criteria(**criteria)
        
        results = []
        for group_key in matching_groups.keys():
            try:
                stack, metadata = self.get_combined_image_stack(group_key)
                results.append((stack, metadata))
            except Exception as e:
                print(f"Warning: Could not create stack for group {group_key}: {e}")
        
        return results
    
    def get_files_by_criteria(self, **criteria) -> List[Dict]:
        """
        Get files matching specified criteria using keyword arguments.
        
        Args:
            **criteria: Any metadata field as keyword argument
            
        Examples:
            get_files_by_criteria(timepoint=1, channel=2)
            get_files_by_criteria(prefix="experiment1", userdefine="condition_A")
            get_files_by_criteria(type="intensity")
        
        Returns:
            List[Dict]: List of file metadata dictionaries
        """
        results = []
        
        for filepath, metadata in self.files.items():
            match = True
            
            for criterion_key, criterion_value in criteria.items():
                if criterion_key not in metadata or metadata[criterion_key] != criterion_value:
                    match = False
                    break
            
            if match:
                results.append(metadata)
        
        return results
    
    def get_image_by_index(self, index: int, file_type: Optional[str] = None) -> Tuple[np.ndarray, Dict]:
        """
        Get image by index.
        
        Args:
            index (int): Index of the file
            file_type (str, optional): Filter by file type ('intensity', 'mask', etc.)
        
        Returns:
            Tuple[np.ndarray, Dict]: Image array and metadata
        """
        if file_type:
            files = [metadata for metadata in self.files.values() if metadata.get('type') == file_type]
        else:
            files = list(self.files.values())
        
        if index >= len(files):
            raise IndexError(f"Index {index} out of range (max: {len(files)-1})")
        
        metadata = files[index]
        image = self.load_image(metadata['filepath'])
        return image, metadata
    
    def get_images_by_criteria(self, **criteria) -> List[Tuple[np.ndarray, Dict]]:
        """
        Get images matching specified criteria.
        
        Args:
            **criteria: Any metadata field as keyword argument
        
        Returns:
            List[Tuple[np.ndarray, Dict]]: List of (image_array, metadata) tuples
        """
        matching_files = self.get_files_by_criteria(**criteria)
        
        results = []
        for metadata in matching_files:
            image = self.load_image(metadata['filepath'])
            results.append((image, metadata))
        
        return results
    
    def load_image(self, filepath: str) -> np.ndarray:
        """
        Load image from filepath.
        
        Args:
            filepath (str): Path to image file
        
        Returns:
            np.ndarray: Image array
        """
        filepath = Path(filepath)
        
        try:
            if filepath.suffix.lower() in ['.tiff', '.tif']:
                return tifffile.imread(str(filepath))
            else:
                return np.array(Image.open(filepath))
        except Exception as e:
            raise ValueError(f"Could not load image {filepath}: {e}")
    
    def get_unique_values(self, field: str) -> List[Any]:
        """
        Get all unique values for a specific metadata field.
        
        Args:
            field (str): Metadata field name
            
        Returns:
            List: Sorted list of unique values
        """
        if field in self.unique_values:
            values = list(self.unique_values[field])
            try:
                return sorted(values)
            except TypeError:
                # If values can't be sorted (mixed types), return as list
                return values
        return []
    
    def get_field_names(self) -> List[str]:
        """Get all available metadata field names."""
        return list(self.unique_values.keys())
    
    def get_group_keys(self) -> List[str]:
        """Get all available image group keys."""
        return list(self.image_groups.keys())
    
    def filter_files(self, **criteria) -> 'MicroscopyDataset':
        """
        Create a new dataset instance with files matching the criteria.
        
        Args:
            **criteria: Any metadata field as keyword argument
            
        Returns:
            MicroscopyDataset: New filtered dataset instance
        """
        # Create a new instance with the same directory and patterns
        filtered_dataset = MicroscopyDataset.__new__(MicroscopyDataset)
        filtered_dataset.directory = self.directory
        filtered_dataset.intensity_pattern = self.intensity_pattern
        filtered_dataset.mask_pattern = self.mask_pattern
        filtered_dataset.external_metadata_csv = self.external_metadata_csv
        filtered_dataset.intensity_regex = self.intensity_regex
        filtered_dataset.intensity_fields = self.intensity_fields
        filtered_dataset.mask_regex = self.mask_regex
        filtered_dataset.mask_fields = self.mask_fields
        
        # Filter files based on criteria
        matching_files = self.get_files_by_criteria(**criteria)
        filtered_dataset.files = {metadata['filepath']: metadata for metadata in matching_files}
        
        # Rebuild indices
        filtered_dataset.unique_values = {}
        filtered_dataset.files_by_metadata = {}
        filtered_dataset.external_metadata = self.external_metadata
        
        for metadata in filtered_dataset.files.values():
            filtered_dataset._update_unique_values(metadata)
        filtered_dataset._organize_files()
        filtered_dataset._build_image_groups()
        
        return filtered_dataset
    
    def summary(self) -> Dict:
        """
        Get summary statistics of the dataset.
        
        Returns:
            Dict: Summary information
        """
        summary_data = {
            'total_files': len(self.files),
            'total_groups': len(self.image_groups),
            'directory': str(self.directory),
            'intensity_pattern': self.intensity_pattern,
            'mask_pattern': self.mask_pattern,
            'external_metadata_csv': self.external_metadata_csv,
            'field_names': self.get_field_names(),
            'file_types': {},
            'groups_summary': {}
        }
        
        # Count files by type
        for metadata in self.files.values():
            file_type = metadata.get('type', 'unknown')
            if file_type not in summary_data['file_types']:
                summary_data['file_types'][file_type] = 0
            summary_data['file_types'][file_type] += 1
        
        # Add unique values for each field
        for field_name in self.get_field_names():
            if field_name not in ['extension', 'auto_detected']:  # Skip technical fields
                summary_data[f'unique_{field_name}'] = self.get_unique_values(field_name)
        
        # Summarize groups
        intensity_counts = []
        mask_counts = []
        for group_data in self.image_groups.values():
            intensity_counts.append(len(group_data['intensities']))
            mask_counts.append(len(group_data['masks']))
        
        if intensity_counts:
            summary_data['groups_summary'] = {
                'avg_intensities_per_group': np.mean(intensity_counts),
                'avg_masks_per_group': np.mean(mask_counts),
                'max_intensities_per_group': max(intensity_counts),
                'max_masks_per_group': max(mask_counts),
                'min_intensities_per_group': min(intensity_counts),
                'min_masks_per_group': min(mask_counts)
            }
        
        return summary_data
    
    def __len__(self) -> int:
        """Return total number of files."""
        return len(self.files)
    
    def __str__(self) -> str:
        """String representation of the dataset."""
        summary = self.summary()
        file_type_str = ", ".join([f"{k}: {v}" for k, v in summary['file_types'].items()])
        
        return (f"MicroscopyDataset(\n"
                f"  Directory: {summary['directory']}\n"
                f"  Total files: {summary['total_files']}\n"
                f"  Total groups: {summary['total_groups']}\n"
                f"  File types: {file_type_str}\n"
                f"  Available fields: {', '.join(summary['field_names'])}\n"
                f"  Intensity pattern: {summary['intensity_pattern']}\n"
                f"  Mask pattern: {summary['mask_pattern']}\n"
                f"  External metadata: {summary['external_metadata_csv']}\n"
                f")")
    

x1 = MicroscopyDataset("/media/hao/Data/Project_r-secretase/2025-05-16_sgRNA_validation/images")
print(x1)
x1.column_to_use
x1.external_metadata
x1.files[1:10]

x1.files_by_metadata[1:10]
x1.filter_files(position=10)
x1.image_groups
x1[1]

x2 = MicroscopyDataset("/media/hao/Data/Project_r-secretase/2025-05-16_sgRNA_validation/images",
                       external_metadata_csv="/media/hao/Data/Project_r-secretase/2025-05-16_sgRNA_validation/position.csv")
print(x2)