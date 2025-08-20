
#%%
import pandas as pd
import re, os
from glob import glob
from os.path import join, dirname, basename, exists
from natsort import natsorted


class ImageDataSet:
    def __init__(
        self,
        data_dir: str = None,
        subset_pattern: str = None,
        image_subdir: str = "images",
        image_suffix: str = '.tiff',
        mask_suffix: str = '.png',
        remove_na_row: bool = True,
        position_metadata: str = None,
        pos_index_col: str = "P Index",
        image_extractor: str = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_ch?(?P<channel>[0-9]{1,})(?P<self_generated>.*)',
        mask_extractor: str = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_ch?(?P<channel>[0-9]{1,})_(?P<mask_name>.*)'
    ):
        """
        Initializes the ImageDataSet class.

        Args:
            data_dir: str, path to the working directory.
            subset_pattern: str, pattern to subset the images.
            image_subdir: str, path to the image directory.
            image_suffix: str, suffix of the image file.
            mask_suffix: str, suffix of the mask file (None for no use).
            remove_na_row: bool, remove the rows with any NA values.
            position_metadata: str, path to the position metadata file (None for no usage).
            pos_index_col: str, which column in position_metadata to use to merge with the dataframe.
            image_extractor: str, regex pattern to extract metadata from image filenames.
            mask_extractor: str, regex pattern to extract metadata from mask filenames.
        """
        self.data_dir = data_dir
        self.subset_pattern = subset_pattern
        self.image_subdir = image_subdir
        self.image_suffix = image_suffix
        self.mask_suffix = mask_suffix
        self.remove_na_row = remove_na_row
        self.position_metadata = position_metadata
        self.pos_index_col = pos_index_col
        self.image_extractor = image_extractor
        self.mask_extractor = mask_extractor

        self.df = None
        self.metadata_colnames = []
        self.intensity_colnames = []
        self.mask_colnames = []
        self.df = pd.DataFrame()

        self._process_images_to_dataset()
        print(self.__repr__())

    def _process_images_to_dataset(self):
        """
        Internal method to process images and build the dataset dataframe.
        """
        self._validate_extractors()

        if not exists(self.data_dir):
            print("data_dir does not exist.")
            return

        image_df, image_index_cols = self._load_and_process_images()
        if image_df is None:
            return

        mask_df = self._load_and_process_masks(image_index_cols)

        self._merge_dataframes(image_df, mask_df)
        
        if self.remove_na_row:
            self.df = self.df.dropna()

        self.df = self.df.reset_index()
        # print(f"{len(self.df)} groups after merging images and masks")
        
        self._add_position_metadata()
        self._cast_dtypes()

    def _validate_extractors(self):
        """Validates the image and mask regex extractors."""
        image_extract_cols = re.findall(r'\?P<([^>]+)>', self.image_extractor)
        if not ('channel' in image_extract_cols or 'self_generated' in image_extract_cols):
            raise ValueError("image_extractor must contain 'channel' or 'self_generated'")
        
        if self.mask_suffix is not None:
            mask_extract_cols = re.findall(r'\?P<([^>]+)>', self.mask_extractor)
            if 'mask_name' not in mask_extract_cols:
                raise ValueError("mask_extractor must contain 'mask_name'")

    def _load_and_process_images(self):
        """Loads and processes image files into a DataFrame."""
        image_paths = glob(join(self.data_dir, self.image_subdir, '**/*') + self.image_suffix, recursive=True)
        
        if self.subset_pattern is not None:
            image_paths = [f for f in image_paths if re.search(self.subset_pattern, f)]
        
        if len(image_paths) == 0:
            print("no valid images found, check directory or subset_pattern")
            return None, None

        image_df = pd.DataFrame({
            "directory": [dirname(x) for x in image_paths],
            "filename": [basename(x) for x in image_paths]
        })

        image_extract_cols = re.findall(r'\?P<([^>]+)>', self.image_extractor)
        image_meta = image_df["filename"].str.extract(self.image_extractor + f'{self.image_suffix}')
        image_meta["channel"] = image_meta["channel"].astype("str") + image_meta["self_generated"]
        image_meta.pop('self_generated')
        
        image_df = pd.concat([image_df, image_meta], axis=1)
        
        image_index_cols = ['directory'] + [x for x in image_extract_cols if x not in ['channel', 'self_generated']]
        image_df = image_df.set_index(image_index_cols)
        
        image_df["channel"] = "ch" + image_df["channel"]
        self.intensity_colnames = natsorted(image_df['channel'].unique().tolist())
        
        image_df = image_df.pivot(columns="channel", values="filename")
        # print(f"{len(image_df)} grouped intensity images: {list(image_df)}")
        
        self.metadata_colnames = image_index_cols
        return image_df, image_index_cols

    def _load_and_process_masks(self, image_index_cols):
        """Loads and processes mask files into a DataFrame."""
        if self.mask_suffix is None:
            self.mask_colnames = None
            return None
        
        mask_paths = glob(join(self.data_dir, self.image_subdir, '**/*') + self.mask_suffix, recursive=True)
        if len(mask_paths) == 0:
            self.mask_colnames = None
            return None
        else:
            mask_df = pd.DataFrame({
                "directory": [dirname(x) for x in mask_paths],
                "filename": [basename(x) for x in mask_paths]
            })
            
            mask_meta = mask_df["filename"].str.extract(self.mask_extractor + f'{self.mask_suffix}')
            mask_df = pd.concat([mask_df, mask_meta], axis=1)
            
            mask_df = mask_df.set_index(image_index_cols)
            self.mask_colnames = mask_df['mask_name'].unique().tolist()
            
            mask_df = mask_df.pivot(columns="mask_name", values="filename")
            return mask_df

    def _merge_dataframes(self, image_df, mask_df):
        """Merges image and mask dataframes."""
        if mask_df is None:
            self.df = image_df
        else:
            self.df = pd.merge(image_df, mask_df, how='left', left_index=True, right_index=True)

    def _add_position_metadata(self):
        """Adds position metadata if a file is provided."""
        if self.position_metadata is not None and exists(join(self.data_dir, self.position_metadata)):
            pos_meta = pd.read_csv(join(self.data_dir, self.position_metadata))
            pos_meta = pos_meta.dropna(subset=[self.pos_index_col])
            pos_meta = pos_meta.drop_duplicates(subset=['prefix', 'Time [s]', self.pos_index_col], keep='first')
            
            if not pd.isna(pos_meta.iloc[0]["Position Name"]):
                print("add position well matadata")
                pos_meta_well = pos_meta["Position Name"].str.extract("(?P<row>[A-Z])(?P<column>[0-9]+)#(?P<field>.+)")
                pos_meta_well["well"] = pos_meta_well["row"] + pos_meta_well["column"]
                
                if "T Index" in pos_meta.columns:
                    pos_meta["T Index"] = pos_meta["T Index"].fillna(1)
                    pos_meta = pd.concat([pos_meta["prefix"],
                                          pos_meta[self.pos_index_col].rename("position"),
                                          pos_meta["T Index"].rename("timepoint"),
                                          pos_meta_well], axis=1)
                    pos_meta = pos_meta.drop(["row", "column"], axis=1)
                    pos_meta["prefix"] = pos_meta["prefix"].astype("str")
                    self.df = pd.merge(pos_meta, self.df, how="right", on=["prefix", "position", "timepoint"])
                else:
                    pos_meta = pd.concat([pos_meta["prefix"],
                                          pos_meta[self.pos_index_col].rename("position"),
                                          pos_meta_well], axis=1)
                    pos_meta = pos_meta.drop(["row", "column"], axis=1)
                    pos_meta["prefix"] = pos_meta["prefix"].astype("str")
                    self.df = pd.merge(pos_meta, self.df, how="right", on=["prefix", "position"])

                self.metadata_colnames.extend(["well", "field"])

    def _cast_dtypes(self):
        """Casts column data types for consistency."""
        if 'prefix' in self.df.columns:
            self.df['prefix'] = self.df['prefix'].astype("str")
        if 'position' in self.df.columns:
            self.df['position'] = self.df['position'].astype("int64")
        if 'timepoint' in self.df.columns:
            self.df['timepoint'] = self.df['timepoint'].astype("int64")
            
    def export_to_cellprofiler_dataloader(self, filename: str = 'cp_dataloader.csv'):
        """
        Converts the internal DataFrame to a CellProfiler-style format.

        Returns:
            pd.DataFrame: The DataFrame formatted for CellProfiler.
        """
        if self.df is None:
            print("DataFrame is not initialized.")
            return None

        cp_df = self.df.copy()

        # Image paths and filenames
        for ch in self.intensity_colnames:
            cp_df[f'Image_PathName_{ch}'] = cp_df['directory']
            cp_df = cp_df.rename(columns={f'{ch}': f'Image_FileName_{ch}'})
        
        # Mask paths and filenames
        if self.mask_colnames is not None:
            for mask_colname in self.mask_colnames:
                cp_df[f'Image_ObjectsPathName_mask_{mask_colname}'] = cp_df['directory']
                cp_df = cp_df.rename(columns={f'{mask_colname}': f'Image_ObjectsFileName_mask_{mask_colname}'})
        
        # Metadata
        for meta in self.metadata_colnames:
            cp_df = cp_df.rename(columns={f'{meta}': f'Metadata_{meta}'})
        
        # save to file
        if filename:
            cp_df.to_csv(os.path.join(self.data_dir, filename), index=False)

        return cp_df


    def __repr__(self):
        info = (f"ImageDataSet\ndirectory: [{self.data_dir}]\nintensity: {self.intensity_colnames} \nmask: {self.mask_colnames} \ngroups: [{len(self.df)}] \nmetadata: {self.metadata_colnames}")
        return info


#%%
if __name__ == '__main__':

    dataset = ImageDataSet('/media/hao/Data/Project_MYC/2024-10-31_MYCN_saturation_curve_adding_truncation')
    # dataset.export_to_cellprofiler_dataloader()
    