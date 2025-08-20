
#%%
import numpy as np
import pandas as pd
import re
from os.path import join, dirname, basename, exists, splitext
from skimage.io import imread
from skimage.transform import rescale
from tqdm import tqdm
import torch
from cellpose import models, io as cp_io
from skimage.morphology import closing


class CellposeSegmentation:
    def __init__(
        self,
        channel_names: list[list[str]] = [['ch1','ch2'], ['ch3']],
        channel_weights: list[list[float]] = [[1,2], [1]],
        channel_merge: str = "mean",
        model_name: str = "cpsam",
        diameter: int = 0,
        normalize: dict = {'percentile': [0.5, 99.5]},
        resize_factor: float = 1,
        mask_name: str = 'cell',
        reduce_mask_name: bool = True,
        overwrite_mask: bool = True,
        flow_threshold: float = 0.4,
        cellprob_threshold: float = 0,
        device: str = None,
        gpu_batch_size: int = 64,
    ):
        """
        Initializes the CellposeSegmentation class.

        Args:
            channel_names: list of channel names, 1st layer list means channels to use, 2nd layer list means merged channel.
            channel_weights: list of channel weights.
            channel_merge: str, method to merge channels ('mean' or 'sum').
            model_name: str, name of the Cellpose model to use.
            diameter: int, diameter of cells in pixels (0 or None for auto-size).
            normalize: dict, normalization settings for the model.
            resize_factor: float, factor to resize the image before segmentation.
            mask_name: str, name to append to the saved mask file.
            reduce_mask_name: bool, if True, simplifies channel name to ch0 in the mask filename.
            overwrite_mask: bool, if True, overwrites existing mask files.
            flow_threshold: float, flow threshold for segmentation.
            cellprob_threshold: float, cell probability threshold for segmentation.
            device: str, device to use for the model ('cpu', 'gpu', or None for auto-detect).
            gpu_batch_size: int, batch size for GPU processing, decrease it when accurring out of memory.
        """
        self.channel_names = channel_names
        self.channel_weights = channel_weights
        self.channel_merge = channel_merge
        self.diameter = diameter
        self.model_name = model_name
        self.normalize = normalize
        self.resize_factor = resize_factor
        self.mask_name = mask_name
        self.reduce_mask_name = reduce_mask_name
        self.overwrite_mask = overwrite_mask
        self.flow_threshold = flow_threshold
        self.cellprob_threshold = cellprob_threshold
        self.gpu_batch_size = gpu_batch_size
        self.device = self._get_device(device)
        self.model = self._initialize_model()

    def _get_device(self, device_input: str) -> torch.device:
        """Determines the device (CPU/GPU) for the model."""
        if device_input is None:
            return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        elif device_input == "gpu":
            return torch.device('cuda')
        elif device_input == "cpu":
            return torch.device('cpu')
        else:
            raise ValueError("device must be 'cpu', 'gpu', or None")

    def _initialize_model(self) -> models.CellposeModel:
        """Initializes the Cellpose model."""
        print('------- CellposeSegmentation --------')
        print(f'segmentation model: {self.model_name}')
        print(f'segmentation channel: {self.channel_names}')
        print(f'segmentation weight: {self.channel_weights}')
        print(f'segmentation mask_name: {self.mask_name}')
        print(f'segmentation device: {self.device}')
        return models.CellposeModel(device=self.device, pretrained_model=self.model_name)
    
    def __repr__(self):
        info = f'CellposeSegmention\nmodel: {self.model_name}\nchannel: {self.channel_names}\nweight: {self.channel_weights}\nmask_name: {self.mask_name}\ndevice: {self.device}'
        return info

    def _read_img_from_dataset(self, df: pd.DataFrame, row_index: int) -> np.ndarray:
        """
        Reads and preprocesses an image from the dataset for a given row.
        """
        assert len(self.channel_names) <= 3, 'channel_names must be less/equal than 3'
        list_flat = sum(self.channel_names, [])
        assert all([x in df.columns.to_list() for x in list_flat]), 'all channel_names must be in df columns'
        assert len(self.channel_names) == len(self.channel_weights), 'different length of channel names and weights'
        
        img_mutli_chan = []
        for idx, channs in enumerate(self.channel_names):
            if len(channs) > 1:
                img_paths = [join(df['directory'].iloc[row_index], x) 
                             for x in df[channs].iloc[row_index].values.flatten().tolist()]
                
                imgs = [imread(j) for j in img_paths]
                
                if self.channel_merge == "mean":
                    norm_weights = [i / sum(self.channel_weights[idx]) for i in self.channel_weights[idx]]
                    img = np.mean([img_data * norm_weights[i] for i, img_data in enumerate(imgs)], axis=0).astype('uint16')
                elif self.channel_merge == "sum":
                    img = np.sum([img_data * self.channel_weights[idx][i] for i, img_data in enumerate(imgs)], axis=0).astype('uint16')
                else:
                    raise ValueError("channel_merge must be one of ('mean', 'sum')")
            else:
                img_path = join(df['directory'].iloc[row_index], df[channs].iloc[row_index].values[0])
                img = imread(img_path)
            
            img_mutli_chan.append(img)
        
        img_mutli_chan = np.stack(img_mutli_chan)

        if self.resize_factor != 1:
            img_mutli_chan = rescale(img_mutli_chan, scale=self.resize_factor, channel_axis=0, order=1)

        return img_mutli_chan


    def run(self, df: pd.DataFrame, test: int = None) -> None:
        """
        Performs segmentation on all images in the dataset.
        Parameters:
        df: pd.DataFrame
        test: int if int > 0, it will only run test mode on N=test segmentation
        """
        diameter_val = None if self.diameter is None or self.diameter == 0 else int(self.diameter * self.resize_factor)
        
        if test is not None and test > 0:
            run_df = df.sample(n=test, random_state=42)
        else:
            run_df = df

        for idx in tqdm(range(len(run_df))):
            try:
                fname = join(run_df['directory'].iloc[idx], run_df.iloc[idx][self.channel_names[0][0]])
                save_stem = join(dirname(fname), f'{splitext(basename(fname))[0]}') + '_cp_masks'
                
                if self.reduce_mask_name:
                    save_stem = re.sub('_ch\\d+_', '_ch0_', save_stem)
                
                save_fname = f'{save_stem}_{self.mask_name}.png'
                
                if not self.overwrite_mask and exists(save_fname):
                    continue

                img = self._read_img_from_dataset(run_df, idx)
                
                mask, flow, style, *_ = self.model.eval(
                    img, 
                    batch_size=self.gpu_batch_size,
                    channel_axis=0,
                    normalize=self.normalize,
                    diameter=diameter_val,
                    flow_threshold=self.flow_threshold,
                    cellprob_threshold=self.cellprob_threshold)

                # Filter if no masks found
                if len(np.unique(mask)) > 1:
                    if self.resize_factor != 1:
                        mask = rescale(mask, 1/self.resize_factor, order=0)
                    cp_io.save_masks(
                        img, closing(mask), flow, file_names=save_stem, suffix=f'_{self.mask_name}')

            except Exception as e:
                print(f"Error processing row {idx}: {e}")



#%%
if __name__ == '__main__':

    from .dataset import ImageDataSet
    dataset = ImageDataSet('/media/hao/Data/Project_MYC/2024-10-31_MYCN_saturation_curve_adding_truncation')

    segmenter = CellposeSegmentation(
        df=dataset.df,
        channel_names=[['ch1'], ['ch2']],
        channel_weights=[[1], [1]],
        model_name="cpsam",
        diameter=0,
        resize_factor=0.5,
        mask_name="cell",
        gpu_batch_size=32)
    segmenter
    # segmenter.run(test=2)
