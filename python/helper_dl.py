#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 23:06:28 2025

@author: hao
"""

import os
import numpy as np
import re
from natsort import natsorted
import pandas as pd
from tqdm import tqdm
import imageio.v3 as iio
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.transform import resize

# from skimage.transform import resize
# from skimage.color import gray2rgb
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from torch.utils.data import random_split
from umap import UMAP

import torch
from torch.utils.data import Dataset
# from torch import nn
# from torchvision import transforms as T
# import pytorch_lightning as pl
from monai.inferers import sliding_window_inference

import pickle
import base64
from io import BytesIO
from PIL import Image

import dash
from dash import dcc, html, Input, Output, State
from dash.dependencies import Input, Output
from flask import Flask
import plotly.express as px


class SingCell(Dataset):
    def __init__(self, 
                 data_dir: str = None,
                 include_pattern: str = None, # only works if data_dir is not None
                 exclude_pattern: str = None, # only works if data_dir is not None
                 image_paths: list[str] = None, # used if data_dir is None
                 transform = None, 
                 normalize: str = "percentile",
                 img_suffix: str = ".tiff",
                 max_pixel_value: int = 65535
                 ):
        super().__init__()
        
        self.normalize = normalize
        self.transform = transform
        self.max_pixel_value = max_pixel_value

        # recursively find all TIFF images within the root directory
        if image_paths is None and data_dir is not None:
            image_paths = []
            for dirpath, dirnames, filenames in os.walk(data_dir, followlinks=False):
                for filename in filenames:
                    if filename.endswith(img_suffix):
                        image_paths.append(os.path.join(dirpath, filename))
            if include_pattern is not None:
                print(f"include files by {include_pattern}")
                image_paths = [f for f in image_paths if re.match(include_pattern, f)]
            if exclude_pattern is not None:
                print(f"exclude files by {exclude_pattern}")
                image_paths = [f for f in image_paths if not re.match(exclude_pattern, f)]
            if not image_paths:
                raise FileNotFoundError(f"No TIFF images found in directory: {data_dir}")
                
        self.image_paths = image_paths
        print(f"Dataset initialized with {len(image_paths)} images")
        
        # get image parent dirname as class label
        self.labels = [os.path.basename(os.path.dirname(f)) for f in image_paths]
        uni_labels = natsorted(list(set(self.labels)))
        label_to_idx = {l: idx for idx, l in enumerate(uni_labels)}
        self.idxs = [label_to_idx[l] for l in self.labels]
        self.idx_to_label = {idx: l for l, idx in label_to_idx.items()}
        print(f"Total {len(self.idx_to_label)} parent directory categories found")
        
    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # # singele channel
        # image_path = "/media/hao/Data1/single_cell_dataset/single_channel_128/myc_condensate/60x-293t-24h-96w-100ng-300ng-656-130-149-150__1/p1_c1_ch1.tiff"
        # # multichannel
        # image_path = "/media/hao/Data/Project_barcode/2024-12-16_BC_test/cropped_cell/BC100/M8__B16/p2_c1.tiff"
        
        image_path = self.image_paths[idx]
        img = iio.imread(image_path) # CHW for 2-channel
        # print(img.shape)

        # # expand channel dim
        # img = img[None,:,:] if img.ndim == 2 else img
        # # resize
        # img = resize(img, (img.shape[0], self.size, self.size), preserve_range=True)
        
        # scale to [0, 1]
        img = (img / self.max_pixel_value).astype("float32")
        
        # normalized image
        if self.normalize:
            img = normalize_img(img, method=self.normalize, channel_index=0)
            # print(img.shape)
                
        # apply transforms if any (currently defined on numpy array)
        if self.transform:
            # albumentations method, albumentations accept HWC, toTensor output CHW
            img = self.transform(image=img.transpose(1,2,0))['image']
            # print(img.shape)
        
            
        return [img, self.idxs[idx], image_path]
    


class SingCellInstSeg(Dataset):
    def __init__(self, 
                 data_dir: str = None,
                 include_pattern: str = None, # only works if data_dir is not None
                 exclude_pattern: str = None, # only works if data_dir is not None
                 image_paths: list[str] = None, # used if data_dir is None
                 transform = None, 
                 normalize: str = "percentile",
                 img_suffix: str = ".tiff",
                 max_pixel_value: int = 65535
                 ):
        super().__init__()
        
        self.normalize = normalize
        self.transform = transform
        self.max_pixel_value = max_pixel_value

        # recursively find all TIFF images within the root directory
        if image_paths is None and data_dir is not None:
            image_paths = []
            for dirpath, dirnames, filenames in os.walk(data_dir, followlinks=False):
                for filename in filenames:
                    if filename.endswith(img_suffix):
                        image_paths.append(os.path.join(dirpath, filename))
            if include_pattern is not None:
                print(f"include files by {include_pattern}")
                image_paths = [f for f in image_paths if re.match(include_pattern, f)]
            if exclude_pattern is not None:
                print(f"exclude files by {exclude_pattern}")
                image_paths = [f for f in image_paths if not re.match(exclude_pattern, f)]
            if not image_paths:
                raise FileNotFoundError(f"No TIFF images found in directory: {data_dir}")
                
        self.image_paths = image_paths
        print(f"Dataset initialized with {len(image_paths)} images")
        
        # get image parent dirname as class label
        self.labels = [os.path.basename(os.path.dirname(f)) for f in image_paths]
        uni_labels = natsorted(list(set(self.labels)))
        label_to_idx = {l: idx for idx, l in enumerate(uni_labels)}
        self.idxs = [label_to_idx[l] for l in self.labels]
        self.idx_to_label = {idx: l for l, idx in label_to_idx.items()}
        print(f"Total {len(self.idx_to_label)} parent directory categories found")
        
    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # # singele channel
        # image_path = "/media/hao/Data1/single_cell_dataset/single_channel_128/myc_condensate/60x-293t-24h-96w-100ng-300ng-656-130-149-150__1/p1_c1_ch1.tiff"
        # # multichannel
        # image_path = "/media/hao/Data/Project_barcode/2024-12-16_BC_test/cropped_cell/BC100/M8__B16/p2_c1.tiff"
        
        image_path = self.image_paths[idx]
        img = iio.imread(image_path) # CHW for 2-channel
        # print(img.shape)

        # # expand channel dim
        # img = img[None,:,:] if img.ndim == 2 else img
        # # resize
        # img = resize(img, (img.shape[0], self.size, self.size), preserve_range=True)
        
        # scale to [0, 1]
        img = (img / self.max_pixel_value).astype("float32")
        
        # normalized image
        if self.normalize:
            img = normalize_img(img, method=self.normalize, channel_index=0)
            # print(img.shape)
                
        # apply transforms if any (currently defined on numpy array)
        if self.transform:
            # albumentations method, albumentations accept HWC, toTensor output CHW
            img = self.transform(image=img.transpose(1,2,0))['image']
            # print(img.shape)
        
            
        return [img, self.idxs[idx], image_path]


    
class OpenCellDataSet(Dataset):
    def __init__(self, 
                 data_dir: str = None, 
                 size = 300,
                 transform = None, 
                 normalize: str = "percentile",
                 img_suffix = ".tiff",
                 max_pixel_value: int = 65535):
        super().__init__()
        
        self.data_dir = data_dir
        self.size = size
        self.normalize = normalize
        self.transform = transform
        self.max_pixel_value = max_pixel_value

        # recursively find all TIFF images within the root directory
        image_paths = []
        for dirpath, dirnames, filenames in os.walk(data_dir, followlinks=False):
            for filename in filenames:
                if filename.endswith(img_suffix):
                    image_paths.append(os.path.join(dirpath, filename))
        if not image_paths:
            raise FileNotFoundError(f"No TIFF images found in directory: {data_dir}")
        self.image_paths = image_paths
        print(f"Dataset initialized with {len(image_paths)} images")
        
        # get image parent dirname as class label
        self.labels = [os.path.basename(os.path.dirname(f)) for f in image_paths]
        uni_labels = natsorted(list(set(self.labels)))
        label_to_idx = {l: idx for idx, l in enumerate(uni_labels)}
        self.idxs = [label_to_idx[l] for l in self.labels]
        self.idx_to_label = {idx: l for l, idx in label_to_idx.items()}
        print(f"Total {len(self.idx_to_label)} parent directory categories found")
        
    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # # singele channel
        # image_path = "/media/hao/Data1/single_cell_dataset/single_channel_128/myc_condensate/60x-293t-24h-96w-100ng-300ng-656-130-149-150__1/p1_c1_ch1.tiff"
        # # multichannel
        # image_path = "/media/hao/Data/Project_barcode/2024-12-16_BC_test/cropped_cell/BC100/M8__B16/p2_c1.tiff"
        
        image_path = self.image_paths[idx]
        img = iio.imread(image_path)
        # expand channel dim
        img = img[None,:,:] if img.ndim == 2 else img
        # resize
        img = resize(img, (img.shape[0], self.size, self.size), preserve_range=True)
        # scale to [0, 1]
        img = (img / self.max_pixel_value).astype("float32")
        # plt.imshow(image[0,:,:])
        
        # normalized image
        if self.normalize:
            img = normalize_img(img, method=self.normalize)
                   
        # to torch
        img = torch.from_numpy(img)
        
        # apply transforms if any (currently defined on numpy array)
        if self.transform:
            img = self.transform(img)

        return [img, self.idxs[idx], image_path]
    
    
class HPASingleCell(Dataset):
    def __init__(self, 
                 data_dir: str = None, 
                 size = 128,
                 normalize: str = None,
                 transform = None, 
                 suffix = ".tiff",
                 max_pixel_value: int = 255):
        super().__init__()
        
        self.data_dir = data_dir
        self.size = size
        self.normalize = normalize
        self.transform = transform
        self.max_pixel_value = max_pixel_value

        # recursively find all TIFF images within the root directory
        image_paths = []
        for dirpath, dirnames, filenames in os.walk(data_dir, followlinks=False):
            for filename in filenames:
                if filename.endswith(suffix):
                    image_paths.append(os.path.join(dirpath, filename))
        if not image_paths:
            raise FileNotFoundError(f"No TIFF images found in directory: {data_dir}")
        self.image_paths = image_paths
        print(f"Dataset initialized with {len(image_paths)} images")
        
        # get image parent dirname as class label
        self.labels = [os.path.basename(os.path.dirname(f)) for f in image_paths]
        uni_labels = natsorted(list(set(self.labels)))
        label_to_idx = {l: idx for idx, l in enumerate(uni_labels)}
        self.idxs = [label_to_idx[l] for l in self.labels]
        self.idx_to_label = {idx: l for l, idx in label_to_idx.items()}
        print(f"Total {len(self.idx_to_label)} parent directory categories found")
        
    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # # singele channel
        # image_path = "/media/hao/Data1/single_cell_dataset/single_channel_128/myc_condensate/60x-293t-24h-96w-100ng-300ng-656-130-149-150__1/p1_c1_ch1.tiff"
        # # multichannel
        # image_path = "/media/hao/Data/Project_barcode/2024-12-16_BC_test/cropped_cell/BC100/M8__B16/p2_c1.tiff"
        
        image_path = self.image_paths[idx]
        img = iio.imread(image_path)
        # only keep 3rd green channel for training
        img = img[None,2,:,:]
        img = resize(img, (img.shape[0], self.size, self.size), preserve_range=True)
        # scale to [0, 1]
        img = (img / self.max_pixel_value)
        # plt.imshow(img)
        
        # normalized img
        if self.normalize:
            img = normalize_img(img, method=self.normalize)
                   
        # to torch
        img = torch.from_numpy(img)
        
        # apply transforms if any (currently defined on numpy array)
        if self.transform:
            img = self.transform(img)

        return [img, self.idxs[idx], image_path]




def normalize_img(img: np.array,
                  method: str = "percentile",
                  percentile_value: tuple = (0.1, 99.9),
                  channel_index: int = None):
    """
    Normalizes an image, optionally by channel.

    Args:
        img: The input NumPy array representing the image.
        method: The normalization method ("scale", "percentile", "minmax", or None).
        percentile_value: The percentile values for the "percentile" method.
        channel_index: The index of the channel to normalize, or None to normalize the entire image.
    """
    if method not in ["scale", "percentile", "minmax", None]:
        quit('method should one of ["scale","percentile","minmax",None]')

    if channel_index is not None:
        if channel_index >= img.shape[-1] or channel_index < 0:
            raise ValueError(f"channel_index {channel_index} is out of range for image with shape {img.shape}")
        
        channel_img = img[..., channel_index]
        mask = channel_img > 0
        nonzero_values = channel_img[mask]
        
        if len(nonzero_values) == 0:
            return img  # If all zeros in channel, return original
        
        if method == "scale":
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
                img_norm_channel = (channel_img - p_low) / (p_high - p_low)
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

        if method == "scale":
            mean = np.mean(nonzero_values)
            std = np.std(nonzero_values)
            if std == 0: std = 1.0
            img_norm = (img - mean) / std

        elif method == "percentile":
            p_low = np.percentile(nonzero_values, percentile_value[0])
            p_high = np.percentile(nonzero_values, percentile_value[1])
            if p_high - p_low > 0:
                img_norm = (img - p_low) / (p_high - p_low)
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
    
    
    # Keep zero regions unchanged
    img_norm[~mask] = 0  

    return img_norm



# from sklearn.model_selection import train_test_split
def split_dataset(dataset, 
                  train_ratio=0.8, 
                  val_ratio=0.2, 
                  test_ratio=0, 
                  random_seed=None):
    """Splits a dataset into training, validation, and test sets.
    
    Args:
        dataset: The dataset to split.
        train_ratio: The ratio of the training set.
        val_ratio: The ratio of the validation set.
        test_ratio: The ratio of the test set.
        random_seed: Random seed for reproducibility.
    
    Returns:
        A tuple of three datasets: (train_dataset, val_dataset, test_dataset).
    """
    
    if sum([train_ratio, val_ratio, test_ratio]) != 1:
        raise ValueError("train_ratio, val_ratio, and test_ratio must sum to 1")

    total_length = len(dataset)
    train_length = int(train_ratio * total_length)
    val_length = int(val_ratio * total_length)
    test_length = total_length - train_length - val_length

    # Ensure the splits add up to the total length
    lengths = [train_length, val_length, test_length]

    # Use a generator for reproducibility
    generator = torch.Generator()
    if random_seed is not None:
        generator.manual_seed(random_seed)
    
    dataset_train, dataset_val, dataset_test = random_split(dataset, lengths, generator)
    
    dataset_val.dataset.transform = None
    dataset_test.dataset.transform = None
    
    return dataset_train, dataset_val, dataset_test



def visualize_dataloader(batch: torch.Tensor, 
                         labels: torch.Tensor = None, 
                         num_images: int = 4, 
                         cmap: str = "jet", 
                         show_outline: bool = True, 
                         figsize_scale: float = 2, 
                         colorbar_mode: str = "each", 
                         vmin: float = None, 
                         vmax: float = None):
    """
    Visualizes samples from a PyTorch dataloader with optional labels and color bars.
    
    Args:
        batch (torch.Tensor): Input tensor of shape (B, C, H, W)
        labels (torch.Tensor, optional): Labels corresponding to each batch sample
        num_images (int): Number of images to display from the batch
        cmap (str): Colormap to use for displaying images (default: "jet")
        show_outline (bool): Whether to show the subplot rectangle outline (default: True)
        figsize_scale (float): Scale factor for figure size (default: 3)
        colorbar_mode (str): Determines color bar usage - "each" (default), "row", or "column"
        vmin (float, optional): Minimum value for color scaling across all images
        vmax (float, optional): Maximum value for color scaling across all images
    """
    B, C, H, W = batch.shape
    num_images = min(B, num_images)  # Ensure we don't exceed batch size
    
    fig, axes = plt.subplots(num_images, C, 
                             figsize=(C * figsize_scale, num_images * figsize_scale), 
                             gridspec_kw={'hspace': 0.5, 'wspace': 0.5})
    if C == 1:
        axes = axes[:, None]  # Ensure 2D indexing for single channel
    elif num_images == 1:
        axes = axes[None, :]  # Ensure 2D indexing for single row batch
    
    vmin_values, vmax_values = {}, {}
    
    if colorbar_mode == "row":
        for i in range(num_images):
            nonzero_data = batch[i].cpu().numpy()
            nonzero_data = nonzero_data[nonzero_data > 0]
            if nonzero_data.size > 0:
                vmin_values[i] = vmin if vmin is not None else np.percentile(nonzero_data, 0.05)
                vmax_values[i] = vmax if vmax is not None else np.percentile(nonzero_data, 99.5)
            else:
                vmin_values[i], vmax_values[i] = 0, 1
    elif colorbar_mode == "column":
        for j in range(C):
            nonzero_data = batch[:, j].cpu().numpy()
            nonzero_data = nonzero_data[nonzero_data > 0]
            if nonzero_data.size > 0:
                vmin_values[j] = vmin if vmin is not None else np.percentile(nonzero_data, 0.05)
                vmax_values[j] = vmax if vmax is not None else np.percentile(nonzero_data, 99.5)
            else:
                vmin_values[j], vmax_values[j] = 0, 1
    
    for j in range(C):
        for i in range(num_images):
            ax = axes[i, j]
            img = batch[i, j].cpu().numpy()
            img[img == 0] = np.nan  # Omit zero pixels
            
            if colorbar_mode == "each":
                nonzero_data = img[img > 0]
                if nonzero_data.size > 0:
                    vmin_calc = np.percentile(nonzero_data, 0.05)
                    vmax_calc = np.percentile(nonzero_data, 99.5)
                    vmin_final = vmin if vmin is not None else vmin_calc
                    vmax_final = vmax if vmax is not None else vmax_calc
                else:
                    vmin_final, vmax_final = 0, 1
            elif colorbar_mode == "row":
                vmin_final, vmax_final = vmin_values[i], vmax_values[i]
            elif colorbar_mode == "column":
                vmin_final, vmax_final = vmin_values[j], vmax_values[j]
            
            im = ax.imshow(img, cmap=cmap, vmin=vmin_final, vmax=vmax_final)
            
            if not show_outline:
                ax.set_xticks([])
                ax.set_yticks([])
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
            else:
                ax.axis('on')
            
            if labels is not None and j == 0:
                ax.set_title(f"Label: {labels[i].item()}")
            
            # Add color bar next to each image
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    plt.show()
# Example usage:
# Assuming `dataloader` is a PyTorch DataLoader that returns (images, labels)
# for images, labels in dataloader:
#     visualize_dataloader_samples(images, labels)
#     break

    

def embed_umap(array, 
               color=None, 
               plot=False, 
               scale=True, 
               dim_thres_for_pca=256, 
               pca_components=100, 
               **kwargs):
    """
    Perform 2D dimensionality reduction using UMAP and optionally visualize results.

    Parameters:
        array (ndarray): Input data where rows are samples and columns are features.
        scale (bool): Whether to standardize features before applying UMAP.
        dim_thres_for_pca (int): If feature dimension exceeds this, apply PCA first.
        pca_components (int): Number of components to keep in PCA (if applied).
        color (array-like or None): Optional labels for coloring points in the scatter plot.
        plot (bool): If True, displays the UMAP visualization.
        **kwargs: Additional parameters for UMAP.

    Returns:
        DataFrame: Contains UMAP dimensions and labels (if provided).
    """
    if torch.is_tensor(array):
        array = array.detach().numpy()
        
    if np.any(np.isnan(array)):
        raise ValueError("Input array contains NaN values. Please handle missing values before processing.")

    if scale:
        array = StandardScaler().fit_transform(array)

    # Apply PCA if dimensionality is too high
    if array.shape[1] > dim_thres_for_pca:
        array = PCA(n_components=min(pca_components, array.shape[1])).fit_transform(array)

    # Perform UMAP
    umap_results = UMAP(n_components=2, random_state=42, **kwargs).fit_transform(array)

    # Create DataFrame with UMAP results
    df = pd.DataFrame(umap_results, columns=["UMAP_1", "UMAP_2"])
    
    if color is not None:
        df["Label"] = color

    # Plot if requested
    if plot:
        plt.figure(figsize=(8, 6))
        if color is not None:
            sns.scatterplot(x=df["UMAP_1"], y=df["UMAP_2"], hue=df["Label"], palette="tab10", alpha=0.75)
            plt.legend(title="Label")
        else:
            plt.scatter(df["UMAP_1"], df["UMAP_2"], alpha=0.75)
        plt.title("UMAP Projection")
        plt.xlabel("UMAP_1")
        plt.ylabel("UMAP_2")
        plt.show()

    return df



def extract_feature_timm(
    model, # timm model, last layer removed
    dataloader, # torch dataloader,
    device: torch.device = torch.device("cuda"),
    output_umap_enbed: bool = False,
    plot: bool = False, # display umap
    output_encode_img = False
    ):
    """
    extract feature using trained timm model
    """
    
    # prepare model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    model.eval()
    
    feats = []
    paths = []
    labels = []
    img_encoder = []
    
    with torch.no_grad():
        for img, label, img_path in tqdm(dataloader):
            # feature
            feat = model(img.to(device)).cpu().numpy()
            feats.append(feat)
            
            # path
            paths.extend(list(img_path))
            # label
            labels.extend(label.cpu().numpy().tolist())
            
            # img encoded
            if output_encode_img:
                img_pil = img.cpu().numpy() # BCHW
                
                for b in img_pil:
                    # to PIL format
                    b = b.transpose(1, 2, 0)  # CHW → HWC
                    if b.dtype != np.uint8:
                        b = (b * 255).clip(0, 255).astype(np.uint8)
                    if b.shape[-1] == 1:
                        b = b.squeeze(-1) 
                    b = Image.fromarray(b)
                    
                    # Encode to base64
                    buffered = BytesIO()
                    b.save(buffered, format="PNG")
                    encoded_string = base64.b64encode(buffered.getvalue()).decode("utf-8")
                    img_encoder.append(f'data:image/png;base64,{encoded_string}')
    
    feats = np.vstack(feats)
    
    if output_umap_enbed:
        umap_embed = embed_umap(
            feats, color=labels, plot=plot, scale=True, 
            dim_thres_for_pca=256, pca_components=100,
            n_neighbors=5, min_dist=0.01)
    
    # output df
    df = pd.DataFrame()
    df['Path'] = paths
    df['Label'] = labels
    
    if output_umap_enbed:
        df[['UMAP_1', 'UMAP_2']] = umap_embed[['UMAP_1', 'UMAP_2']]
        
    if output_encode_img:
        df['Image'] = img_encoder
    
    return [feats, df]



def extract_feature_semantic(
    model, # semantic_segmentation_pytorch model
    dataloader, # torch dataloader
    device: torch.device = torch.device("cuda"),
    output_umap_enbed: bool = False,
    plot: bool = False, # display umap
    output_encode_img = False
    ):
    """
    extract feature using trained model
    """
    
    # # prepare model
    # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # model.to(device)
    # model.eval()
    
    feats = []
    paths = []
    labels = []
    img_encoder = []
    
    with torch.no_grad():
        for img, label, img_path in tqdm(dataloader):
            # feature
            feat = model(img.to(device)).cpu().numpy()
            feats.append(feat)
            
            # path
            paths.extend(list(img_path))
            # label
            labels.extend(label.cpu().numpy().tolist())
            
            # img encoded
            if output_encode_img:
                img_pil = img.cpu().numpy() # BCHW
                
                for b in img_pil:
                    # to PIL format
                    b = b.transpose(1, 2, 0)  # CHW → HWC
                    if b.dtype != np.uint8:
                        b = (b * 255).clip(0, 255).astype(np.uint8)
                    if b.shape[-1] == 1:
                        b = b.squeeze(-1) 
                    b = Image.fromarray(b)
                    
                    # Encode to base64
                    buffered = BytesIO()
                    b.save(buffered, format="PNG")
                    encoded_string = base64.b64encode(buffered.getvalue()).decode("utf-8")
                    img_encoder.append(f'data:image/png;base64,{encoded_string}')
    
    feats = np.vstack(feats)
    
    if output_umap_enbed:
        umap_embed = embed_umap(
            feats, color=labels, plot=plot, scale=True, 
            dim_thres_for_pca=256, pca_components=100,
            n_neighbors=5, min_dist=0.01)
    
    # output df
    df = pd.DataFrame()
    df['Path'] = paths
    df['Label'] = labels
    
    if output_umap_enbed:
        df[['UMAP_1', 'UMAP_2']] = umap_embed[['UMAP_1', 'UMAP_2']]
        
    if output_encode_img:
        df['Image'] = img_encoder
    
    return [feats, df]



def vis_feature(data_path, plot_size=800, image_size=200):
    
    # Load data
    with open(data_path, 'rb') as file:
        _, df = pickle.load(file)

    # Initialize Dash App
    server = Flask(__name__)
    app = dash.Dash(__name__, server=server)

    app.layout = html.Div([
        html.H1("Interactive UMAP Visualization"),
        
        # Scatter plot with image overlay container
        html.Div([
            dcc.Graph(id="umap-scatter", config={"displayModeBar": False}),
            html.Img(id="hover-image", style={
                "position": "absolute",
                "width": f"{image_size}px", # image size
                # "height": "150px", # image size
                "display": "none",
                "border": "2px solid black",
                "border-radius": "5px",
                "box-shadow": "2px 2px 10px rgba(0,0,0,0.5)",
                "right": "20px",  # Move image to right side
                "top": "20px"    # Move image to top side
            })
        ], style={"position": "relative", "width": f"{plot_size}px", "margin": "auto"}),
    ], style={"textAlign": "center"})

    # Callback to update the scatter plot (forces square shape)
    @app.callback(
        Output("umap-scatter", "figure"),
        Input("umap-scatter", "hoverData")
        )
    def update_scatter_plot(hoverData):
        fig = px.scatter(
            df, x="UMAP_1", y="UMAP_2",
            hover_name=df["Label"],
            title="Hover over a point to see the image",
        )
        
        # Force square aspect ratio
        fig.update_layout(
            width=plot_size, height=plot_size,
            xaxis=dict(scaleanchor="y"),  # Ensures equal scale for both axes
            yaxis=dict(scaleanchor="x")
        )

        return fig

    # Callback to update hover image position and display
    @app.callback(
        [Output("hover-image", "src"),
          Output("hover-image", "style")],
        [Input("umap-scatter", "hoverData")]
    )
    def update_hover_image(hoverData):
        if hoverData is None:
            return "", {"display": "none"}

        # Extract hover data
        point_index = hoverData["points"][0]["pointIndex"]
        img_src = df.iloc[point_index]["Image"]

        img_style = {
            "position": "absolute",
            "right": "800px",  
            "top": "300px",  
            "width": f"{image_size}px",  
            # "height": f"{image_size}px",
            "display": "block",
            "border": "2px solid black",
            "border-radius": "5px",
            "box-shadow": "2px 2px 10px rgba(0,0,0,0.5)",
        }

        return img_src, img_style

    # Run Dash App
    app.run_server(debug=True)



def sliding_window_inference(
        model, 
        image, 
        roi_size=512, 
        sw_batch_size=16, 
        overlap=0.0,
        device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
        ):
    """
    Perform sliding window inference using MONAI's implementation
    
    Args:
        model: Trained PyTorch model
        image: Input image as numpy array (H,W) 
        roi_size: Size of patches to process
        sw_batch_size: Batch size for inference
        overlap: Overlap fraction between patches
    
    Returns:
        Full size prediction mask
    """

    # model.to(device)
    # model.eval()
    
    # Prepare image
    assert image.ndim == 2, f"Image has {image.ndim} dimensions, expected 2"
    # if image.ndim == 2:
        # image = (image / 65535).astype(np.float32)
        # # Add batch and channel dimensions
    image = torch.from_numpy(image).unsqueeze(0).unsqueeze(0)
    image = image.to(device)
    
    with torch.no_grad():
        # MONAI expects input shape: (B,C,H,W)
        # roi_size needs to be a sequence for 2D: (H,W)
        outputs = sliding_window_inference(
            inputs=image,
            roi_size=(roi_size, roi_size),
            sw_batch_size=sw_batch_size,
            predictor=model,
            overlap=overlap,
            mode="constant",
            padding_mode="constant")
    
    # Convert back to numpy array and remove batch/channel dims
    prediction = outputs.cpu().numpy()[0, 0]
    
    return prediction
