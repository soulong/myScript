#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 22:44:18 2025

@author: hao
"""

import os
from os.path import basename, dirname, exists, join, splitext
from sys import exit
from subprocess import Popen, PIPE, DEVNULL
from multiprocessing import Process
import re
from glob import glob
import time
import random
from random import sample
import numpy as np
import pandas as pd
from natsort import natsorted
import nd2
import cv2
from tifffile import imread, imwrite
import imageio.v3 as iio
from tqdm import tqdm
import torch
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
import albumentations as A
from albumentations.pytorch import ToTensorV2
from cellpose import models
from cellpose import io as cp_io
import math
from skimage.morphology import closing, disk
from skimage.measure import regionprops, regionprops_table
from skimage.transform import resize, rescale
from skimage.segmentation import clear_border
# from scipy.ndimage import rotate
import scipy.ndimage as ndi
import sqlite3
import gc



#%% functions
def timeit(f):
    def timed(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        print('func:%r \ntook: %2.4f sec' % (f.__name__, te-ts))
        return result
    return timed


def fix_filepath(paths: list[str] = None) -> list[str]:
    """
    Fix file path, remove " " and "." in file name
    """
    if isinstance(paths, str): paths = [paths]
    
    fixed_paths = []
    for path in paths:
        if not exists(path): 
            print(f'path not exist, skip it \n {path}')
        
        if " " in path or "." in path:
            base_dir = dirname(path)
            base_name = basename(path)
            base_name_fixed = base_name.replace(" ", "-")
            base_name_fixed = base_name_fixed.replace(".", "-")
            base_name_fixed = base_name_fixed.replace("-nd2", ".nd2")
            path_fixed = join(base_dir, base_name_fixed)
            fixed_paths.append(path_fixed)
            
    return fixed_paths



def set_img_channel(img: np.array = None, # 
                    order='HWC') -> np.ndarray: # 3D or 4D array [[B],order]
    """
    Set image channel order
    Args:
        img: np.ndarray, image array,  accept 2D: [HW], 3D: [HWC,CHW], 4D: [BCHW, BHWC]
        order: str, order of the image channel, accept [HWC,CHW]
    Returns:
        np.ndarray, image array with the specified channel order
    """
    assert isinstance(img, np.ndarray), 'img must be numpy array'

    if img.ndim == 2:
        if order == 'CHW':
            img = img[None,:,:]
        else:
            img = img[:,:, None]
            
    elif img.ndim == 3:
        if order == 'CHW' and np.argmin(img.shape) == 2:
            img = img.transpose(1,2,0)
        elif order == 'HWC' and np.argmin(img.shape) == 0:
            img = img.transpose(1,2,0)

    elif img.ndim == 4:
        if order == 'CHW' and img.shape[1] > img.shape[2]: # bhwc -> bchw
            img = img.transpose(0,3,1,2)
        elif order == 'HWC' and img.shape[1] < img.shape[2]: # bchw -> bhwc
            img = img.transpose(0,2,3,1)

    else:
        exit('img dim must be 2 or 3, or img size < channel size')
        
    return img
# x = np.random.random([8,500,400])
# set_img_channel(x, order='HWC').shape


def nd2_to_tiff(data_dir: str,
                in_dir: str =".",
                out_dir: str ="images",
                resize: list[int, int] = None,
                crop_N: int = None,
                only_meta: bool = False
                ) -> None:
    """
    Convert nd2 files to tiff files
    Args:
        data_dir: str, path to the data directory, this is working directory
        in_dir: str, path to the input directory, search for nd2 files in this directory
        out_dir: str, path to the output directory, save the tiff files in this directory
        resize: list[int, int], resize the image to [H,W]
        crop_N: int, split the image to NxN small separated images
        only_meta: bool, only write metadata
    """
    if not exists(data_dir):
        exit(f"{data_dir} not found, check it")

    dir_in = join(data_dir, in_dir)
    dir_out = join(data_dir, out_dir)
    
    # if exists(dir_out): 
    #     exit(f"{dir_out} exist, skip")
    os.makedirs(dir_out, exist_ok=True)
    
    nd2_files = glob(dir_in + "/**/*.nd2", recursive=True)
    nd2_files = natsorted(nd2_files)
    if len(nd2_files) == 0:
        exit("no nd2 files found, check it")

    # rename nd2 files to remove " ", in case of some unintended error
    # f_name = "V1-MCHERRY-MIFP-EBFP2-100-2.nd2"
    # f_path = join(dir_in, f_name)
    print("check and correct file name")
    fix_filepath(nd2_files)

    # extract position metadata
    print("extract and merge position metadata")
    position_metadata_f = join(data_dir, "position.csv")
    position_list = []
    for f_path in nd2_files:
        # pass
        # get metadata
        with nd2.ND2File(f_path) as file:
            position = pd.DataFrame(file.events())
            # to avoid index probelm
            position["P Index"] = position["P Index"] + 1
            if "T Index" in position.columns:
                position["T Index"] = position["T Index"] + 1
        # add filename prefix
        fname = splitext(basename(f_path))[0]
        # print(fname)
        position.insert(0, "prefix", fname)
        position_list.append(position)
    # merge to one file
    position_df = pd.concat(position_list, axis=0)
    position_df['Time [s]'] = round(position_df['Time [s]'])
    # save
    if exists(position_metadata_f):
        position_df_old = pd.read_csv(position_metadata_f)
        position_df = pd.concat([position_df_old, position_df], ignore_index=True)
        position_df = position_df.drop_duplicates(
            subset=['prefix','Time [s]','P Index'], keep='first')
        position_df.to_csv(position_metadata_f, index=False)
    else:
        position_df.to_csv(position_metadata_f, index=False)

    # extract other metadata
    print("extract channel/resolution metadata")
    metadata_f = join(data_dir, "metadata.txt")
    for f_path in nd2_files:
        try:
            with nd2.ND2File(f_path) as file:
                info = pd.DataFrame(file.frame_metadata(0).channels)
                string = f'\n\n{f_path}:' + f"""
                    -> Objective: {info["microscope"].iloc[0]["objectiveMagnification"]}
                    -> Resize image: {resize}
                    """
                for c_idx in range(len(info["channel"])):
                    if resize is None:
                        string = string + f"""
                        -> Channel: {int(info["channel"].iloc[c_idx]["index"]) + 1}
                        Channel Name: {info["channel"].iloc[c_idx]["name"]}
                        Axis X Pixels: {info["volume"].iloc[0]["voxelCount"][0]}
                        Axis X Resolution per Pixel: {info["volume"].iloc[0]["axesCalibration"][0]} \n"""
                    else:
                        string = string + f"""
                        -> Channel: {int(info["channel"].iloc[c_idx]["index"]) + 1}
                        Channel Name: {info["channel"].iloc[c_idx]["name"]}
                        Axis X Pixels: {resize[0]}
                        Axis X Resolution per Pixel: {info["volume"].iloc[0]["axesCalibration"][0] * (info["volume"].iloc[0]["voxelCount"][0] / resize[0])} \n"""
                # print(string)
                with open(metadata_f, "a") as f:
                    f.write(string)
        except Exception as e:
            print(e)
            print(f"error in extracting {f_path}, skip metadata pharsing")
            continue

    # lopp process
    if not only_meta:
        for f_path in nd2_files:
            fname = splitext(basename(f_path))[0]
            print(f">> {fname}")
            with nd2.ND2File(f_path) as file:
                print(file.sizes)
                if resize is not None:
                    print(f"resized to {resize}")
                # > z-stack, using max projection
                max_projection = True if "Z" in list(file.sizes.keys()) else False
        
                # read data
                dat = file.asarray()
        
                # full recover array to [T, P, Z, C, Y, X]
                if "C" not in list(file.sizes.keys()):
                    dat = np.expand_dims(dat, -3)
                if "Z" not in list(file.sizes.keys()):
                    dat = np.expand_dims(dat, -4)
                if "P" not in list(file.sizes.keys()):
                    dat = np.expand_dims(dat, -5)
                if "T" not in list(file.sizes.keys()):
                    dat = np.expand_dims(dat, -6)
        
                # save to indivual tiff file
                for ti, t in enumerate(dat):
                    for pi, p in enumerate(tqdm(t, desc=f"T{ti+1}: >>P")):
                        # max projection
                        if max_projection:
                            p = np.max(p, axis=0)
                            p = np.expand_dims(p, 0)
                        for zi, z in enumerate(p):
                            for ci, c in enumerate(z):
                                # resize
                                if resize is not None:
                                    c = cv2.resize(c, dsize=resize, interpolation=cv2.INTER_CUBIC)
                                # split into pieces
                                if crop_N is None:
                                    # save
                                    imwrite(join(dir_out, f"{fname}__t{ti+1}_p{pi+1}_z{zi+1}_ch{ci+1}.tiff"), c, compression="zlib")
                                else:
                                    M = c.shape[0]//crop_N
                                    N = c.shape[1]//crop_N
                                    tiles = [c[x:x+M, y:y+N] for x in range(0, c.shape[0], M) for y in range(0, c.shape[1], N)]
                                    for xi, x in enumerate(tiles):
                                        imwrite(join(dir_out, f"{fname}__t{ti+1}_p{pi+1}000{xi}_z{zi+1}_ch{ci+1}.tiff"), x, compression="zlib")



def images_to_dataset(
        data_dir: str = None,
        subset_pattern: str = None,
        image_subdir: str = "images",
        image_suffix: str = '.tiff',
        mask_suffix: str = '.png',  # None for no use
        remove_na_row: bool = True,
        cellprofiler_style: bool = False,
        position_metadata: str = "position.csv",  # None for no usage
        pos_index_col: str = "P Index", # which col in position_metadata to use to merge with df
        image_extractor: str = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_ch?(?P<channel>[0-9]{1,})(?P<self_generated>.*)',
        mask_extractor: str = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_ch?(?P<channel>[0-9]{1,})_(cp_masks_)?(?P<mask_name>.*)'
) -> dict:
    """
    Convert images to dataset
    Args:
        data_dir: str, path to the working directory
        subset_pattern: str, pattern to subset the images
        image_subdir: str, path to the image directory
        image_suffix: str, suffix of the image file
        mask_suffix: str, suffix of the mask file
        remove_na_row: bool, remove the rows with any NA values
        cellprofiler_style: bool, whether to convert to cellprofiler style
        position_metadata: str, path to the position metadata file
    Returns:
        dict, {'df': pd.DataFrame, 'image_colnames': list, 'mask_colnames': list}
    """

    image_extract_cols = re.findall(r'\?P<([^>]+)>', image_extractor)
    if not ('channel' in image_extract_cols or 'self_generated' in image_extract_cols):
        exit("image_extractor must contain [self_generated, channel]")
    mask_extract_cols = re.findall(r'\?P<([^>]+)>', mask_extractor)
    if not 'mask_name' in mask_extract_cols:
        exit("mask_extractor must contain [mask_name]")

    # data_dir = "/media/hao/Data/2024-04-19_barcode_test"
    if not exists(data_dir):
        return None

    print("\n>>>>> processing: " + data_dir)

    # intensity files -------------------------------------
    image = glob(join(data_dir, image_subdir, '**/*') + image_suffix, recursive=True)
    
    # # remove common output files
    # image = [f for f in image if not re.search("Probabilities", f)]
    # image = [f for f in image if not re.search("_pred", f)]
    
    # filter images
    if subset_pattern is not None:
        image = [f for f in image if re.search(subset_pattern, f)]
    
    if len(image) == 0:
        print("no valid images found, check directory or subset_pattern")
        return None

    # extract image info
    # create df
    image_df = pd.DataFrame({
        "directory": [dirname(x) for x in image],
        "filename": [basename(x) for x in image]})
    # images_df = images_df.iloc[:10]
    # get metadata
    image_meta = image_df["filename"].str.extract(image_extractor + f'{image_suffix}')
    if "channel" in image_meta.columns:
        image_meta["channel"] = image_meta["channel"].astype("str") + image_meta["self_generated"]
    else:
        image_meta["channel"] = '0'
    
    image_meta.pop('self_generated')
    # print(image_meta)
    
    # merge with directory, filename
    image_df = pd.concat([image_df, image_meta], axis=1)

    # set index
    # images_df = images_df.reset_index()
    image_index_cols = ['directory'] + [x for x in image_extract_cols if x not in ['channel', 'self_generated']]
    image_df = image_df.set_index(image_index_cols)
    # add prefix to channel
    image_df["channel"] = "ch" + image_df["channel"]
    # get unique channel names
    image_colnames = natsorted(image_df['channel'].unique().tolist())
    # reshape channel to wide
    image_df = image_df.pivot(columns="channel", values="filename")
    print(f"{len(image_df)} grouped intensity images: {list(image_df)}")

    # mask files (cellpose masks) -------------------------------------
    if mask_suffix is None:
        mask_colnames = None
        df = image_df
    else:
        mask = glob(join(data_dir, image_subdir, '**/*') + mask_suffix, recursive=True)
        if len(mask) == 0:
            mask_colnames = None
            df = image_df
        else:
            mask_df = pd.DataFrame({
                "directory": [dirname(x) for x in mask],
                "filename": [basename(x) for x in mask]})
            # masks_df = masks_df.iloc[:10]
            # get metadata
            mask_meta = mask_df["filename"].str.extract(mask_extractor + f'{mask_suffix}')
            
            if "channel" not in mask_meta.columns:
                mask_meta["channel"] = '0'
            # print(mask_meta)
            
            # merge
            mask_df = pd.concat([mask_df, mask_meta], axis=1)
            # set index
            # masks_df = masks_df.reset_index()
            mask_df = mask_df.set_index(image_index_cols)
            # get unique mask names
            mask_colnames = mask_df['mask_name'].unique().tolist()
            # mask_df = mask_df.drop(columns=['channel'])
            # reshape channel to wide
            mask_df = mask_df.pivot(columns="mask_name", values="filename")
            print(f"{len(mask_df)} grouped object masks: {list(mask_df)}")

            # merge image and mask
            # print(image_df)
            # print(mask_df)
            df = pd.merge(image_df, mask_df, how='left', left_index=True, right_index=True)
            # print(df)

    # keep NA rows
    if remove_na_row:
        df = df.dropna()
    # reset index
    df = df.reset_index()
    print(f"{len(df)} groups after merging images and masks")

    # # add well data
    # row_abc = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    # # row number to abc: '02' to 'B', column keep A02, A03 cellprofiler like stype
    # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(str)
    # # row number to abc: '02' to 'B', column to A2, A3 harmony like stype
    # # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(int).astype(str)
    # df.insert(loc = 1, column = 'well', value = well_value)

    # # change to absolute path
    # df['directory'] = [os.path.abspath(f) for f in df['directory']]
    # # to unix style
    # df['directory'] = ['/'.join(x.split('\\')) for x in df['directory']]
    # # to windows style
    # # df['directory'] = df['directory'].replace('\\\\', '\\\\\\\\', regex=True)

    # metadata colnames
    metadata_colnames = image_index_cols

    # order images
    if 'prefix' in df.columns.to_list():
        df['prefix'] = df['prefix'].astype("str")
    if 'position' in df.columns.to_list():
        df['position'] = df['position'].astype("int64")
    if 'timepoint' in df.columns.to_list():
        df['timepoint'] = df['timepoint'].astype("int64")
    # df = df.sort_values(by=['prefix', 'position', 'timepoint'])

    # add position metadata
    if position_metadata is not None:
        # add position_metadata for nikon nd2 exported project
        if exists(join(data_dir, position_metadata)):
            pos_meta = pd.read_csv(join(data_dir, position_metadata))
            # remove NA rows which were interupted
            pos_meta = pos_meta.dropna(subset=[pos_index_col])
            # remove duplicate
            pos_meta = pos_meta.drop_duplicates(
                subset=['prefix','Time [s]', pos_index_col], keep='first')
            # print(pos_meta[[pos_index_col,"Position Name"]])

            # check if first row Position Name is not None, then do following steps
            if not pd.isna(pos_meta.iloc[0]["Position Name"]):
                print("add position well matadata")
                # split well and field
                pos_meta_well = pos_meta["Position Name"].str.extract(
                    "(?P<row>[A-Z])(?P<column>[0-9]+)#(?P<field>.+)")
                # pos_meta_well["column"] = pos_meta_well["column"].str.pad(2, "left", "0")
                pos_meta_well["well"] = pos_meta_well["row"] + pos_meta_well["column"]
                if "T Index" in pos_meta.columns:
                    # in case of partial NaN
                    pos_meta["T Index"] = pos_meta["T Index"].fillna(1)
                    pos_meta = pd.concat([pos_meta["prefix"],
                                          pos_meta[pos_index_col].rename("position"),
                                          pos_meta["T Index"].rename("timepoint"),
                                          pos_meta_well], axis=1)
                    pos_meta = pos_meta.drop(["row", "column"], axis=1)
                    pos_meta["prefix"] = pos_meta["prefix"].astype("str")
                    # pos_meta['timepoint'] = pos_meta['timepoint'].astype("int64")
                    # pos_meta['position'] = pos_meta['position'].astype("int64")
                    # merge to df
                    df = pd.merge(pos_meta, df, how="right", on=["prefix", "position", "timepoint"])
                    
                else:
                    pos_meta = pd.concat([pos_meta["prefix"],
                                          pos_meta[pos_index_col].rename("position"),
                                          pos_meta_well], axis=1)
                    pos_meta = pos_meta.drop(["row", "column"], axis=1)
                    pos_meta["prefix"] = pos_meta["prefix"].astype("str")
                    # pos_meta['position'] = pos_meta['position'].astype("int64")
                    # merge to df
                    df = pd.merge(pos_meta, df, how="right", on=["prefix", "position"])

                # add metadata cols
                metadata_colnames.extend(["well", "field"])

    # format df to cellprofiler dataloader
    # refer to: https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-4.2.4/modules/fileprocessing.html
    if cellprofiler_style:
        # image
        for ch in image_colnames:
            df[f'Image_PathName_{ch}'] = df['directory']
            df = df.rename(columns={f'{ch}': f'Image_FileName_{ch}'})
        # mask
        if mask_colnames is not None:  # if masks existed or not
            for mask_colname in mask_colnames:
                df[f'Image_ObjectsPathName_mask_{mask_colname}'] = df['directory']
                df = df.rename(
                    columns={f'{mask_colname}': f'Image_ObjectsFileName_mask_{mask_colname}'})
        # metadata
        for meta in metadata_colnames:
            df = df.rename(columns={f'{meta}': f'Metadata_{meta}'})
            
    return {'df': df,
            'metadata_colnames': metadata_colnames,
            'intensity_colnames': image_colnames,
            'mask_colnames': mask_colnames}



def read_img_from_dataset(
    dataset,
    row_index: int = 0,
    channel_names: list[list[str]] = [['ch1','ch2'], ['ch3']],
    channel_weights: list[list[float]] = [[1,2], [1]],
    channel_merge: str = "mean",
    resize_factor: float = 1
    ) -> np.ndarray:
    """
    Read image from dataset
    Args:
        dataset: result from images_to_dataset
        row_index: index of the dataset['df'] row
        channel_names: list of channel names, 1st layer list means channels to use, 2nd layer list means merged channel (method by channel_merge)
        channel_weights: list of channel weights, 1st layer list means channels to use, 2nd layer list means merged channel (method by channel_merge)
        channel_merge: str, method to merge channels
        resize_factor: float, resize the image by multipy factor
    Returns:
        np.ndarray, 2D image array [HW]
    """
    
    df = dataset['df']
    assert len(channel_names) <= 3, 'channel_names must be less/equal than 3'
    list_flat = sum(channel_names, [])
    assert all([x in df.columns.to_list() for x in list_flat]), 'all channel_names must be in df columns'
    assert len(channel_names) == len(channel_weights), 'different length of channel names and weights'
    
    img_mutli_chan = []
    for idx, channs in enumerate(channel_names):
        if len(channs) > 1:
            img_paths = [join(df['directory'][row_index], x) 
                         for x in df[channs].iloc[row_index].values.flatten().tolist()]
            if channel_merge == "mean":
                norm_weights = [i / sum(channel_weights[idx]) for i in channel_weights[idx]]
                imgs = [imread(j) * norm_weights[i] for i, j in enumerate(img_paths)]
                img = np.mean(imgs, axis=0).astype('uint16')
            elif channel_merge == "sum":
                imgs = [imread(j) * channel_weights[idx][i] for i, j in enumerate(img_paths)]
                img = np.sum(imgs, axis=0).astype('uint16')
            else:
                print("channel_merge must be one of (mean, sum)")
        else:
            img_path = join(df['directory'][row_index], 
                            df[channs].iloc[row_index].values[0])
            img = imread(img_path)
        
        # merge to list
        img_mutli_chan.append(img)
    
    # list to array
    img_mutli_chan = np.stack(img_mutli_chan)
    # resize
    if resize_factor != 1:
        img_mutli_chan = rescale(img_mutli_chan, scale=resize_factor, channel_axis=0, order=1)

    return img_mutli_chan



def cellpose_segment_dataset(
        dataset,
        channel_names: list[list[str]] = [['ch1','ch2'], ['ch3']],
        channel_weights: list[list[float]] = [[1,2], [1]],
        channel_merge: str = "mean",
        model_name: str = "cpsam",
        diameter: int = 30, # 0 or None mean auto size
        normalize = {'percentile': [0.5, 99.5]},
        resize_factor: float = 1,
        mask_name: str = 'cell',
        reduce_mask_name: bool = True,
        overwrite_mask: bool = False,
        flow_threshold: float = 0.4,
        cellprob_threshold: float = 0,
        device: str = None,
        gpu_batch_size: int = 64,
        ) -> None:
    """
    Segment dataset using cellpose model
    Args:
        dataset: result from images_to_dataset
        channel_names: list of channel names
        channel_weights: list of channel weights
        channel_merge: str, method to merge channels
        model_name: str, model name, any models stored in cellpose models folder '~/cellpose/models'
        diameter: int, diameter
        normalize: True
        resize_factor: float, resize the image
        mask_name: str, mask name
        reduce_mask_name: bool, reduce mask name to ch0 from any channel inputs
        overwrite_mask: bool, overwrite mask
        flow_threshold: float, flow threshold
        cellprob_threshold: float, cell probability threshold
    """
    df = dataset['df']

    # initialize model
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    elif device == "gpu":
        device = torch.device('cuda')
    elif device == "cpu":
        device = torch.device('cpu')

    model = models.CellposeModel(device=device, pretrained_model=model_name)
    diameter = None if diameter is None or diameter == 0 else int(diameter * resize_factor)
    
    print(f'segmentation model: {model_name}')
    print(f'segmentation channel: {channel_names}')
    print(f'segmentation diameter: {diameter}')
    print(f'segmentation mask_name: {mask_name}')
    
    clean_up_memory()
    
    for idx in tqdm(range(len(df))):
        fname = join(df['directory'].iloc[idx], df.iloc[idx][channel_names[0][0]])
        save_stem = join(dirname(fname), f'{splitext(basename(fname))[0]}') + '_cp_masks'
        # rename chN to ch1 to reduce any reductancy
        # if need more mask, give different {mask_name}
        if reduce_mask_name:
            save_stem = re.sub('_ch\\d+_', '_ch0_', save_stem)
        save_fname = f'{save_stem}_{mask_name}.png'
        if not overwrite_mask and exists(save_fname):
            continue

        try:
            img = read_img_from_dataset(
                dataset, idx, channel_names, channel_weights, channel_merge, resize_factor)
            
            mask, flow, style, *_ = model.eval(
                img, 
                batch_size = gpu_batch_size,
                channel_axis = 0,
                normalize = normalize,
                diameter = diameter,
                flow_threshold = flow_threshold,
                cellprob_threshold = cellprob_threshold)

            # # check mask
            # from cellpose import plot
            # from matplotlib.pyplot import figure
            # plot.show_segmentation(figure(figsize=(7,7)), img, mask, flow)

            # filter if no masks found (in case of cellprofiler no object error)
            if len(np.unique(mask)) > 1:
                if resize_factor != 1:
                    mask = rescale(mask, 1/resize_factor, order=0)
                cp_io.save_masks(
                    img, closing(mask), flow, file_names=save_stem, suffix=f'_{mask_name}')

        except Exception as e: 
            print(e)
            clean_up_memory()


    return None


def clean_up_memory():
    if torch.cuda.is_available():
        print('clean up GPU memory')
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
    gc.collect()
    return None


def ilastik_sinlge(images,
                   ilastik_proj,  # load saved checkpoint
                   # operetta exported measurement directory
                   ilastik_exec="D:/Software/Ilastik/ilastik-1.4.0b21-gpu/ilastik.exe"
                   ):
    # device = torch.device("cuda:" + str(gpu_id) if torch.cuda.is_available() else 'cpu')
    # print('using device:', device)
    # if torch.cuda.is_available():
    #     print('device name:', torch.cuda.get_device_name(gpu_id))
    #     torch.cuda.empty_cache()

    # cmd = f'\""{ilastik_exec}\"" --headless --project=\""{ilastik_proj}\"" {images}'
    # cmd = f'nohup \""{ilastik_exec}\"" --headless --project=\""{ilastik_proj}\"" {images} &'
    # cmd = f'LAZYFLOW_THREADS=2 LAZYFLOW_TOTAL_RAM_MB=2000 \""{ilastik_exec}\"" --readonly --headless --project=\""{ilastik_proj}\"" {images}'
    cmd = f'\""{ilastik_exec}\"" --headless --project=\""{ilastik_proj}\"" {images}'
    # subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process = Popen(cmd, shell=True, stdin=PIPE, stdout=DEVNULL, stderr=DEVNULL, text=True)
    out, err = process.communicate()
    return None


def ilastik_multiple(
    measurement: str = None,
    ilastik_proj: str = None,
    subset_pattern: str = None,
    image_subdir: str = "images",
    image_suffix: str = ".tiff",
    output_suffix: str = "Probabilities",
    ilastik_exec: str = "/home/hao/ilastik/run_ilastik.sh",
    thread_num: int = 16,
    overwrite_mask: bool = False,
    test: bool = False
):

    # enter image dir
    cur_dir = os.getcwd()
    os.chdir(os.path.join(measurement, image_subdir))
    # get images
    files = glob("*" + image_suffix)
    # remove any files containing output suffix
    files = [f for f in files if not re.search(output_suffix + ".tiff", f)]
    print("found %d intensity files" % (len(files)))

    # subset_pattern = "r0[1]c0[1-2]f.+p.+-"
    # subset_pattern = None
    if subset_pattern is not None:
        files = [f for f in files if re.search(subset_pattern, f)]
        # print(files[0])
        print("found %d files, filtered by %s" % (len(files), subset_pattern))

    # remove founded segmentated images
    if overwrite_mask:
        print(f"overwrite any existed masks, will process {len(files)} images")
    else:
        files_save_name = [f.replace(".tiff",  "_" + output_suffix + ".tiff") for f in files]
        # print(files_save_name[0:5])
        segmentated_files = []
        for idx, f in enumerate(files_save_name):
            if os.path.isfile(f):
                segmentated_files.append(files[idx])
        print("ignore %d predicted masks from total %d images contain %s" %
              (len(segmentated_files), len(files), subset_pattern))
        files = list(set(files) - set(segmentated_files))
    files = natsorted(files)
    # print(files[0:5])

    if len(files) == 0:
        print("no valid files, check if already done or inputs")
        return None

    # test mode
    if test is True:
        print("enable test mode, process random 4 files")
        files = random.sample(files, 4)

    # batch in each cpu core
    batch_images = []
    batch_size = math.ceil(len(files) / thread_num)
    print(f'using {thread_num} cores with each batch size of {batch_size}')
    for i in range(0, len(files), batch_size):
        batch = files[i:i+batch_size]
        # ' '.join([str(elem) for elem in batch])
        batch_images.append(' '.join(f'\"{w}\""' for w in batch))
    # print(batch_images)

    # run in parallel
    # partial_func = partial(ilastik_sinlge,
    #                        ilastik_exec=ilastik_exec,
    #                        ilastik_proj=ilastik_proj)

    start_time = time.time()
    jobs = []
    for batch in batch_images:
        process = Process(target=ilastik_sinlge, args=(batch, ilastik_proj, ilastik_exec,))
        jobs.append(process)

    [j.start() for j in jobs]
    [j.join() for j in jobs]
    [j.close() for j in jobs]
    
    os.chdir(cur_dir)

    # print("elapsed time: " + str(timedelta(seconds=int(time.time() - start_time))))
    print("done")






def softmax(x):
    """
    Computes the softmax function for a given input array.

    Args:
        x (numpy.ndarray): The input array (logits).

    Returns:
        numpy.ndarray: The softmax output (probabilities).
    """
    # Subtract the maximum value for numerical stability
    e_x = np.exp(x - np.max(x, axis=-1, keepdims=True))
    return e_x / np.sum(e_x, axis=-1, keepdims=True)



def inference_sliding_window(
        model, 
        image: np.ndarray, 
        window_size: int = 512, 
        overlap: float = 0.0, 
        batch_size: int = 16, 
        device: str = 'cuda'
        ) -> np.ndarray:
    """
    Perform sliding window inference on large images using a trained segmentation model.
    
    Args:
        model: Trained PyTorch segmentation model
        image: Input np.ndarray image of shape (H,W)
        window_size: Size of sliding window (assumes square window)
        overlap: Overlap ratio between windows (0-1)
        batch_size: Batch size for inference
        device: Device to run inference on ('cuda' or 'cpu')
    Returns:
        np.ndarray: segmentation mask, same shape with input image shape (H,W)
    """
    if image.ndim == 2:
        image = image.unsqueeze(0).unsqueeze(0)  # Add batch, channel dimension

    _, _, height, width = image.shape
    stride = int(window_size * (1 - overlap))
    
    # Pad image if needed
    pad_h = (stride - height % stride) % stride
    pad_w = (stride - width % stride) % stride
    image = F.pad(image, (0, pad_w, 0, pad_h), mode='reflect')
    
    windows = []
    coords = []
    
    # Extract windows
    for y in range(0, height + pad_h - window_size + 1, stride):
        for x in range(0, width + pad_w - window_size + 1, stride):
            window = image[:, :, y:y+window_size, x:x+window_size]
            windows.append(window)
            coords.append((y, x))
    # print(len(windows))

    # Run inference in batches
    predictions = []
    with torch.no_grad():
        for i in range(0, len(windows), batch_size):
            batch = torch.cat(windows[i:i+batch_size]).to(device)
            pred = torch.sigmoid(model(batch))
            predictions.extend(pred.cpu())
    
    # Initialize output mask
    output = torch.zeros((1, 1, height + pad_h, width + pad_w))
    count = torch.zeros((1, 1, height + pad_h, width + pad_w))
    
    # Reconstruct full image
    for pred, (y, x) in zip(predictions, coords):
        pred = pred.unsqueeze(0)  # Add batch dimension back
        output[:, :, y:y+window_size, x:x+window_size] += pred
        count[:, :, y:y+window_size, x:x+window_size] += 1
    
    # Average overlapping regions
    output = output / count
    
    # Remove padding
    output = output[:, :, :height, :width]
    
    return output


def spot_segment_dataset(
        dataset, # result from images_to_dataset
        channel_names: list[str] = ['ch1'],
        pred_names: list[str] = ["spot_pred.png"],
        overlap: float = 0.0,
        threshold: float = 0.5,
        model_path: str = "/media/hao/Data1/spot_segmentation_traindata/logs_v3/model.pth",
        model_input_size: int = 512,
        model_batch_size: int = 16,
        overwrite_spot: bool = False
        ) -> None:
    """
    Segment dataset using spot segmentation model
    Args:
        dataset: result from images_to_dataset
        channel_names: list of channel names
        pred_names: list of prediction name
        overlap: float, overlap ratio between windows
        threshold: float, threshold for segmentation
        model_path: str, path to the model, model should have a dict with key 'model'
        model_input_size: int, input size of the model
        model_batch_size: int, batch size for the model
        overwrite_spot: bool, overwrite spot
    Returns:
        None
        it will write prediction to dataset['df']['directory']
    """
    df = dataset['df']
    for ch in channel_names:
        assert ch in dataset['intensity_colnames'], 'channel_names must be in intensity_colnames'  
    
    # load model
    ckpt = torch.load(model_path, weights_only=False)
    model = ckpt['model']
    model.eval()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    
    try:
        for idx, data in tqdm(df.iterrows(), total=len(df)):
            for ch_idx, chan in enumerate(channel_names):
                save_name = join(
                    data['directory'], 
                    f'{splitext(basename(data[chan]))[0]}_{pred_names[ch_idx]}')
                # print(save_name)
                if not overwrite_spot and exists(save_name):
                    continue

                img = iio.imread(join(data['directory'], data[chan]))
                img = img / 65535
                img = torch.from_numpy(img).unsqueeze(0).unsqueeze(0)
                img = img.to(device, dtype=torch.float)
                
                with torch.no_grad():
                    outputs = inference_sliding_window(
                        model, img, model_input_size, overlap, model_batch_size, device)
                    outputs = outputs.cpu().numpy().squeeze(0).squeeze(0)
                    if threshold is not None:
                        outputs[outputs > threshold] = 1
                        outputs = outputs.astype(np.uint8)
                iio.imwrite(save_name, outputs)
                
    except Exception as e:
        print(e)
        


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
    if exists(save_path):
        print(f'{save_path} exist, exit!')
        return None
    elif save_path:
        conn = sqlite3.connect(save_path)
    else:
        print(f'save measurement: {save_path}')
    
    props_dict = {}
    # measure each mask
    for mask_name in mask_colnames:
        try:
            props = []
            for idx, data in tqdm(df.iterrows(), total=len(df)):
            
                img = [imread(join(data['directory'], data[x])) for x in intensity_colnames]
                img = np.stack(img, axis=2) # skimage regionprops accept imgs as HWC
                mask = iio.imread(join(data['directory'], data[mask_name]))
                prop = measure_image(mask, img, **kwargs)

                # rename channel
                # some columns need to be renamed first in case of naming error
                prop.rename(columns={'area_filled':'area', 
                                      'equivalent_diameter_area':'diameter'}, inplace=True)
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



class SegmentationDatasetFromArrayList(Dataset):
    """Dataset class for semantic segmentation that takes lists of numpy arrays as input.
    
    Args:
        img_list (list): List of images as numpy arrays
        mask_list (list): List of masks as numpy arrays
        ids (list): List of ids
        image_size (int): Size to resize images to. Defaults to 256.
    """
    def __init__(self, 
                 img_list: list, 
                 mask_list: list = None, 
                 ids: list = None, 
                 target_size: int = 256) -> None:
        self.img_list = img_list
        self.mask_list = mask_list if mask_list else [None] * len(img_list)
        self.ids = ids if ids else [None] * len(img_list)
        self.target_size = target_size if target_size else 256
        self.transform = A.Compose([
            A.Resize(self.target_size, self.target_size),
            # A.PadIfNeeded(min_height=target_size, min_width=target_size, 
            #               border_mode=cv2.BORDER_CONSTANT, value=0),
            ToTensorV2()])
    
    def __len__(self) -> int:
        return len(self.img_list)

    def __getitem__(self, idx: int) -> tuple:
        # normalize to [0,1]
        img = (self.img_list[idx] / 65535).astype(np.float32)
        mask = self.mask_list[idx]
        if mask is None:
            mask = np.zeros_like(img)

        if self.transform:
            augmented = self.transform(image=img, mask=mask)
            img, mask = augmented["image"], augmented["mask"]

        return img, mask, self.ids[idx]


def measure_dataset_by_model(
        res, # resulf from images_to_dataset
        model, # should be a model class
        model_eval_colname: str = None, # run model channel
        intensity_colnames: list[str] = None, # measurnent channel, if None, use res['intensity_colnames']
        mask_colname: str = None, # cell mask name, if None, use res['mask_colnames'][0]
        img_size: int = 256, # image size resize to
        save_path: str = None, # file path, if None, not save to db
        table_name: str = 'spot', # save db table name
        ) -> pd.DataFrame:
    
    df = res['df']
    mask_colname = mask_colname if mask_colname else res['mask_colnames'][0]
    intensity_colnames = intensity_colnames if intensity_colnames else res['intensity_colnames'] 
    intensity_colnames = list(set(intensity_colnames + [model_eval_colname]))
    
    # config db
    # fname = join(dirname(df['directory'].iloc[0]), save_db_name)
    if save_path:
        if exists(save_path):
            conn = sqlite3.connect(save_path)
            tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", conn)
            print(tables)
            if table_name in tables:
                conn.close()
                print(f'{table_name} exist in {save_path}, check it!')
                return None
        else:
            conn = sqlite3.connect(save_path)

    try:
        props_dataset = []
        
        # prepare model
        if isinstance(model, str):
            model = torch.load(model)
            if isinstance(model, dict):
                model = model['model']
        model.eval()
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model.to(device)
        
        props_dataset = []
        for idx, data in tqdm(df.iterrows(), total=len(df)):
            mask_whole_image = iio.imread(join(data['directory'], data[mask_colname]))
            
            # props_multi_channel = []
            intensity_imgs = [
                iio.imread(join(data['directory'], data[x])) 
                    for x in intensity_colnames]
            intensity_imgs = np.stack(intensity_imgs, axis=2)
                        
            # prepate single cell list
            cell_ids, cropped_cells = crop_cell(
                mask_whole_image, intensity_imgs, target_size=img_size, 
                pad_square=True, rotate=None, 
                exclude_border=False, max_cell_per_image=None)
            # print(len(cropped_cells))
            
            # extract channel for model eval from cropped cells
            ch_use_index = intensity_colnames.index(model_eval_colname)
            cropped_cells_single_channel = [
                x[:,:,ch_use_index] for x in cropped_cells]
            # print(cropped_cells_single_channel[0].shape)
            # print(cropped_cells_single_channel[0].dtype)
            
            # prepare dataset
            dataset = SegmentationDatasetFromArrayList(
                img_list=cropped_cells_single_channel, 
                mask_list=None, ids=cell_ids, target_size=img_size)
            # print(dataset[0][0].shape)
            dataloader = DataLoader(dataset, batch_size=64, shuffle=False,
                                    drop_last=False)
            # run model
            masks_pred = []
            ids = []
            with torch.no_grad():
                for img, mask, label in dataloader:
                    img = img.to(device)
                    pred = model(img)
                    pred = torch.sigmoid(pred).cpu().numpy()
                    # pred = (pred > 0.5).astype(np.uint8)
                    masks_pred.append(pred)
                    ids.append(label.cpu().numpy())
            masks_pred = np.concatenate(masks_pred, axis=0)
            ids = list(np.concatenate(ids, axis=0))
            # print(masks_pred.shape)
            # plt.imshow(masks_pred[0,0,:,:])
            
            # reshape to BHWC for skimage measure
            masks_pred = set_img_channel(masks_pred, order='HWC')
            # probability to mask, remove mask channel
            masks_pred = (masks_pred > 0.5)[:,:,:,0].astype(np.uint8)
            
            # measure intensity channels          
            props = []
            for b_idx in range(masks_pred.shape[0]): # enumerate(masks_pred_list):
                # print(np.unique(masks_pred[b_idx,:,:]))
                prop = regionprops_table(
                    masks_pred[b_idx,:,:], cropped_cells[b_idx],
                    properties=['label','area_filled','intensity_mean'])
                if len(prop) > 0:
                    prop = pd.DataFrame(prop)
                    prop['label'] = ids[b_idx]
                    props.append(prop)
            props = pd.concat(props)
            
            # rename intensity channel
            # some columns need to be renamed first in case of naming error
            props.rename(columns={'area_filled':'area', 
                                  'equivalent_diameter_area':'diameter'}, inplace=True)
            props.rename(columns=lambda x: re.sub('bbox-','bbox_', x), inplace=True)
            props.rename(columns=lambda x: re.sub('eigvals-','eigvals_', x), inplace=True)
            for ch_idx, ch in enumerate(intensity_colnames):
                props.rename(columns=lambda x: x.replace(f'-{ch_idx}', f'-{ch}'), inplace=True)
            
            # insert metadata
            meta = data.to_dict()
            [meta.pop(col, None) for col in res['intensity_colnames'] + res['mask_colnames']]
            props = props.assign(**meta)


            # save to db
            if save_path:
                props.to_sql(table_name, conn, if_exists='append', index=False)

            # to whole dataset
            props_dataset.append(props)

        # to long format
        props_dataset = pd.concat(props_dataset)

    except Exception as e: 
        print(e)
    finally: 
        pass

    # close db
    if save_path is not None: 
        conn.close()

    return props_dataset



def crop_cell(
        mask: np.ndarray, # HW
        img: np.ndarray, # HWC, HW, output will preserve original ndim
        target_size: int = None,
        pad_square: bool = True,
        rotate: chr = None,
        exclude_border: bool = True, 
        max_cell_per_image: int = None,
        ) -> list[list[int], list[np.ndarray]]:
    """
    Crop cell
    Args:
        mask: np.ndarray, shape HW
        img: np.ndarray, shape HWC, HW, output will preserve original ndim
        target_size: int, target size
        pad_square: bool, pad square
        rotate: chr, rotate direction
        exclude_border: bool, exclude border
        max_cell_per_image: int, max cell per image
    Returns:
        list[list[int], list[np.ndarray]], cell ids and cropped cells
    """
    # check input
    if img.ndim == 2:
        raw_ndim = 2 # record original ndim
        img = np.expand_dims(img, 2) # HWC
    assert img.shape[:2] == mask.shape, "mask and img have different shape"
    
    # Clear border cell
    if exclude_border:
        mask = clear_border(mask)
    
    cell_ids = np.unique(mask)
    cell_ids = cell_ids[cell_ids != 0].ravel().tolist() # exclude background 0

    # Sample cell
    if max_cell_per_image and len(cell_ids) > max_cell_per_image:
        cell_ids = sample(cell_ids, max_cell_per_image)

    cropped_cells =[]
    for cell_id in cell_ids:
        cell_mask = (mask == cell_id).astype(np.uint8)
        
        # Get bounding box of the cell
        coords = np.argwhere(cell_mask)
        y_min, x_min = coords.min(axis=0)
        y_max, x_max = coords.max(axis=0)
        
        # Crop
        cropped_mask = cell_mask[y_min:y_max+1, x_min:x_max+1]
        cropped_cell = img[y_min:y_max+1, x_min:x_max+1, :] * cropped_mask[:,:,None]
        # imshow(cropped_cell[0,:,:])
        
        # Rotate mask and cell
        if rotate:
            cropped_mask, cropped_cell = rotate_cell(
                cropped_mask, cropped_cell, direction=rotate)
            # update coords
            coords = np.argwhere(cropped_mask)
            y_min, x_min = coords.min(axis=0)
            y_max, x_max = coords.max(axis=0)
        
        if pad_square:
            # Calculate size of the square
            height = y_max - y_min + 1
            width = x_max - x_min + 1
            square_size = max(height, width)
            # Calculate padding
            pad_y = (square_size - height) // 2
            pad_x = (square_size - width) // 2
            
            # Pad to square
            cropped_cell = np.pad(
                cropped_cell, ((pad_y, square_size - height - pad_y), (pad_x, square_size - width - pad_x), (0, 0)),
                mode='constant', constant_values=0)
            # cropped_mask = np.pad(
            #     cropped_mask, ((pad_y, square_size - height - pad_y), (pad_x, square_size - width - pad_x)),
            #     mode='constant', constant_values=0)
        
        # Resize to target size
        if target_size:
            cropped_cell = resize(
                cropped_cell, (target_size, target_size, cropped_cell.shape[2]), 
                preserve_range=True, anti_aliasing=True)
            # # Smooth resized mask edge
            # smoothed_mask = gaussian(cropped_padded_resized[0,:,:], sigma=1)
            # imshow(smoothed_mask)
            # smoothed_mask = (smoothed_mask > 0.5).astype(int)
            # imshow(smoothed_mask)
            # cropped_padded_resized[0,:,:] = smoothed_mask
            # imshow(cropped_padded_resized[0, :, :])
        
        # restore ndim
        if 'raw_ndim' in locals():
            cropped_cell = cropped_cell[:,:,0]

        cropped_cells.append(cropped_cell)

    return [cell_ids, cropped_cells]



def rotate_cell(
        cropped_mask: np.ndarray,
        cropped_cell: np.ndarray,
        direction: str = "horizonal"
        ) -> list[np.ndarray, np.ndarray]:
    """
    Rotate cell
    Args:
        cropped_mask: np.ndarray, shape HW
        cropped_cell: np.ndarray, shape HWC, HW
        direction: str, direction to rotate
    Returns:
        list[np.ndarray, np.ndarray], rotated mask and cell
    """
    # check input
    assert cropped_cell.shape[:2] == cropped_mask.shape, "mask and cell have different shape"

    # Get object properties
    props = regionprops(cropped_mask)
    
    # For one object, get its orientation
    orientation = props[0].orientation  # Orientation in radians

    # Calculate padding to prevent cropping during rotation
    max_dim = max(cropped_cell.shape)
    pad_width = max_dim
    if cropped_cell.ndim == 3:
        padded_img = np.pad(cropped_cell, ((pad_width, pad_width), (pad_width, pad_width), (0,0)), mode='constant', constant_values=0)
    else:
        padded_img = np.pad(cropped_cell, ((pad_width, pad_width), (pad_width, pad_width)), mode='constant', constant_values=0)
    padded_mask = np.pad(cropped_mask, ((pad_width, pad_width), (pad_width, pad_width)), mode='constant', constant_values=0)

    # Rotate image to align long axis horizontally
    angle_degrees = -np.degrees(orientation)
    if direction == "horizonal":
        angle_degrees = angle_degrees + 90
    rotated_img = ndi.rotate(padded_img, angle_degrees, reshape=False, mode='nearest', axes=(0, 1))
    rotated_mask = ndi.rotate(padded_mask, angle_degrees, reshape=False, mode='constant', cval=0)

    # Crop back to original image dimensions
    coords = np.argwhere(rotated_mask)
    y_min, x_min = coords.min(axis=0)
    y_max, x_max = coords.max(axis=0)

    if cropped_cell.ndim == 3:
        rotated_img = rotated_img[y_min:y_max+1, x_min:x_max+1, :]
    else:
        rotated_img = rotated_img[y_min:y_max+1, x_min:x_max+1]
    rotated_mask = rotated_mask[y_min:y_max+1, x_min:x_max+1]

    return [rotated_mask, rotated_img]




def crop_cell_from_dir(
        data_dir: str,
        destination_dir: str = None,
        subset_well: list[str] = None,
        subet_position: list[int] = None,
        mask_column_name: str = "cp_masks_cell",
        intensity_channel_name: list[str] = None,
        target_size: int = None,
        rotation: str = "horizonal",
        exclude_border: bool = True,
        max_image_per_well: int = 2,
        max_cell_per_image: int = 100,
        save_multi_channel_tiff: bool = True
        ) -> None:
    """
    Crop cell from directory
    Args:
        data_dir: str, path to the data directory
        destination_dir: str, path to save the cropped cells
        subset_well: list[str], list of well names to subset
        subet_position: list[int], list of position names to subset
        mask_column_name: str, mask column name
        intensity_channel_name: list[str], list of intensity channel names
        target_size: int, target size
        rotation: str, rotation direction
        exclude_border: bool, exclude border
        max_image_per_well: int, max image per well
        max_cell_per_image: int, max cell per image
    """
    # measurement_dir = "/media/hao/Data/Project_barcode/2024-12-06_BC_test"
    # destination_dir = "/media/hao/Data1/single_cell_dataset/293T_barcode"
    if destination_dir is None:
        destination_dir = os.path.join(data_dir, "cropped_cell")
    
    # get dataset
    res = images_to_dataset(
        data_dir,
        subset_pattern=None,
        image_subdir="images",
        image_suffix='.tiff',
        mask_suffix='.png',  # None for no use
        remove_na_row=True,
        cellprofiler_style=False,
        position_metadata="position.csv")
    
    if intensity_channel_name is None:
        intensity_channel_name = res["intensity_colnames"]
    df = res["df"][res["metadata_colnames"] + intensity_channel_name + [mask_column_name]]
    # subset image
    if subset_well is not None:
        df = df[df['well'].isin(subset_well)]
    if subet_position is not None:
        df = df[df['position'].isin(subet_position)]
    # sample image
    if max_image_per_well is not None and "well" in df.columns.tolist():
        df = df.groupby(["prefix","well"]).sample(n=max_image_per_well, replace=False)

    print(f'cropping images: {intensity_channel_name}')
    print(f'cropping mask: {mask_column_name}')
    print(f'cropping size: {target_size}')
    # get cell
    for row_idx in tqdm(range(len(df))):
        # create dir for cells from same well/position
        if "well" in res["df"].columns.tolist():
            save_dir = os.path.join(
                destination_dir, 
                df["prefix"].values[row_idx] + "__" +  df["well"].values[row_idx])
        else:
            save_dir = os.path.join(
                destination_dir, 
                df["prefix"].values[row_idx] + "__p" +  str(df["position"].values[row_idx]))
        os.makedirs(save_dir, exist_ok=True)
        
        # read mask img
        img = [iio.imread(join(data_dir, "images", df[ch].values[row_idx])) 
               for ch in intensity_channel_name]
        img = np.stack(img, 2) # HWC
        mask = iio.imread(join(data_dir, "images", df[mask_column_name].values[row_idx]))

        # crop cell
        cell_ids, cropped_cells = crop_cell(
            mask, img, target_size=target_size, pad_square=True, rotate=rotation, 
            exclude_border=exclude_border, max_cell_per_image=max_cell_per_image)
  
        # preserve cell dim == 3
        if cropped_cells.ndim == 2:
            cropped_cells = [x[:,:, None] for x in cropped_cells]
        
        # save cell
        for id, cell in zip(cell_ids, cropped_cells):
            if save_multi_channel_tiff:
                f_name = os.path.join(
                    save_dir, f"p{df['position'].values[row_idx]}_c{id}.tiff")
                iio.imwrite(f_name, cell)
            else:
                for ch_idx, ch in enumerate(intensity_channel_name):
                    f_name = os.path.join(
                        save_dir, f"p{df['position'].values[row_idx]}_c{id}_{ch}.tiff")
                    iio.imwrite(f_name, cell[:,:,ch_idx])

# crop_cell_from_dir('/media/hao/Data/Project_barcode/2025-01-06_BC_test', 
#                    '/media/hao/Data1/single_cell_dataset/single_channel_128/2025-01-06_BC_test',
#                    intensity_channel_name=['ch1','ch4'],
#                    target_size=128)





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
        image = resize(image, (target_size , target_size), anti_aliasing=True, preserve_range=True)
        mask = resize(mask, (target_size, target_size), order=0, anti_aliasing=False, preserve_range=True)

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
    # cv2.distanceTransform don't accept bool type
    if mask.dtype == bool:
        mask = mask.astype(np.uint8)
    dist_img = cv2.distanceTransform(mask, distanceType=cv2.DIST_L2, maskSize=5).astype(np.float32)
    center_y, center_x = np.where(dist_img==dist_img.max())
    center_y, center_x = center_y.mean(), center_x.mean() # there are sometimes cases where there are multiple values returned for the visual center
    # print(center_y, center_x)
    
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



def generate_test_image_and_mask(width=200, height=200, mask_shape="star"):
    """
    Generates a test image with a radial gradient and a mask with a specified shape.

    Parameters:
    - width: Image width.
    - height: Image height.
    - mask_shape: "star", "ellipse", "round", or "rectangle".
    """
    
    import numpy as np
    import skimage.draw
    import matplotlib.pyplot as plt
    
    # Generate radial gradient image
    center_y, center_x = height // 2, width // 2
    y, x = np.ogrid[:height, :width]
    distances = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    max_distance = np.max(distances)
    gradient = (distances / max_distance)
    image = (gradient * 255).astype(np.uint8)

    # Generate mask based on mask_shape
    mask = np.zeros((height, width), dtype=np.uint8)
    outer_radius = min(width, height) // 3

    if mask_shape == "star":
        # Five-pointed star
        angles = np.linspace(np.pi / 2, 2 * np.pi + np.pi / 2, 6)[:-1]
        x_points = center_x + outer_radius * np.cos(angles)
        y_points = center_y - outer_radius * np.sin(angles)
        rr, cc = skimage.draw.polygon(y_points, x_points, shape=mask.shape)
        mask[rr, cc] = 1

    elif mask_shape == "ellipse":
        # Ellipse
        rr, cc = skimage.draw.ellipse(center_y, center_x, outer_radius, outer_radius * 0.7, shape=mask.shape)
        mask[rr, cc] = 1

    elif mask_shape == "round":
        # Round
        rr, cc = skimage.draw.disk((center_y, center_x), outer_radius, shape=mask.shape)
        mask[rr, cc] = 1

    elif mask_shape == "rectangle":
        # Rectangle
        rr, cc = skimage.draw.rectangle(
            (center_y - outer_radius // 2, center_x - outer_radius),
            extent=(outer_radius, outer_radius * 2),
            shape=mask.shape,
        )
        mask[rr, cc] = 1

    # mask = mask.astype(bool)
    
    # Display the image and mask
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(image, cmap='gray')
    plt.title("Radial Gradient Image")
    plt.subplot(1, 2, 2)
    plt.imshow(mask, cmap='gray')
    plt.title("Mask Shape: Ellipse")
    plt.show()

    return image, mask
# image, mask = generate_test_image_and_mask(mask_shape='rectangle')
# radial_distribution(mask, image, num_bins=5, display=True)




#%% test
if False:
    import os
    import sys
    sys.path.append('/home/hao/Documents/GitHub/myScript/python_functions')
    import torch
    
    
    indi_dir = '/media/hao/Data/Others/2025-03-02_Calp-50x_Chan-I/images'
    img = iio.imread(indi_dir+'/60x large 31-138+31-137+Ca2+ ionomycin 1uM__t93_p1_z1_ch1.tiff')
    img2 = np.repeat(img[:,:,None], 2, axis=2)
    mask = iio.imread(indi_dir+'/60x large 31-138+31-137+Ca2+ ionomycin 1uM__t93_p1_z1_ch2_cp_masks_cell.png')
    ids, cells = crop_cell(mask, img2, target_size=128)
    cells[0].shape
    x0 = measure_image(mask, img2)
    x0 = measure_image(mask, img2, extra_properties=[
        granularity, radial_distribution,])
    x0.columns


    data_dir = "/media/hao/Data/Others/2025-03-02_Calp-50x_Chan-I"
    dataset = images_to_dataset(data_dir, subset_pattern=None, remove_na_row=True)
    dataset['df'] = dataset['df'].iloc[0:2]
    
    model_path = "/media/hao/Data1/single_cell_dataset/cell_spot_cnn_training_data/logs/v3/model.pth"
    model = torch.load(model_path, weights_only=False)
    print(model)
    # model = model['val_transform']
    trained_model = model['model']
    
    x1 = measure_dataset(dataset, save_path=None, extra_properties=None)
    x2 = measure_dataset(dataset, save_path=None, 
                         extra_properties=[granularity, radial_distribution,]
                         )
    x3 = measure_dataset_by_model(
        dataset, trained_model, 
        model_eval_colname="ch1",
        intensity_colnames=["ch1","ch2"],
        mask_colname='cp_masks_cell',
        img_size=256,
        save_path=None, 
        table_name='spot')
    x3.columns