# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 13:29:41 2023

@author: haohe
"""


import argparse
import os, re, random, time
import warnings
warnings.filterwarnings('ignore')

from datetime import timedelta
from shutil import copy
# import glob
# import pandas as pd
import numpy as np
from tifffile import imread, imwrite
# from functools import partial
# from tqdm.contrib.concurrent import process_map
from tqdm import tqdm
# from natsort import natsorted
import pickle

import jax
jax.config.update('jax_platform_name', 'cpu')
from basicpy import BaSiC

from functions import image_mask_to_df_nd2


# # for nd2 exported format
# def image_to_df_nd2(
#         measurement = ".",
#         subset_pattern="",
#         subdir='images',
#         image_suffix='.tiff',
#         # mask_suffix='.png',
#         keep_na_rows=True
#         ):
#     # measurement = "D:\\Others\\2023-11-13_NMYC_spark_value"
#     ## get all images
#     if subdir is None:
#         files_all = glob.glob(os.path.join(measurement, "*" ))
#     else:
#         files_all = glob.glob(os.path.join(measurement, subdir, "*"))
        
#     # subset images
#     # remove Probabilities files
#     files_all = [f for f in files_all if not re.search("_Probabilities", f)]
    
#     # print(subset_pattern)
#     if subset_pattern is not None:
#         files_all = [f for f in files_all if re.search(subset_pattern, f)]
        
#     # sort files
#     files_all = natsorted(files_all)
#     # check if any masks exists
#     if len(files_all) == 0:
#         print("no valid images/masks found, check directory or subset_pattern")
#         return None

#     ## intensity files
#     images = [f for f in files_all if re.search(f'{image_suffix}', f)]
#     if len(images) == 0:
#         print("no valid images found, please check it")
#         return None
#     # create df
#     image_df = pd.DataFrame({  # "path": files,
#         "directory": [os.path.dirname(x) for x in images],
#         "filename": [os.path.basename(x) for x in images]})
#     # images_df = images_df.iloc[:10]
#     # get metadata
#     extractor = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_c(?P<channel>[0-9]{1,})(?P<self_generated>.*)' + f'{image_suffix}'
#     image_meta = image_df["filename"].str.extract(extractor)
#     image_meta["channel"] = image_meta["channel"].astype("str") + image_meta["self_generated"]
#     image_meta.pop('self_generated')
#     # merge with directory, filename
#     image_df = pd.concat([image_df, image_meta], axis=1)
    
#     # set index
#     # images_df = images_df.reset_index()
#     image_df = image_df.set_index(["directory", "prefix", "timepoint", "position", "stack"])
#     # add prefix to channel
#     image_df["channel"] = "ch" + image_df["channel"]
#     # get unique channel names
#     image_colnames = image_df['channel'].unique().tolist()
#     # reshape channel to wide
#     image_df = image_df.pivot(columns="channel", values="filename")
#     print(f"they're total {len(image_df)} grouped intensity images: {list(image_df)}")

#     # keep NA rows
#     if not keep_na_rows:
#         image_df = image_df.dropna()
#     # reset index
#     image_df = image_df.reset_index()
#     print(f"they're total {len(image_df)} groups after merging images and masks")

#     # # add well data
#     # row_abc = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
#     # # row number to abc: '02' to 'B', column keep A02, A03 cellprofiler like stype
#     # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(str)
#     # # row number to abc: '02' to 'B', column to A2, A3 harmony like stype
#     # # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(int).astype(str)
#     # df.insert(loc = 1, column = 'well', value = well_value)

#     # change to absolute path
#     image_df['directory'] = [os.path.abspath(f) for f in image_df['directory']]
#     # to unix style
#     image_df['directory'] = ['/'.join(x.split('\\')) for x in image_df['directory']]
#     # to windows style
#     # df['directory'] = df['directory'].replace('\\\\', '\\\\\\\\', regex=True)

#     # metadata colnames
#     metadata_colnames = ["directory", "prefix", "timepoint", "position", "stack"]
    
#     return {'df':image_df, 
#             'metadata_colnames':metadata_colnames, 
#             'intensity_colnames':image_colnames}


# basic =BaSiC(fitting_mode='approximate', get_darkfield=True, mu_coef=12.5, working_size=256)
# print(basic.Config)
def _basic_fit(
        images,
        n_image = 50,
        enable_dfd = True,
        resize = 128
        ):
    images = images if isinstance(images, list) else [images]
    if len(images) > n_image:
        images = random.sample(images, k = n_image)
    # read images
    imgs = np.array([imread(f) for f in images])
    # fit basic model
    basic = BaSiC(get_darkfield = enable_dfd,
                  smoothness_flatfield = 1,
                  smoothness_darkfield = 1,
                  working_size = resize,
                  max_workers = 16)
    basic.fit(imgs)

    return basic


def _basic_fit_df(
        df,
        process_channels,
        target_dir = None, # like: xx1/xx2/xx3_Measurement 1, if None, basic_dir_name will be located under original Measurement directory
        n_image = 100,
        enable_dfd = True,
        resize = 128
        ):
    # df = pd.read_csv('/media/hao/Data1/HCI/image_analysis_scripts/test__2023-Measurement 1/images.csv')
    # process by intensity channels
    all_column_names = df.columns.values.tolist()
    intensity_colnames = [i for i in all_column_names if re.search('^ch[0-9]+$', i)]
    channels = intensity_colnames if process_channels is None else process_channels
    channels = [chan for chan in channels if chan in intensity_colnames]
    if len(channels) == 0:
        quit('provided channels are not found, check it')

    # set basic save path
    if target_dir is not None:
        basic_save_path = os.path.join(target_dir, os.path.basename(os.path.dirname(df['directory'].to_list()[0])), "BaSiC_model")
    else:
        basic_save_path = os.path.join(os.path.dirname(df['directory'].to_list()[0]), "BaSiC_model")
    os.makedirs(basic_save_path, exist_ok = True)

    for chan in channels:
        fname = f'{basic_save_path}/model_{chan}.pkl'
        print(f'fit BaSiC model on channel: {chan}')
        images = [os.path.join(df['directory'].to_list()[0], i) for i in df[chan].to_list()]
        model = _basic_fit(images, n_image, enable_dfd, resize)

        # save profile
        imwrite(f'{basic_save_path}/model_{chan}_dfd.tiff', model.darkfield, compression = 'zlib')
        imwrite(f'{basic_save_path}/model_{chan}_ffd.tiff', model.flatfield, compression = 'zlib')
        with open(fname, 'wb') as f:
            pickle.dump(model, f)

        # apply on random images
        test_images = random.sample(images, k = 2)
        test_images_transformed = _basic_transform(test_images, model, target_dir, timelapse=False, test_run=True)
        for i in range(test_images_transformed.shape[0]):
            copy(test_images[i], f'{basic_save_path}/{os.path.basename(test_images[i])}')
            imwrite(f'{basic_save_path}/{os.path.basename(test_images[i])}_cor.tiff',
                    test_images_transformed[i, :, :], compression = 'zlib')
    return None


def _basic_transform(
        images,
        fitted_model,
        target_dir = None, # like: xx1/xx2/xx3_Measurement 1
        timelapse = False,
        test_run = False
        ):
    # read images
    images = images if isinstance(images, list) else [images]
    imgs = [imread(f) for f in images]
    imgs = np.array(imgs)  # print(imgs.shape)
    # transform
    imgs_transformed = fitted_model.transform(imgs, timelapse = timelapse)
    # in case of over-range transform error
    if imgs.dtype == 'uint16':
        imgs_transformed[imgs_transformed > 65535] = 65535
        imgs_transformed[imgs_transformed < 0] = 0
    elif imgs.dtype == 'uint8':
        imgs_transformed[imgs_transformed > 255] = 255
        imgs_transformed[imgs_transformed < 0] = 0
    else:
        quit('image dtype error, not uint16 or uint8')
    imgs_transformed = imgs_transformed.astype(imgs.dtype)

    # overwrite raw images with corrected images
    if test_run:
        return imgs_transformed
    else:
        if target_dir is not None:
            for idx, f in enumerate(images):
                imwrite(os.path.join(target_dir, os.path.basename(os.path.dirname(os.path.dirname(f))),
                    "images", os.path.basename(f)),
                    imgs_transformed[idx, :, :], compression = 'zlib')
        else:
            for idx, f in enumerate(images):
                imwrite(f, imgs_transformed[idx, :, :], compression = 'zlib')
        return None


def _basic_transform_df(
        df,
        process_channels,
        target_dir = None # like: xx1/xx2/xx3_Measurement 1
        ):
    all_column_names = df.columns.values.tolist()
    intensity_colnames = [i for i in all_column_names if re.search('^ch[0-9]+$', i)]
    channels = intensity_colnames if process_channels is None else process_channels
    channels = [chan for chan in channels if chan in intensity_colnames]
    if len(channels) == 0:
        quit('provided channels are not found, check it')

    # use original Measurement dir is target_dir is None
    image_dirname = df['directory'].to_list()[0] # xx1/xx2/xx3_Measurement 1/Images
    if target_dir is None:
        model_dir = os.path.dirname(image_dirname)
    else:
        model_dir = os.path.join(target_dir, os.path.basename(os.path.dirname(image_dirname)))
        # mkdir for transformed images
        os.makedirs(os.path.join(model_dir, "images"), exist_ok = True)
    for chan in channels:
        fname = os.path.join(model_dir, "BaSiC_model", f"model_{chan}.pkl")
        print(f'transform BaSiC model on channel: {chan}')
        if not os.path.exists(fname):
            quit(f'BaSiC model not found: {fname}')
        else:
            images = [os.path.join(image_dirname, i) for i in df[chan].to_list()]
            with open(fname, 'rb') as f:
                model = pickle.load(f)
                
            # if test:
            #     images = random.sample(images, k=4)
            #     image_transformed = _basic_transform(images, model, overwrite_image=False)
            #     # save test images
            #     for i in range(image_transformed.shape[0]):
            #         copy(images[i], f'{save_dir}/{os.path.basename(images[i])}')
            #         imwrite(f'{save_dir}/{os.path.basename(images[i])} + "_cor.tiff"',
            #             image_transformed[i,:,:], compression='zlib')
            # else:
                
            for image in tqdm(images):
                _basic_transform(image, model, target_dir, timelapse=False, test_run=False)
                    
            # # parallel run
            # partial_func = partial(_basic_transform, fitted_model=model, 
            #                        target_dir=target_dir, timelapse=False)
            # _ = process_map(partial_func, images, max_workers = 1)

    return None




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'fit and transform shading correction with pyBASIC model')
    parser.add_argument('--dir', type = str, nargs = "+", default = '.', metavar = ".",
                        help = 'directories (parent or sub-folder) exported from operetta')
    parser.add_argument('--subset_pattern', type = str, default = None, metavar = 'None',
        help = 'use regrex to subset images, eg. "r0[1]c0[1-2]f.+p.+-" used to subset A1 and A2 wells')
    parser.add_argument('--image_suffix', type=str, default=".tiff", metavar=".tiff", help='image suffix')
    parser.add_argument('--process_channels', type = str, nargs = "+", default = "ch1", metavar = '["ch1",]',
                        help = 'process channels, if None, will use all intensity channels (chxx), DPC channel is usually not be processed')
    parser.add_argument('--target_dir', type = str, default = None, metavar = 'None', help = 'where the fit and transform data save, if not provided, use original dir')
    parser.add_argument('--n_image', type = int, default = 50, metavar = 50, help = 'fit model with n random images')
    parser.add_argument('--resize', type = int, default = 128, metavar = 128, help = 'resize image to smooth')
    parser.add_argument('--enable_dfd', action = 'store_true', help = 'enable darkfield correction')
    parser.add_argument('--run_type', type = str, default = 'f', metavar = '"f" or "t"',
                        help = 'which mode to run, default "f", f -> fit images, t -> apply transform')
    # parser.add_argument('--test', action = 'store_true', help = 'test transform on random indicated number of images')
    args = parser.parse_args()
    # print(args.process_channels)
    
    measurements = args.dir
    subset_pattern = args.subset_pattern
    image_suffix = args.image_suffix
    process_channels = args.process_channels
    target_dir = args.target_dir
    n_image = args.n_image
    resize = args.resize
    enable_dfd = args.enable_dfd
    run_type = args.run_type
    
    
    # image parameters
    # measurements = ["D:/Postdoc_Data/2024-01-06_condesate_reporter_LJY"]
    # subset_pattern = "100ng"
    # subset_pattern = None
    # subdir = 'images'
    # image_suffix = '.tiff'
    # keep_na_rows=False

    # fit parameters
    # process_channels = ["ch1","ch2","ch3","ch4"]
    # target_dir = None # or specific directory under the measurement dirctory", None will replace original images
    # n_image = 100
    # enable_dfd = False
    # resize = 128
    # run_type = "f"
    
    if measurements is None:
        quit('no valid measurement found, check it!')

    for measurement in measurements:
        print("\n>>>>> processing: " + measurement)
        start_time = time.time()
    
        # create basic_dir
        # if not os.path.exists(args.basic_dir):
        #    os.makedirs(args.basic_dir, exist_ok = True)
    
        # measurement = "/media/hao/Data/HCI/image_analysis_scripts/test__2023-Measurement 1"
        image_set = image_mask_to_df_nd2(measurement, 
                                         subset_pattern, 
                                         "images",
                                         image_suffix=image_suffix,
                                         mask_suffix=None,
                                         keep_na_rows=True,
                                         cellprofiler_style=False,
                                         position_metadata=None)
    
        # save fit and transform data to target_dir
        # if args.target_dir is not None:
        #     target_dir_measurement = os.path.join(args.target_dir, os.path.basename(measurement))
        # else:
        #     target_dir_measurement = measurement
        # print(target_dir_measurement)
        
        # run here
        # run_type = input("run type: \n[f] -> fit or [t] -> transform:\n>>")
        if run_type == "f":
            _basic_fit_df(image_set['df'], process_channels, target_dir, n_image, enable_dfd, resize)
        elif run_type == "t":
            _basic_transform_df(image_set['df'], process_channels, target_dir)
            
            # move images not transformed to target dir, if target is not None 
            if target_dir is not None:
                no_process_channels = set(image_set['intensity_colnames']).difference(set(process_channels))
                if len(no_process_channels) > 0:
                  image_dir = image_set['df']["directory"].values.tolist()[0]
                  image_targe_dir = os.path.join(target_dir, os.path.basename(measurement), "Images")
                  for copy_ch in no_process_channels:
                      print(f"copy {copy_ch} images to target dir")
                      copy_ch_images = [os.path.join(image_dir, fname) for fname in image_set['df'][copy_ch].values.tolist()]
                      _ = [copy(f, image_targe_dir) for f in tqdm(copy_ch_images)]
        else:
            quit(f'error mode: {run_type}')
    
        print("elapsed time: " + str(timedelta(seconds = int(time.time() - start_time))))
    
    print("\n----------------- done -----------------")
    
    
    
    
    
    
    
    
    
    
