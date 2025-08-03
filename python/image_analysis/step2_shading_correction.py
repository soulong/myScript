# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 13:29:41 2023

@author: haohe
"""


from functions import measurement_to_df_nd2
from basicpy import BaSiC
import jax
import pickle
from tqdm import tqdm
from tifffile import imread, imwrite
import numpy as np
from shutil import copy
from datetime import timedelta
import argparse
import os
import re
import random
import time
import warnings
warnings.filterwarnings('ignore')

# import glob
# import pandas as pd
# from functools import partial
# from tqdm.contrib.concurrent import process_map
# from natsort import natsorted

jax.config.update('jax_platform_name', 'cpu')


# basic =BaSiC(fitting_mode='approximate', get_darkfield=True, mu_coef=12.5, working_size=256)
# print(basic.Config)
def _basic_fit(
        images,
        n_image=50,
        enable_dfd=True,
        resize=128
):
    images = images if isinstance(images, list) else [images]
    if len(images) > n_image:
        images = random.sample(images, k=n_image)
    # read images
    imgs = np.array([imread(f) for f in images])
    # fit basic model
    basic = BaSiC(get_darkfield=enable_dfd,
                  smoothness_flatfield=1,
                  smoothness_darkfield=1,
                  working_size=resize,
                  max_workers=16)
    basic.fit(imgs)

    return basic


def basic_fit_df(
        df,
        process_channels,
        target_dir=None,  # like: xx1/xx2/xx3_Measurement 1, if None, basic_dir_name will be located under original Measurement directory
        n_image=50,
        enable_dfd=True,
        resize=128
):
    # df = pd.read_csv('/media/hao/Data1/HCI/image_analysis_scripts/test__2023-Measurement 1/images.csv')
    # process by intensity channels
    all_column_names = df.columns.values.tolist()
    intensity_colnames = [i for i in all_column_names if re.search('^ch[0-9]+$', i)]
    channels = intensity_colnames if process_channels is None else process_channels
    channels = [chan for chan in channels if chan in intensity_colnames]
    if len(channels) == 0:
        print('provided channels are not found, check it')
        return None

    # set basic save path
    if target_dir is not None:
        basic_save_path = os.path.join(target_dir, os.path.basename(
            os.path.dirname(df['directory'].to_list()[0])), "BaSiC_model")
    else:
        basic_save_path = os.path.join(os.path.dirname(df['directory'].to_list()[0]), "BaSiC_model")
    os.makedirs(basic_save_path, exist_ok=True)

    for chan in channels:
        fname = f'{basic_save_path}/model_{chan}.pkl'
        print(f'fit BaSiC model on channel: {chan}')
        images = [os.path.join(df['directory'].to_list()[0], i) for i in df[chan].to_list()]
        model = _basic_fit(images, n_image, enable_dfd, resize)

        # save profile
        imwrite(f'{basic_save_path}/model_{chan}_dfd.tiff', model.darkfield, compression='zlib')
        imwrite(f'{basic_save_path}/model_{chan}_ffd.tiff', model.flatfield, compression='zlib')
        with open(fname, 'wb') as f:
            pickle.dump(model, f)

        # apply on random images
        test_images = random.sample(images, k=2)
        test_images_transformed = _basic_transform(
            test_images, model, target_dir, timelapse=False, test_run=True)
        for i in range(test_images_transformed.shape[0]):
            copy(test_images[i], f'{basic_save_path}/{os.path.basename(test_images[i])}')
            imwrite(f'{basic_save_path}/{os.path.basename(test_images[i])}_cor.tiff',
                    test_images_transformed[i, :, :], compression='zlib')
    return None


def _basic_transform(
        images,
        fitted_model,
        target_dir=None,  # like: xx1/xx2/xx3_Measurement 1
        timelapse=False,
        test_run=False
):
    # read images
    images = images if isinstance(images, list) else [images]
    imgs = [imread(f) for f in images]
    imgs = np.array(imgs)  # print(imgs.shape)
    # transform
    imgs_transformed = fitted_model.transform(imgs, timelapse=timelapse)
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
                        imgs_transformed[idx, :, :], compression='zlib')
        else:
            for idx, f in enumerate(images):
                imwrite(f, imgs_transformed[idx, :, :], compression='zlib')
        return None


def basic_transform_df(
        df,
        process_channels,
        target_dir=None  # like: xx1/xx2/xx3_Measurement 1
):
    all_column_names = df.columns.values.tolist()
    intensity_colnames = [i for i in all_column_names if re.search('^ch[0-9]+$', i)]
    channels = intensity_colnames if process_channels is None else process_channels
    channels = [chan for chan in channels if chan in intensity_colnames]
    if len(channels) == 0:
        print('provided channels are not found, check it')
        return None

    # use original Measurement dir is target_dir is None
    image_dirname = df['directory'].to_list()[0]  # xx1/xx2/xx3_Measurement 1/Images
    if target_dir is None:
        model_dir = os.path.dirname(image_dirname)
    else:
        model_dir = os.path.join(target_dir, os.path.basename(os.path.dirname(image_dirname)))
        # mkdir for transformed images
        os.makedirs(os.path.join(model_dir, "images"), exist_ok=True)
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


def image_correction(
        image_set: dict = None,
        run_type: str = "f",  # t
        process_channels: list[str] = ["ch1"],
        n_image: int = 50,
        resize: int = 128,
        enable_dfd:  bool = True,
        target_dir: str = None
):
    if run_type == "f":
        basic_fit_df(image_set['df'], process_channels,
                     target_dir, n_image, enable_dfd, resize)

    elif run_type == "t":
        basic_transform_df(image_set['df'], process_channels, target_dir)

        # move images not transformed to target dir, if target is not None
        if target_dir is not None:
            no_process_channels = set(
                image_set['intensity_colnames']).difference(set(process_channels))
            if len(no_process_channels) > 0:
                image_dir = image_set['df']["directory"].values.tolist()[0]
                image_targe_dir = os.path.join(
                    target_dir, os.path.basename(measurement), "Images")
                for copy_ch in no_process_channels:
                    print(f"copy {copy_ch} images to target dir")
                    copy_ch_images = [os.path.join(image_dir, fname)
                                      for fname in image_set['df'][copy_ch].values.tolist()]
                    _ = [copy(f, image_targe_dir) for f in tqdm(copy_ch_images)]
    else:
        quit(f'error mode: {run_type}')

        return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='fit and transform shading correction with pyBASIC model')
    parser.add_argument('--measurements', type=str, default='.', metavar=".",
                        help='directories (parent or sub-folder) exported from operetta')
    parser.add_argument('--subset_pattern', type=str, default=None, metavar='None',
                        help='use regrex to subset images, eg. "r0[1]c0[1-2]f.+p.+-" used to subset A1 and A2 wells')
    parser.add_argument('--image_suffix', type=str, default=".tiff",
                        metavar=".tiff", help='image suffix')
    parser.add_argument('--process_channels', type=str, nargs="+", default="ch1", metavar='["ch1",]',
                        help='process channels, if None, will use all intensity channels (chxx), DPC channel is usually not be processed')
    parser.add_argument('--target_dir', type=str, default=None, metavar='None',
                        help='where the fit and transform data save, if not provided, use original dir')
    parser.add_argument('--n_image', type=int, default=50, metavar=50,
                        help='fit model with n random images')
    parser.add_argument('--resize', type=int, default=128,
                        metavar=128, help='resize image to smooth')
    parser.add_argument('--enable_dfd', action='store_true', help='enable darkfield correction')
    parser.add_argument('--run_type', type=str, default='f', metavar='"f" or "t"',
                        help='which mode to run, default "f", f -> fit images, t -> apply transform')
    # parser.add_argument('--test', action = 'store_true', help = 'test transform on random indicated number of images')
    args = parser.parse_args()
    # print(args.process_channels)

    measurements = args.measurements
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
        image_set = measurement_to_df_nd2(measurement,
                                          subset_pattern=subset_pattern,
                                          image_suffix=image_suffix,
                                          mask_suffix=None,
                                          keep_na_rows=True,
                                          cellprofiler_style=False,
                                          position_metadata=None)

        # run
        image_correction(image_set, run_type=run_type, process_channels=process_channels,
                         n_image=n_image, resize=resize, enable_dfd=enable_dfd, target_dir=target_dir)

        print("elapsed time: " + str(timedelta(seconds=int(time.time() - start_time))))

    print("\n----------------- done -----------------")
