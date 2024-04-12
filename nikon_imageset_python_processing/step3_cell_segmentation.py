# -*- coding: utf-8 -*-

import argparse
# import glob, re
import os, math, time
from datetime import timedelta
import cv2
import numpy as np
# import pandas as pd
# from natsort import natsorted
import torch
from cellpose import models, io
from tifffile import imread
from tqdm import tqdm
from skimage.morphology import closing

import logging
import sys
logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings('ignore')
# warnings.filterwarnings('ignore', '.*UserWarning: Applied workaround for CuDNN issue.*', )

from functions import image_mask_to_df_nd2


# #############################
# # for single channel, simplified usage
# if False:
#     model = models.CellposeModel(
#         device = torch.device("cuda:0"), 
#         pretrained_model = "cyto2")
#     resize = 0.5
#     diameter = 250
#     measurement = "D:/Postdoc_Data/2023-12-06_RTK_PPI_mApple_Ctrl_ZQ"
#     files = glob.glob(os.path.join(measurement, "images", "*.tiff" ))
#     files = [f for f in files if re.search("beas2b.*c002.tif", f)]
#     # files = files[0:5] 
    
#     for f in tqdm(files):
#         img = imread(f)
#         img = cv2.resize(img, (0, 0), fx = resize, fy = resize)
#         mask, flow, _ = model.eval(
#             img,
#             channels = [0, 0],
#             diameter = diameter * resize,  # None, auto size
#             resample = True,
#             net_avg = True,
#             flow_threshold = 0.4,
#             cellprob_threshold = 0,
#             min_size = int(3.14 * math.pow(diameter/2, 2) / 10)
#             )
#         if len(np.unique(mask)) > 1:
#             if resize != 1:
#                 mask = cv2.resize(mask, (0, 0), fx=(1/resize), fy=(1/resize), interpolation=cv2.INTER_NEAREST)
#             io.save_masks(img, closing(mask), flow,
#                 # file_names="test2.tiff".split(".tiff")[0],
#                 file_names = os.path.splitext(f)[0],
#                 suffix = '_' + "cell",
#                 png = True, save_txt = False, save_ncolor = False,
#                 save_flows = False, save_outlines = False)
# ##########################



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



def _prepare_model(
        model_name: str = 'cyto2',
        gpu_id: int = 0
        ):
    # set model
    device = torch.device("cuda:" + str(gpu_id) if torch.cuda.is_available() else 'cpu')
    print('using device:', device)
    # if torch.cuda.is_available():
    #     print('device name:', torch.cuda.get_device_name(gpu_id))
    if model_name in ["cyto", "cyto2", "nuclei"]:
        model = models.Cellpose(device = device, model_type = model_name)
    else:
        model = models.CellposeModel(device = device, pretrained_model = model_name)

    return model


def _prepare_img(
        df,
        row_index: int = 0,
        channel_names: list = ['ch1'],
        channel_weights: list = [1],
        resize: int = 0.5
        ):
    channels_error = [chan for chan in channel_names if chan not in df.columns.values.tolist()]
    if len(channels_error) > 0:
        quit(f'channels used for segmentation not exist: {channels_error}')

    # calculate weights
    if len(channel_names) > 1:
        channel_weights = [i / sum(channel_weights) for i in channel_weights]
        # print("compute mask using channels: %s" % channel_names)
        # print("compute mask using weights: %s" % channel_weights)
        # read img
        image_dir = df['directory'][row_index] + '/'
        image_filenames = df.iloc[[row_index]][channel_names].values.flatten().tolist()
        images = [image_dir + x for x in image_filenames]
        imgs = []
        for idx, item in enumerate(images):
            img = imread(item)
            imgs.append(img * channel_weights[idx])
        img = sum(imgs).astype('uint16')
        # print(np.max(img), np.min(img))

    else:
        image_dir = df['directory'][row_index] + '/'
        image_filename = df.iloc[[row_index]][channel_names].values.flatten().item()
        img = imread(image_dir + image_filename)

    if resize != 1:
        img = cv2.resize(img, (0, 0), fx = resize, fy = resize)

    return img


def _run_cellpose(
        model,
        img,
        diameter: int = 100,
        flow_threshold: float = 0.4,
        cellprob_threshold: float = 0
        ):
    res = model.eval(
        img,
        channels = [0, 0],
        diameter = diameter,  # None, auto size
        resample = True,
        net_avg = True,
        flow_threshold = flow_threshold,
        cellprob_threshold = cellprob_threshold,
        min_size = int(3.14 * math.pow(diameter/2, 2) / 10)
        )

    return [res[0], res[1]]


def _segment_df(
        df,
        channel_names: list = ['ch1'],
        channel_weights: list = [1],
        model_name: str = "cyto2",
        diameter: int = 80,
        resize: float = 1,
        output_suffix: str = 'cell',
        overwrite_mask: bool = False,
        flow_threshold: float = 0.4,
        cellprob_threshold: float = 0
        ):
    # initialize model
    model = _prepare_model(model_name)

    for idx in tqdm(range(len(df))):
        fname = df[['directory']].iloc[idx].values.flatten().item() + '/' + \
            df.iloc[[idx]][channel_names[0]].values.item()

        # skip image with existed corresponding mask
        if not overwrite_mask:
            final_fname = os.path.splitext(fname)[0] + "_cp_masks_" + output_suffix + ".png"
            if os.path.exists(final_fname): continue

        try:
            img = _prepare_img(df, idx, channel_names, channel_weights, resize)
            mask, flow = _run_cellpose(model, img, int(diameter * resize), flow_threshold, cellprob_threshold)
            # filter if no masks found (in case of cellprofiler no object error)
            if len(np.unique(mask)) > 1:
                if resize != 1:
                    mask = cv2.resize(mask, 
                                      (0, 0), fx=(1/resize), fy=(1/resize), 
                                      interpolation=cv2.INTER_NEAREST)
                io.save_masks(img, closing(mask), flow,
                    # file_names="test2.tiff".split(".tiff")[0],
                    file_names = os.path.splitext(fname)[0],
                    suffix = '_' + output_suffix,
                    png = True, save_txt = False, save_ncolor = False,
                    save_flows = False, save_outlines = False)
        except:
           print(f'error skip: {os.path.basename(fname)}')
           pass

    return None




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'CellPose wrapper on images exported from Operetta')
    parser.add_argument('--dir', type = str, nargs = "+", default = '.', metavar = ".",
        help = 'directory list, containing images as subdirectory')
    parser.add_argument('--image_suffix', type = str, default = '.tiff', metavar = '.tiff',
        help = 'image suffix')
    parser.add_argument('--subset_pattern', type = str, default = None, metavar = 'None',
        help = 'use regrex to subset images, eg. "r0[1]c0[1-2]f.+p.+-" used to subset A1 and A2 wells')
    parser.add_argument('--channel_names', type = str, nargs = "+", default = ['ch1'], metavar = "ch1",
        help = 'compute mask use indicated channels, "ch1" means use ch1 images, "ch1 ch2" means merge ch1 and ch2, '
                'channel_names must be one or more of image_set columns')
    parser.add_argument('--channel_weights', type = float, nargs = "+", default = [1], metavar = 1.0,
        help = 'indicate each channel weight, same length as --channel_names')
    parser.add_argument('--model', type = str, default = 'cyto2', metavar = 'cyto2',
        help = 'analysis type, default: cyto2')
    parser.add_argument('--diameter', type = int, default = 50, metavar = 50,
        help = 'mean diameter of raw image object size, typically ')
    # parser.add_argument('--flow_threshold', type = float, default = 0.4, metavar = 0.4,
    #     help = 'cellpose flow threshold, value from 0 to 1')
    # parser.add_argument('--cellpro_threshold', type = float, default = 0, metavar = 0,
    #     help = 'cellpose mask threshold, value from -6 to 6')
    parser.add_argument('--resize', type = float, default = 0.5, metavar = 0.5,
        help = 'scale factor of images before run cellpose, 1 mean no scaling')
    parser.add_argument('--output_suffix', type = str, default = 'cell', metavar = 'cell',
        help = 'output mask file suffix: "_cp_masks_{output_suffix}.png"')
    parser.add_argument('--overwrite_mask', action = "store_true", help = 'overwrite existed masks')
    parser.add_argument('--test', action = 'store_true', help = 'test on random indicated number of images')
    # parser.add_argument('--delay', type=int, default=0, metavar=0, help='delay run, default 0 seconds')
    args = parser.parse_args()

    
    measurements = args.dir
    image_suffix = args.image_suffix
    subset_pattern =args.subset_pattern
    channel_names = args.channel_names
    channel_weights = args.channel_weights
    model = args.model
    diameter = args.diameter
    resize = args.resize
    output_suffix = args.output_suffix
    overwrite_mask = args.overwrite_mask
    test = args.test
    
    if True:
        measurements=["D:\\Postdoc_Data\\2024-02-19_src_project"]
        image_suffix = '.tiff'
        # subset_pattern = ".*SURF-V3[+]RBD-CPSURF.*_p50_.*"
        subset_pattern = None
        channel_names = ["ch1"]
        channel_weights = [1]
        # model = "cyto2"
        # model = "D:/Postdoc_Data/Cellpose_refine_dataset/cyto2_with_puncta_myc/models/cyto2_nuclei_puncta_myc"
        model = "C:\\Users\\haohe\\.cellpose\\models\\cyto2_with_emlalk"
        diameter = 60
        resize = 1
        output_suffix = "cell"
        overwrite_mask = False
        test = False


    if measurements is None:
        quit('no valid measurement found, check it!')
            
    for measurement in measurements:
        print("\n>>>>> processing: " + measurement)
        start_time = time.time()

        # measurement = "/media/hao/Data1/2023-09-01_293T_Panel_Cmpds_TC__2023-09-01T10_08_42-Measurement 1"
        image_set = image_mask_to_df_nd2(measurement, 
                                         subset_pattern, 
                                         "images",
                                         image_suffix=image_suffix,
                                         mask_suffix=None,
                                         keep_na_rows=True,
                                         cellprofiler_style=False,
                                         position_metadata=None)
        # # write images.csv file
        # image_set['df'].to_csv(os.path.join(measurement, 'images.csv'))

        # check input channel_names
        channels_error = [x for x in channel_names if x not in image_set['intensity_colnames']]
        if len(channels_error) > 0:
            quit(f'{channels_error} are not in {image_set["intensity_colnames"]}, check it!')
        
        if test:
            _segment_df(
                image_set['df'].sample(n=4, ignore_index=True), 
                channel_names, 
                channel_weights, 
                model, 
                diameter,
                resize,
                output_suffix,
                overwrite_mask)
        else:
            _segment_df(
                image_set['df'], 
                channel_names, 
                channel_weights, 
                model, 
                diameter,
                resize, 
                output_suffix,
                overwrite_mask)

        print("elapsed time: " + str(timedelta(seconds = int(time.time() - start_time))))

    print("----------------- done -----------------")
    

