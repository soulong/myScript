# -*- coding: utf-8 -*-
"""
search images under directory, then assemble to df used for cellprofiler dataloader

Created on Wed Nov  8 17:03:05 2023
By Hao He
"""

import argparse
import os
# import glob, re
# from natsort import natsorted
# import itertools
# import pandas as pd

from functions import image_mask_to_df_nd2



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'CellPose wrapper on images exported from Operetta')
    parser.add_argument('--dir', type = str, nargs = "+", default = '.', metavar = ".",
        help = 'directory list, containing images as subdirectory')
    parser.add_argument('--subset_pattern', type = str, default = None, metavar = 'None',
        help = 'use regrex to subset images, eg. "r0[1]c0[1-2]f.+p.+-" used to subset A1 and A2 wells')
    parser.add_argument('--image_suffix', type = str, default = '.tiff', metavar = '.tiff',
        help = 'image suffix')
    parser.add_argument('--mask_suffix', type = str, default = '.png', metavar = '.png',
        help = 'mask suffix')
    parser.add_argument('--keep_na_rows', action = "store_true", help = 'filter out rows with any NA')
    parser.add_argument('--cellprofiler_style', action = "store_true", help = 'output cellprofiler style')
    parser.add_argument('--position_metadata', type = str, default = 'position_metadata.csv', 
        metavar = 'position_metadata.csv',
        help = 'nd2 exported position medata file under measurement directory, if none no use position metadata')
    args = parser.parse_args()

    
    measurements = args.dir
    subset_pattern =args.subset_pattern
    image_suffix = args.image_suffix
    mask_suffix = args.mask_suffix
    keep_na_rows = args.keep_na_rows
    cellprofiler_style = args.cellprofiler_style
    position_metadata = args.position_metadata
    position_metadata = None if position_metadata == "None" else position_metadata

    if True:
        measurements=["D:\\Postdoc_Data\\Project_cmpd_screen\\2024-04-10_r-secretase_cmpd_validate"]
        # subset_pattern = "20x-plate"
        subset_pattern = None
        image_suffix = ".tiff"
        mask_suffix = ".png"
        keep_na_rows=False
        cellprofiler_style=True
        position_metadata = "position_metadata.csv"
        # position_metadata = None
    
    if measurements is None:
        quit('no valid measurement found, check it!')
    
    for measurement in measurements:
        print(f">>>> {measurement}")
        res = image_mask_to_df_nd2(measurement, 
                                   subset_pattern, 
                                   "images", 
                                   image_suffix, 
                                   mask_suffix, 
                                   keep_na_rows, 
                                   cellprofiler_style,
                                   position_metadata)
        print(f'metadata_colnames: {res["metadata_colnames"]}')
        res["df"].to_csv(os.path.join(measurement, "cp_dataloader.csv"), index=False)
        
    print("\n---------- Done ---------")
        


