# -*- coding: utf-8 -*-
"""
A script to convert Nikon .nd2 file to indivual tiff files
Seperate position and channel to suffix

Create at 2023-11-08
By Hao He
"""

import argparse
import nd2
import os
from sys import exit
from glob import glob
from tifffile import imwrite
import cv2
import numpy as np
import pandas as pd
from tqdm import tqdm
from time import sleep



def nd2_to_images(measurement,
                 in_dir="nd2",
                 out_dir="images",
                 metadata_file="position_metadata.csv",
                 resize=None, # (2048, 2048) # resize whole original image
                 crop_N=None # 4, # split big tiled image to NxN small separated images
                 ):
    
    if not os.path.exists(measurement):
        exit(f"{measurement} not found, check it")
        
    dir_in = os.path.join(measurement, in_dir)
    dir_out = os.path.join(measurement, out_dir)
    os.makedirs(dir_out, exist_ok=True)
    
    nd2_files = glob(dir_in + "/**/*.nd2", recursive=True)
    if len(nd2_files) == 0:
        exit("no nd2 files found, check it")
        
    # f_name = "V1-MCHERRY-MIFP-EBFP2-100-2.nd2"
    # f_path = os.path.join(dir_in, f_name)
    
    # rename nd2 files to remove " ", in case of some unintended error
    print("\ncheck and correct file name")
    for idx, f_path in enumerate(nd2_files):
        if " " in f_path:
            base_dir = os.path.dirname(f_path)
            base_name = os.path.basename(f_path)
            base_name_new = base_name.replace(" ", "-")
            f_path_new = os.path.join(base_dir, base_name_new)
            print("replace blank to '-'")
            print(f"from > {base_name}")
            print(f"to >>> {base_name_new}")
            os.rename(f_path, f_path_new)
            nd2_files[idx] = f_path_new
    
    
    # extract position metadata
    print("\nextract and merge position metadata")
    position_metadata_f = os.path.join(measurement, metadata_file)
    if os.path.exists(position_metadata_f): 
        exit(f"{position_metadata_f} found, please delete it first")
    meta_list = []
    for f_path in nd2_files:
        # get metadata
        with nd2.ND2File(f_path) as ndfile:
            meta = pd.DataFrame(ndfile.events())
            # to avoid index probelm
            meta["P Index"] = meta["P Index"] + 1
        # print(meta)
        
        # add filename prefix
        fname = os.path.splitext(os.path.basename(f_path))[0]
        # print(fname)
        meta.insert(0, "prefix", fname)
        meta_list.append(meta)
        
    # merge to one file
    meta_df = pd.concat(meta_list, axis=0)
    meta_df.to_csv(position_metadata_f, index=False, mode="w")

    
    # lopp process
    for f_path in nd2_files:
        
        fname = os.path.splitext(os.path.basename(f_path))[0]
        print(f"\n>> {fname}")
        
        # > get metadata
        meta = nd2.ND2File(f_path)
        print(meta.sizes)
        # meta.ndim
        # meta.dtype
        
        # > z-stack, using max projection
        max_projection = True if "Z" in list(meta.sizes.keys()) else False
        # max_projection = False
        
        # > channel
        if "C" in list(meta.sizes.keys()):
            channel_names = [meta.metadata.channels[i].channel.name for i in range(meta.sizes["C"])]
            print(channel_names) 
        
        
        # read data
        dat = meta.asarray()
        # dat.shape
        
        # full recover array to [T, P, Z, C, Y, X]
        if "C" not in list(meta.sizes.keys()):
            dat = np.expand_dims(dat, -3)
        if "Z" not in list(meta.sizes.keys()):
            dat = np.expand_dims(dat, -4)
        if "P" not in list(meta.sizes.keys()):
            dat = np.expand_dims(dat, -5)
        if "T" not in list(meta.sizes.keys()):
            dat = np.expand_dims(dat, -6)
    
        
        # save to indivual tiff file
        for ti, t in enumerate(dat): 
            # print("t shape:", ti, "___", t.shape)
            
            for pi, p in enumerate(tqdm(t, desc=f"T{ti+1}: >>P")):
                # print("p shape:", p.shape)
                
                if max_projection:
                    # print("mip")
                    p = np.max(p, axis=0)
                    p = np.expand_dims(p, 0)
                    # print("f shape mip:", f.shape)
                    
                for zi, z in enumerate(p):
                    # print("z shape:", z.shape)
                    
                    for ci, c in enumerate(z):
                        # print("t:",ti+1, " f:",fi+1, " z:",zi+1, " c:",ci+1)
                        
                        # resize
                        if resize is not None:
                            c = cv2.resize(c, dsize=resize, interpolation=cv2.INTER_CUBIC)
                        
                        # split into pieces
                        if crop_N is None:
                            # save
                            imwrite(
                                os.path.join(dir_out, f"{fname}__t{ti+1}_p{pi+1}_z{zi+1}_c{ci+1}.tiff"), 
                                c, compression="zlib")
                            sleep(0.001)
                            pass
                        else:
                            M = c.shape[0]//crop_N
                            N = c.shape[1]//crop_N
                            tiles = [c[x:x+M,y:y+N] for x in range(0,c.shape[0],M) for y in range(0,c.shape[1],N)]
                            
                            for xi, x in enumerate(tiles):
                                imwrite(
                                    os.path.join(dir_out, f"{fname}__t{ti+1}_p{pi+1}000{xi}_z{zi+1}_c{ci+1}.tiff"), 
                                    x, compression="zlib")
                                sleep(0.001)
                                pass
        
        # close file
        meta.close()
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Convert nikon nd2 to tiff files')
    parser.add_argument('--dir', type = str, nargs = "+", default = '.', metavar = ".",
        help = 'directory list, containing images as subdirectory')
    parser.add_argument('--in_dir', type = str, default = 'nd2', metavar = 'nd2',
        help = 'where nd2 files locate')
    parser.add_argument('--out_dir', type = str, default = 'images', metavar = 'images',
        help = 'where tiff files saved')
    parser.add_argument('--metadata_file', type = str, default = 'position_metadata.csv', 
                        metavar = 'position_metadata.csv',
        help = 'nd2 exported position medata file under measurement directory')
    parser.add_argument('--resize', type = int, nargs = "+", default = None, metavar = 'None',
        help = 'resize output images, two int list')
    parser.add_argument('--crop_N', type = int, default = None, metavar = 'None',
        help = 'crop (resized) image to NxN small separated ones')
    args = parser.parse_args()

    
    measurements = args.dir
    in_dir = args.in_dir
    out_dir =args.out_dir
    metadata_file = args.metadata_file
    resize = args.resize
    crop_N = args.crop_N
    
    if True:
        measurements=["D:\\Postdoc_Data\\Project_cmpd_screen\\2024-04-10_r-secretase_cmpd_validate"]
        in_dir="."
        out_dir="images"
        metadata_file="position_metadata.csv"
        # resize=[8340,8340]
        resize=None
        # crop_N = 4
        crop_N = None



    for measurement in measurements:
        print("\n>>>>> processing: " + measurement)
        nd2_to_images(measurement, in_dir, out_dir, metadata_file, 
                      resize=resize, crop_N=crop_N)
    
    print("\n---------- Done ---------")

            
            
            
            
            
            
            
            
            
            
            