# -*- coding: utf-8 -*-
"""
search images under directory, then assemble to df used for cellprofiler dataloader

Created on Wed Nov  8 17:03:05 2023
By Hao He
"""

import os, glob, re
from natsort import natsorted
import pandas as pd



def image_mask_to_df_nd2(
        measurement = ".",
        subset_pattern=None,
        subdir='images',
        image_suffix='.tiff',
        mask_suffix='.png', # None for no use
        keep_na_rows=True,
        cellprofiler_style=True,
        position_metadata="position_metadata.csv" # None for no usage
        ):
    
    # measurement = "D:\\Others\\2023-11-08_LJY_SPARK_spot_intensity"
    ## get all images
    if subdir is None:
        files_all = glob.glob(os.path.join(measurement, "*" ))
    else:
        files_all = glob.glob(os.path.join(measurement, subdir, "*"))
    # subset images
    # print(subset_pattern)
    if subset_pattern is not None:
        files_all = [f for f in files_all if re.search(subset_pattern, f)]
    # sort files
    files_all = natsorted(files_all)
    # check if any masks exists
    if len(files_all) == 0:
        print("no valid images/masks found, check directory or subset_pattern")
        return None

    ## intensity files
    images = [f for f in files_all if re.search(f'{image_suffix}', f)]
    if len(images) == 0:
        print("no valid images found, please check it")
        return None
    # create df
    image_df = pd.DataFrame({  # "path": files,
        "directory": [os.path.dirname(x) for x in images],
        "filename": [os.path.basename(x) for x in images]})
    # images_df = images_df.iloc[:10]
    # get metadata
    extractor = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_c(?P<channel>[0-9]{1,})(?P<self_generated>.*)' + f'{image_suffix}'
    image_meta = image_df["filename"].str.extract(extractor)
    image_meta["channel"] = image_meta["channel"].astype("str") + image_meta["self_generated"]
    image_meta.pop('self_generated')
    # merge with directory, filename
    image_df = pd.concat([image_df, image_meta], axis=1)
    
    # set index
    # images_df = images_df.reset_index()
    image_df = image_df.set_index(["directory", "prefix", "timepoint", "position", "stack"])
    # add prefix to channel
    image_df["channel"] = "ch" + image_df["channel"]
    # get unique channel names
    image_colnames = image_df['channel'].unique().tolist()
    # reshape channel to wide
    image_df = image_df.pivot(columns="channel", values="filename")
    print(f"they're total {len(image_df)} grouped intensity images: {list(image_df)}")

    ## masks (cellpose masks)
    if mask_suffix is None:
            # create df
        # return only images
        mask_colnames = None
        df = image_df
    else:
        mask = [f for f in files_all if re.search(f"{mask_suffix}", f)]
        if len(mask) == 0:
            mask_colnames = None
            df = image_df
        else:
            mask_df = pd.DataFrame({  # "path": files,
                "directory": [os.path.dirname(x) for x in mask],
                "filename": [os.path.basename(x) for x in mask]})
            # masks_df = masks_df.iloc[:10]
            # get metadata
            extractor = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_c(?P<channel>[0-9]{1,})_(?P<mask_name>.*)' + f'{mask_suffix}'
            mask_meta = mask_df["filename"].str.extract(extractor)
            # mask_meta['mask_name'] = 'mask_' + mask_meta["mask_name"] + '_ch' + mask_meta["channel"].astype("str")
    
            # merge
            mask_df = pd.concat([mask_df, mask_meta], axis=1)
            # set index
            # masks_df = masks_df.reset_index()
            mask_df = mask_df.set_index(["directory", "prefix", "timepoint", "position", "stack"])
            # get unique mask names
            mask_colnames = mask_df['mask_name'].unique().tolist()
            # mask_df = mask_df.drop(columns=['channel'])
            # reshape channel to wide
            mask_df = mask_df.pivot(columns="mask_name", values="filename")
            print(f"they're total {len(mask_df)} grouped object masks: {list(mask_df)}")
    
            # merge image and mask
            df = pd.merge(image_df, mask_df, left_index=True, right_index=True)

    # keep NA rows
    if not keep_na_rows:
        df = df.dropna()
    # reset index
    df = df.reset_index()
    print(f"they're total {len(df)} groups after merging images and masks")

    # # add well data
    # row_abc = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    # # row number to abc: '02' to 'B', column keep A02, A03 cellprofiler like stype
    # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(str)
    # # row number to abc: '02' to 'B', column to A2, A3 harmony like stype
    # # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(int).astype(str)
    # df.insert(loc = 1, column = 'well', value = well_value)

    # change to absolute path
    df['directory'] = [os.path.abspath(f) for f in df['directory']]
    # to unix style
    df['directory'] = ['/'.join(x.split('\\')) for x in df['directory']]
    # to windows style
    # df['directory'] = df['directory'].replace('\\\\', '\\\\\\\\', regex=True)

    # metadata colnames
    metadata_colnames = ["directory", "prefix", "timepoint", "position", "stack"]
    
    # order images
    df['prefix'] = df['prefix'].astype("str")
    df['position'] = df['position'].astype("int64")
    df = df.sort_values(by=['prefix','position'])
    
    # add position metadata
    if position_metadata is not None:
        # add position_metadata for nikon nd2 exported project
        if os.path.exists(os.path.join(measurement, position_metadata)):
            pos_meta = pd.read_csv(os.path.join(measurement, position_metadata))
            # remove NA rows which were interupted
            pos_meta = pos_meta.dropna(subset = ["P Index"])
            # print(pos_meta[["P Index","Position Name"]])
            
            # check if first row Position Name is not None, then do following steps
            if not pd.isna(pos_meta.iloc[0]["Position Name"]):
                print("add position well matadata")
                # split well and field
                pos_meta_well = pos_meta["Position Name"].str.extract("(?P<row>[A-Z])(?P<column>[0-9]+)#(?P<field>.+)")
                # pos_meta_well["column"] = pos_meta_well["column"].str.pad(2, "left", "0")
                pos_meta_well["well"] = pos_meta_well["row"] + pos_meta_well["column"]
                pos_meta = pd.concat([pos_meta["prefix"], 
                                      pos_meta["P Index"].rename("position"), 
                                      pos_meta_well], axis=1)
                pos_meta = pos_meta[["prefix", "position", "well", "field"]]
                pos_meta['position'] = pos_meta['position'].astype("int64")
                # merge to df       
                df = pd.merge(pos_meta, df, how="right", on=["prefix", "position"])
                
                # add metadata cols
                metadata_colnames.extend(["well", "field"])
            

    ## format df to cellprofiler dataloader
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
                df = df.rename(columns={f'{mask_colname}': f'Image_ObjectsFileName_mask_{mask_colname}'})
        # metadata
        for meta in metadata_colnames:
            df = df.rename(columns={f'{meta}': f'Metadata_{meta}'})
    
    return {'df':df, 
            'metadata_colnames':metadata_colnames, 
            'intensity_colnames':image_colnames, 
            'mask_colnames':mask_colnames}




# for nd2 exported format
def image_to_df_nd2(
        measurement = ".",
        subset_pattern="",
        subdir='images',
        image_suffix='.tiff',
        mask_suffix='.png',
        keep_na_rows=True,
        position_metadata="position_metadata.csv"
        ):
    # measurement = "D:\\Others\\2023-11-13_NMYC_spark_value"
    ## get all images
    if subdir is None:
        files_all = glob.glob(os.path.join(measurement, "*" ))
    else:
        files_all = glob.glob(os.path.join(measurement, subdir, "*"))
        
    # subset images
    # remove Probabilities files
    files_all = [f for f in files_all if not re.search("_Probabilities", f)]
    
    # print(subset_pattern)
    if subset_pattern is not None:
        files_all = [f for f in files_all if re.search(subset_pattern, f)]
        
    # sort files
    files_all = natsorted(files_all)
    # check if any masks exists
    if len(files_all) == 0:
        print("no valid images/masks found, check directory or subset_pattern")
        return None

    ## intensity files
    images = [f for f in files_all if re.search(f'{image_suffix}', f)]
    if len(images) == 0:
        print("no valid images found, please check it")
        return None
    # create df
    image_df = pd.DataFrame({  # "path": files,
        "directory": [os.path.dirname(x) for x in images],
        "filename": [os.path.basename(x) for x in images]})
    # images_df = images_df.iloc[:10]
    # get metadata
    extractor = '(?P<prefix>.*)__t(?P<timepoint>[0-9]{1,})_p(?P<position>[0-9]{1,})_z(?P<stack>[0-9]{1,})_c(?P<channel>[0-9]{1,})(?P<self_generated>.*)' + f'{image_suffix}'
    image_meta = image_df["filename"].str.extract(extractor)
    image_meta["channel"] = image_meta["channel"].astype("str") + image_meta["self_generated"]
    image_meta.pop('self_generated')
    # merge with directory, filename
    image_df = pd.concat([image_df, image_meta], axis=1)
    
    # set index
    # images_df = images_df.reset_index()
    image_df = image_df.set_index(["directory", "prefix", "timepoint", "position", "stack"])
    # add prefix to channel
    image_df["channel"] = "ch" + image_df["channel"]
    # get unique channel names
    image_colnames = image_df['channel'].unique().tolist()
    # reshape channel to wide
    image_df = image_df.pivot(columns="channel", values="filename")
    print(f"they're total {len(image_df)} grouped intensity images: {list(image_df)}")

    # keep NA rows
    if not keep_na_rows:
        image_df = image_df.dropna()
    # reset index
    image_df = image_df.reset_index()
    print(f"they're total {len(image_df)} groups after merging images and masks")

    # # add well data
    # row_abc = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    # # row number to abc: '02' to 'B', column keep A02, A03 cellprofiler like stype
    # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(str)
    # # row number to abc: '02' to 'B', column to A2, A3 harmony like stype
    # # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(int).astype(str)
    # df.insert(loc = 1, column = 'well', value = well_value)

    # change to absolute path
    image_df['directory'] = [os.path.abspath(f) for f in image_df['directory']]
    # to unix style
    image_df['directory'] = ['/'.join(x.split('\\')) for x in image_df['directory']]
    # to windows style
    # df['directory'] = df['directory'].replace('\\\\', '\\\\\\\\', regex=True)

    # metadata colnames
    metadata_colnames = ["directory", "prefix", "timepoint", "position", "stack"]
    
    # order images
    image_df['prefix'] = image_df['prefix'].astype("str")
    image_df['position'] = image_df['position'].astype("int64")
    image_df = image_df.sort_values(by=['prefix','position'])
    
    
    # add position_metadata for nikon nd2 exported project
    if os.path.exists(os.path.join(measurement, position_metadata)):
        pos_meta = pd.read_csv(os.path.join(measurement, position_metadata))
        # remove NA rows which were interupted
        pos_meta = pos_meta.dropna(subset = ["P Index"])
        # print(pos_meta[["P Index","Position Name"]])
        
        # check if first row Position Name is not None, then do following steps
        if not pd.isna(pos_meta.iloc[0]["Position Name"]):
            print("add position well matadata")
            # split well and field
            pos_meta_well = pos_meta["Position Name"].str.extract("(?P<row>[A-Z])(?P<column>[0-9]+)#(?P<field>.+)")
            # pos_meta_well["column"] = pos_meta_well["column"].str.pad(2, "left", "0")
            pos_meta_well["well"] = pos_meta_well["row"] + pos_meta_well["column"]
            pos_meta = pd.concat([pos_meta["prefix"], 
                                  pos_meta["P Index"].rename("position"), 
                                  pos_meta_well], axis=1)
            pos_meta = pos_meta[["prefix", "position", "well", "field"]]
            pos_meta['position'] = pos_meta['position'].astype("int64")
            # merge to df       
            image_df = pd.merge(pos_meta, image_df, how="right", on=["prefix", "position"])
            
            # add metadata cols
            metadata_colnames.extend(["well", "field"])
            
            

    return {'df':image_df, 
            'metadata_colnames':metadata_colnames, 
            'intensity_colnames':image_colnames}





# for imageJ exported format
def image_to_df_ij(
        measurement = ".",
        subset_pattern="",
        subdir='images',
        image_suffix='.tif',
        # mask_suffix='.png',
        keep_na_rows=True
        ):
    # measurement = "D:\\Others\\2023-11-13_NMYC_spark_value"
    ## get all images
    if subdir is None:
        files_all = glob.glob(os.path.join(measurement, "*" ))
    else:
        files_all = glob.glob(os.path.join(measurement, subdir, "*"))
        
    # subset images
    # remove Probabilities files
    files_all = [f for f in files_all if not re.search("_Probabilities", f)]
    
    # print(subset_pattern)
    if subset_pattern is not None:
        files_all = [f for f in files_all if re.search(subset_pattern, f)]
        
    # sort files
    files_all = natsorted(files_all)
    # check if any masks exists
    if len(files_all) == 0:
        print("no valid images/masks found, check directory or subset_pattern")
        return None

    ## intensity files
    images = [f for f in files_all if re.search(f'{image_suffix}', f)]
    if len(images) == 0:
        print("no valid images found, please check it")
        return None
    # create df
    image_df = pd.DataFrame({  # "path": files,
        "directory": [os.path.dirname(x) for x in images],
        "filename": [os.path.basename(x) for x in images]})
    # images_df = images_df.iloc[:10]
    # get metadata
    extractor = '(?P<prefix>.*)_t(?P<time>[0-9]{3,})_c(?P<channel>[0-9]{3,})(?P<self_generated>.*)' + f'{image_suffix}'
    image_meta = image_df["filename"].str.extract(extractor)
    image_meta["channel"] = image_meta["channel"].astype("str") + image_meta["self_generated"]
    image_meta.pop('self_generated')
    # merge with directory, filename
    image_df = pd.concat([image_df, image_meta], axis=1)
    
    # set index
    # images_df = images_df.reset_index()
    image_df = image_df.set_index(["directory", "prefix", "time"])
    # add prefix to channel
    image_df["channel"] = "ch" + image_df["channel"]
    # get unique channel names
    image_colnames = image_df['channel'].unique().tolist()
    # reshape channel to wide
    image_df = image_df.pivot(columns="channel", values="filename")
    print(f"they're total {len(image_df)} grouped intensity images: {list(image_df)}")

    # keep NA rows
    if not keep_na_rows:
        image_df = image_df.dropna()
    # reset index
    image_df = image_df.reset_index()
    print(f"they're total {len(image_df)} groups after merging images and masks")

    # # add well data
    # row_abc = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    # # row number to abc: '02' to 'B', column keep A02, A03 cellprofiler like stype
    # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(str)
    # # row number to abc: '02' to 'B', column to A2, A3 harmony like stype
    # # well_value = [row_abc[x - 1] for x in df['row'].astype(int)] + df['column'].astype(int).astype(str)
    # df.insert(loc = 1, column = 'well', value = well_value)

    # change to absolute path
    image_df['directory'] = [os.path.abspath(f) for f in image_df['directory']]
    # to unix style
    image_df['directory'] = ['/'.join(x.split('\\')) for x in image_df['directory']]
    # to windows style
    # df['directory'] = df['directory'].replace('\\\\', '\\\\\\\\', regex=True)

    # metadata colnames
    metadata_colnames = ['directory','prefix','time']

    return {'df':image_df, 
            'metadata_colnames':metadata_colnames, 
            'intensity_colnames':image_colnames}


