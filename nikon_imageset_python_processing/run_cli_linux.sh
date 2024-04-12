#! /bin/bash


## initial setup
cd "D:/Script/nikon_imageset_python_processing"
conda activate cellpose

## input here
project_dir="D:/Postdoc_Data/2024-02-04_MYC_swapping"
cellpose_model=""


## step1: nd2 -> tiff
python step1_nd2tiff.py --dir $project_dir

