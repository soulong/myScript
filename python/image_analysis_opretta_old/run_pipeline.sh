#! /bin/bash

## script directory
cd '/mnt/c/Users/haohe/Documents/GitHub/myScript/python/image_analysis_opretta'

## setup
root_dir="/mnt/c/Users/haohe/Desktop/"
ilastik_proj="${root_dir}/p53_spot.ilp"
cp_project="${root_dir}/by_cell_spot.cpproj"

# ilastik_proj="/mnt/d/TP53/HCS/p53_spot.ilp"
# cp_project="/mnt/d/TP53/HCS/by_cell_spot.cpproj"


## ---------------------------------------
conda activate cellpose

# optinal
python run_zprojection.py $root_dir

# optinal
python run_split_frame.py $root_dir -w 2048 -h 2048

# optinal
python run_basic_correction.py $root_dir --mode fit --channels ch1
python run_basic_correction.py $root_dir --mode transform

# cell segmentation
python run_cellpose.py $root_dir --diameter 160 --chan1 ch1

# optinal
python run_ilastik.py $root_dir --ilastik-proj $ilastik_proj


## ---------------------------------------
conda activate cellprofiler

# prepare dataloader
python run_dataloader.py $root_dir

# run cellprfiler
python run_cellprofiler.py $root_dir --cp-project $cp_project


## ---------------------------------------
## then in R, run cellprofiler_analysis.Rmd
