#! /bin/bash

## ---------------------------------------
conda activate cellpose
## script directory
cd '/mnt/c/Users/haohe/Documents/GitHub/myScript/python/image_analysis'
## setup
# root_dir="/mnt/c/Users/haohe/Desktop/"
# root_dir="/mnt/d/LGLab/Project_NMN/2025-12-22_H3K27me3_IF_analysis/20251222 old plate IMR90 NMN__2025-12-22T15_05_24-Measurement 1"
root_dir="/mnt/d/LGLab/Project_NMN/2025-12-22_H3K27me3_IF_analysis"

## ---------------------------------------
# optinal
python run_zprojection.py "$root_dir"
# optinal
python run_split_frame.py "$root_dir" --tile_width 4096 --tile_height 4096

# check dataset
python image_helper.py "$root_dir" --dry_run

# kwargs='{"remove_na_row":False,"image_subdir":"Images","subset_pattern":".*_z[123][13579]_ORG--t","image_suffix":".tif","image_extractor":"(?P<celltype>.*) (?P<treat>.*) SIR DNA_SIM²_z(?P<stack>.*)_ORG--t(?P<field>\\d+)_r\\d+c\\d+_?(?P<self_generated>.*)"}'
kwargs='{"remove_na_row":True, "subset_pattern":".*p02-ch"}'
# kwargs='{"remove_na_row":True, "image_extractor":"(?P<field>\\d+)-(?P<treat>\\d+)", "mask_extractor":"(?P<field>\\d+)-(?P<treat>\\d+)_cp_(?P<mask_name>.*)"}'
python image_helper.py "$root_dir" --dry_run --dataset_kwargs "$kwargs"
# IMR90 Control SIR DNA_SIM²_z01_ORG.tif
# MEF Scramble KO NMN_3D_1_SIM²_z40_ORG.tif

# optinal
python run_basic_correction.py "$root_dir" --mode fit --channels ch1 #--dataset_kwargs "$kwargs"
python run_basic_correction.py "$root_dir" --mode transform #--dataset_kwargs "$kwargs"

# cell segmentation
python run_cellpose.py "$root_dir" --chan1 ch1 --resize 0.5 --dataset_kwargs "$kwargs" #--overwrite
# python run_cellpose.py "$root_dir" --chan1 ch0 --diameter 500 --resize 0.5 --model "cpsam_nuclei_SIM_63x" --normalize "percentile:0.1,99.9" --flow_thresh 0.4 --cellprob_thresh -1 --dataset_kwargs "$kwargs"

# optinal
python run_ilastik.py "$root_dir" --image_suffix ".tif" --subset_pattern ".*_z[123][13579]_ORG--t" --ilastik_proj "${root_dir}/../nulei_sim.ilp"

# prepare dataloader
python image_helper.py $root_dir --cp_dataloder --dataset_kwargs "$kwargs" --overwrite



# histogram
# root_dir="/mnt/d/LGLab/Project_NMN/Raw_images_H3K27me3_IF_NMN/LH 20231107 IMR90 NMN treatment H3K27me3__2023-11-07T15_58_57-Measurement 1"
python analyze_histogram.py "$root_dir" --min_intensity 150 --max_intensity 65535 --normalize mean --n_bins 50 --bin_method geomspace --db_name hist_150_65535_mean_50_geomspace.db --overwrite_db --dataset_kwargs "$kwargs"
python analyze_histogram.py "$root_dir" --min_intensity 150 --max_intensity 20000 --normalize mean --n_bins 50 --bin_method linear --db_name hist_150_20000_mean_50_linear.db --overwrite_db --dataset_kwargs "$kwargs"
python analyze_histogram.py "$root_dir" --min_intensity 150 --max_intensity 20000 --normalize mean --n_bins 50 --bin_method linear --db_name hist_150_20000_mean_50_linear.db --overwrite_db --dataset_kwargs "$kwargs"



## ---------------------------------------
conda activate cellprofiler

# run cellprfiler
python run_cellprofiler.py $root_dir --cp_project "${root_dir}/by_cell_spot.cpproj"


## ---------------------------------------
## then in R, run cellprofiler_analysis.Rmd
