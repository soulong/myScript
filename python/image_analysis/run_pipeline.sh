#! /bin/bash

## ---------------------------------------
conda activate cellpose
## script directory
cd '/mnt/f/GitHub/myScript/python/image_analysis'
## setup
# root_dir="/mnt/c/Users/haohe/Desktop/"
# root_dir="/mnt/d/LGLab/Project_NMN/2025-12-22_H3K27me3_IF_analysis/20251222 old plate IMR90 NMN__2025-12-22T15_05_24-Measurement 1"
root_dir="/mnt/f/workspace/3p51_hoechst" # 2025-12-30_HEK293_p53 # 2025-12-24_H1299_reporter_p53

## ---------------------------------------
# optinal
python run_zprojection.py "$root_dir"
# optinal
python run_split_frame.py "$root_dir" --tile_width 4096 --tile_height 4096

# check dataset
# kwargs='{"remove_na_row":True, "subset_pattern":".*p03-ch"}'
# kwargs='{"remove_na_row":True, "image_extractor":"(?P<field>\\d+)-(?P<treat>\\d+)", "mask_extractor":"(?P<field>\\d+)-(?P<treat>\\d+)_cp_(?P<mask_name>.*)"}'
kwargs='{"remove_na_row":False,"subset_pattern":"z[1234567][159]c","image_suffix":".tif","image_extractor":"(?P<name>.*)_SIM².*_z(?P<stack>.*)c(?P<channel>.*)_ORG(?P<self_generated>.*)"}'

python image_helper.py "$root_dir" --dry_run
python image_helper.py "$root_dir" --dry_run --dataset_kwargs "$kwargs"

# optinal
python run_basic_correction.py "$root_dir" --mode fit --channels ch1 #--dataset_kwargs "$kwargs"
python run_basic_correction.py "$root_dir" --mode transform --channels ch1 #--dataset_kwargs "$kwargs"

# cell segmentation
# python run_cellpose.py "$root_dir" --chan1 ch2 --resize 0.5 --dataset_kwargs "$kwargs" #--overwrite
python run_cellpose.py "$root_dir" --chan1 ch2 --resize 0.5 --model "cpsam_nuclei_SIM_63x" --normalize "percentile:0.1,99.9" --flow_thresh 0.4 --cellprob_thresh -1 --dataset_kwargs "$kwargs"

# optinal
python run_ilastik.py "$root_dir" --image_suffix ".tif" --subset_pattern ".*_z[1234567][159]" --ilastik_proj "${root_dir}/MyProject.ilp"

# prepare dataloader
python image_helper.py $root_dir --cp_dataloder --dataset_kwargs "$kwargs" --overwrite



# histogram analysis
# root_dir="/mnt/d/LGLab/Project_NMN/Raw_images_H3K27me3_IF_NMN/LH 20231107 IMR90 NMN treatment H3K27me3__2023-11-07T15_58_57-Measurement 1"
python analyze_histogram.py "$root_dir" --min_intensity 150 --max_intensity 65535 --normalize mean --n_bins 50 --bin_method geomspace --db_name hist_150_65535_mean_50_geomspace.db --overwrite_db --dataset_kwargs "$kwargs"
python analyze_histogram.py "$root_dir" --min_intensity 100 --max_intensity 6000 --normalize mean --n_bins 50 --bin_method geomspace --db_name hist_100_6000_mean_50_geomspace.db --overwrite_db --dataset_kwargs "$kwargs"



## ---------------------------------------
conda activate cellprofiler

# run cellprfiler
python run_cellprofiler.py $root_dir --cp_project "${root_dir}/by_cell_spot.cpproj"


## ---------------------------------------
## then in R, run cellprofiler_analysis.Rmd
