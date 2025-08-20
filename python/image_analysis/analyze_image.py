#%% import
# env: pytorch
import os, sys, warnings
sys.path.append('/home/hao/Documents/GitHub/myScript/python/image_analysis')
warnings.filterwarnings('ignore')
from image_helper import (
    nd2_to_tiff, images_to_dataset, cellpose_segment_dataset, 
    ilastik_multiple, spot_segment_dataset, 
    measure_dataset, measure_dataset_by_model,
    radial_distribution, granularity, crop_cell_from_dir)

#%% data_dir
data_dir = """
/media/hao/Data/Project_ZM_translocation_reporter/20240924 ZM translocation Casp3+TEV reporter
"""
data_dir = data_dir.strip()

#%% nd2 to tiff
nd2_to_tiff(data_dir, resize=None, crop_N=None, only_meta=False) # resize=[512,512]


# %% flatfield correction
# if False:
#    from step2_shading_correction import image_correction
#    res0 = measurement_to_df_nd2(
#        measurement=measurement_dir,
#        subset_pattern=None,
#        image_subdir="images",
#        image_suffix='.tiff',
#        mask_suffix=None,  # None for no use
#        keep_na_rows=True,
#        cellprofiler_style=False,
#        position_metadata=None)
#    # fit first
#    image_correction(res0, run_type="f", process_channels=["ch2"],
#                     n_image=15, resize=128, enable_dfd=False)
#    # check fitted result
#    # then transform
#    image_correction(res0, run_type="t", process_channels=["ch1"],
#                     target_dir="test")

# %% channel registration
if False:
    from skimage.registration import phase_cross_correlation
    from scipy.ndimage import shift
    from glob import glob
    from re import sub
    from tifffile import imread, imwrite
    import numpy as np
    
    f = glob(data_dir + "/**/*.tif")
    prefix = list(set([sub("\(Alexa.*", "", x) for x in f]))
    channel = list(set([sub(".*\(Alexa", "", x) for x in f]))
    print(channel)
    
    for p in prefix:
        img = [imread(p + "(Alexa" + c) for c in channel]
        img_ref = img[0]
        img_offset = img[1]
        shift_yx, error, diffphase = phase_cross_correlation(img_ref, img_offset)
        aligned = shift(img_offset, shift=shift_yx)
        imwrite(p + "(Alexa" + " aligned" + channel[1], aligned)
        
        # import matplotlib.pyplot as plt
        # plt.figure(figsize=(10, 5))
        # plt.subplot(1, 3, 1)
        # plt.title("Channel 1")
        # plt.imshow(img_ref, cmap='gray')
        # plt.subplot(1, 3, 2)
        # plt.title("Original Channel 2")
        # plt.imshow(img_offset, cmap='gray')
        # plt.subplot(1, 3, 3)
        # plt.title("Aligned Channel 2")
        # plt.imshow(aligned, cmap='gray')
        # plt.show()

#%% cell segmentation
dataset = images_to_dataset(data_dir, subset_pattern=None, remove_na_row=True,
    image_subdir='.', image_suffix=".tif", cellprofiler_style=True,
    image_extractor='(?P<prefix>.*)(?P<self_generated>.*)',
    mask_extractor='(?P<prefix>.*)_cp_masks_(?P<mask_name>.*)')

# dataset = images_to_dataset(data_dir, subset_pattern=None, remove_na_row=True)
# dataset["df"].to_csv(os.path.join(data_dir, "cp_dataloader.csv"), index=False)

# for cell
cellpose_segment_dataset(
    dataset, 
    channel_names = [['ch1'], ], 
    channel_weights = [[1,], ],
    channel_merge = "mean", # mean or sum
    model_name = "cpsam", # cyto3, cyto3_M7, cyto3_spot, cpsam_EM_mito
    diameter = None, 
    normalize = {'percentile': [0.1, 99.9]},
    resize_factor = 1, 
    mask_name = 'cell',
    reduce_mask_name = True, overwrite_mask = False,
    device = 'gpu', gpu_batch_size = 32)

# for nuclei
cellpose_segment_dataset(
    dataset, 
    channel_names = [['ch0',], ], 
    channel_weights = [[1,], ],
    channel_merge = "mean", # mean or sum
    model_name = "cpsam_EM_mito", # cyto3, cyto3_M7, cyto3_spot, cpsam_EM_mito
    diameter = None, 
    normalize = {'percentile': [0.1, 99.9]},
    resize_factor = 0.5, 
    mask_name = 'mito',
    reduce_mask_name = True, overwrite_mask = False,
    device = 'gpu', gpu_batch_size = 16)

#%% ilastik segmentation (probability)
ilastik_multiple(
    measurement=data_dir,
    ilastik_proj=f"{data_dir}/MyProject.ilp",
    # ilastik_proj="/media/hao/Data/Project_dRTK/2024-06-28_dRTK_SPARK2_w-mKO3_add_dALK_dEGFR-Mut/MyProject.ilp",
    subset_pattern=".*_ch[1].tiff",
    image_subdir="images",
    image_suffix=".tiff",
    output_suffix="Probabilities",
    ilastik_exec="/home/hao/software/ilastik/run_ilastik.sh",
    thread_num=36,
    overwrite_mask=False,
    test=False)

#%% spot segmentation (probability)
dataset = images_to_dataset(data_dir, subset_pattern=None, remove_na_row=False)
spot_segment_dataset(
    dataset, 
    channel_names = ['ch1','ch4'], 
    pred_names = ['spot_pred_ch1.png','spot_pred_ch4.png'], 
    overlap = 0.0, threshold = 0.5,
    model_path = "/media/hao/Data1/single_cell_dataset/cell_spot_cnn_training_data/logs/v3/model.pth",
    model_input_size = 512, model_batch_size = 16, 
    overwrite_spot = True)


#%% dataoader for cellprofiler
res_cp = images_to_dataset(data_dir, remove_na_row=True, cellprofiler_style=True)
res_cp["df"].to_csv(os.path.join(data_dir, "cp_dataloader.csv"), index=False)

"""
## run via terminal
# conda activate cellprofiler
f'(echo "sleep for N seconds"; sleep 2) && python /home/hao/Documents/GitHub/myScript/python/image_analysis/step5_run_cellprofiler.py \
    --measurements "{data_dir}" \
    --cp_project_path "{data_dir}/by_cellRegion.cpproj" \
    --thread_num 32'
"""


#%% cell measurement
dataset = images_to_dataset(data_dir, subset_pattern=None, remove_na_row=True)
# dataset['df'] = dataset['df'].iloc[1:20] # for test a few images
_ = measure_dataset(dataset, 
                    save_path = os.path.join(data_dir, 'cell.db'),
                    # extra_properties = [radial_distribution,]
                    # extra_properties = [radial_distribution, granularity,]
                    )



#%% test radial distribution and granularity

# import imageio.v3 as iio
# from matplotlib import pyplot as plt
# img = iio.imread('/media/hao/Data/Project_MYC/2024-10-31_MYCN_saturation_curve_adding_truncation/images/60x-293t-50ng-MYCN+N-146-100ng-30ng-10ng-3ng-nd2__t1_p25_z1_c1.tiff')
# mask = iio.imread('/media/hao/Data/Project_MYC/2024-10-31_MYCN_saturation_curve_adding_truncation/images/60x-293t-50ng-MYCN+N-146-100ng-30ng-10ng-3ng-nd2__t1_p25_z1_c1_cp_masks_cell.png')
# img = img/65535
# plt.imshow(img)
# plt.imshow(mask)


#%% common scripts
if False:  
    # crop and save single cell
    crop_cell_from_dir(data_dir, intensity_channel_name=['ch1'], 
                       target_size=256, rotation=None,
                       max_image_per_well=2, max_cell_per_image=100,
                       save_multi_channel_tiff=False)
