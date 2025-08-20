
#%%
import numpy as np
import pandas
from tqdm import tqdm
from dataset import ImageDataSet
from segmentation import CellposeSegmentation
from analyzer import IntensityAnalyzer, GranularityAnalyzer, RadialDistributionAnalyzer
import tifffile, cv2
import os
import numpy as np
import pandas as pd



#%% initial dataset
data_dir = '/media/hao/Data/Project_MYC/2025-08-01_MYC_LNCaP_P63_P64'

dataset = ImageDataSet(data_dir)


#%% segmentation
segmentor = CellposeSegmentation(
    channel_names=[['ch1'], ['ch2']], channel_weights=[[1], [1]], channel_merge='mean',
    resize_factor=0.5, mask_name='cell', overwrite_mask=False, gpu_batch_size=64)

segmentor.run(dataset.df, test=None)



#%% measurement
channel_names = dataset.intensity_colnames
print(channel_names)
mask_names = dataset.mask_colnames #[x.replace('cp_masks_', '') for x in dataset.mask_colnames]
print(mask_names)



intensity = IntensityAnalyzer(channel_names)
radialdistribution = RadialDistributionAnalyzer(['ch1'])
granularity = GranularityAnalyzer(['ch1'], 
                                  granularity_spectrum_range=[1,3,4,8,12,16],
                                  object_size=30, subsampling_factor=0.25)

per_image =  dataset.df
per_image['image_uid'] = np.nan
per_object = pd.DataFrame()
for idx, row_data in tqdm(dataset.df.iloc[0:6].iterrows()):

    # per_image =  dataset.df.iloc[idx]

    intensity_dict = {ch: tifffile.imread(os.path.join(row_data['directory'], row_data[ch])) for ch in channel_names}
    msk_dict = {ch: cv2.imread(os.path.join(row_data['directory'], row_data[ch]), cv2.IMREAD_UNCHANGED) for ch in mask_names}

    for msk_idx, msk_name in enumerate(mask_names):
        msk = msk_dict[msk_name]
        measure_global = True if msk_idx == 0 else False

        measure_arrays = np.stack([intensity_dict[k] for k in intensity.channel_names], axis=0)
        intensity_result = intensity.measure(measure_arrays, msk, measure_global=measure_global)
        
        measure_arrays = np.stack([intensity_dict[k] for k in radialdistribution.channel_names], axis=0)
        radialdistribution_result = radialdistribution.measure(measure_arrays, msk)
        
        measure_arrays = np.stack([intensity_dict[k] for k in granularity.channel_names], axis=0)
        granularity_result = granularity.measure(measure_arrays, msk)

        if measure_global:
            per_image.loc[idx, list(intensity_result['global'].columns)] = intensity_result['global'].values
            per_image.loc[idx, 'image_uid'] = idx
        
        if 'objects' in list(intensity_result.keys()):
            obj = intensity_result['objects']
        else:
            obj = pd.DataFrame()
        if intensity_result:
            obj = pd.merge(obj, radialdistribution_result['radialdistribution'], on='object_id')
        if granularity_result:
            obj = pd.merge(obj, granularity_result['granularity'], on='object_id')
        
        obj.insert(0, 'image_uid', idx, allow_duplicates=False)

        # write to sqlite
        

        # if idx == 0:
        #     per_object = pd.DataFrame(columns=obj.columns)

        # per_object.loc[idx] = obj.values
        per_object = pd.concat([per_object, obj], axis=0)


