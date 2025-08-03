#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 02:15:32 2025

@author: hao
"""

import os
import sys
sys.path.append('/home/hao/Documents/GitHub/myScript/nikon_imageset_python_processing')
from functions import crop_cell_from_dir
from os.path import basename, dirname
import re
import glob
import shutil


# setup
source_dir = '/media/hao/Data/Project_barcode/2025-01-06_BC_test'

target_dir = '/media/hao/Data1/single_cell_dataset/single_channel_128/2025-01-06_BC_test'

os.chdir(target_dir)



# step1: crop cell
crop_cell_from_dir(source_dir, target_dir,
                   max_image_per_well=4, max_cell_per_image=200,
                   intensity_channel_name=['ch1','ch4'], 
                   target_size=128)



# step2: merge rows and clean
files = glob.glob('./**/*.tiff', recursive=True)
rows = [re.search('[a-zA-z]', basename(dirname(file)).split('__')[1]) for file in files]
rows = [x.group() for x in rows if x]

mapping = {'B': 'BC100', 'C': 'BC101', 'D': 'BC102', 'E': 'BC103',
           'F': 'BC104', 'G': 'BC105', 'H': 'BC110', 'I': 'Mix-1', 'J': 'Mix-2'}

_ = [os.makedirs(mapping[x], exist_ok=True) for x in set(rows)]

for idx, file in enumerate(files):
    shutil.move(file, os.path.join(mapping[rows[idx]], basename(file)))
    pass