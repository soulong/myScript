#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 23:10:07 2025

@author: hao
"""

import os 
os.chdir('/media/hao/Data1/single_cell_dataset/two_channel_293T_barcode/2025-01-06_BC_test')
import numpy as np
import re
from natsort import natsorted
import pickle
import sys
sys.path.append('/home/hao/Documents/GitHub/myScript/python_functions')
from helper_dl import SingCell, visualize_dataloader, embed_umap, extract_feature
import imageio.v3 as iio
import timm
import torch
from torch import nn
import torch.optim as optim
from torch.utils.data import DataLoader
# import torchvision.transforms.v2 as T
import albumentations as A
import pytorch_lightning as pl
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split


#%% transform
img_size =128

train_transform = A.Compose([
    A.LongestMaxSize(max_size=img_size),
    A.PadIfNeeded(min_height=img_size, min_width=img_size),
    A.HorizontalFlip(p=0.5),
    A.Rotate(p=0.5),
    A.Blur(blur_limit=5),
    # A.GaussNoise(),
    A.Affine(scale=(0.9, 1.1), translate_percent=0.1, shear=15),
    # A.RandomBrightnessContrast(p=0.2),
    A.ToTensorV2(),
    ], seed=42, strict=True)

val_transform = A.Compose([
    A.LongestMaxSize(max_size=img_size),
    A.PadIfNeeded(min_height=img_size, min_width=img_size),
    A.ToTensorV2()
    ])


#%% dataset
def prepare_classifier_data(
        data_dir, 
        include_pattern = None,
        exclude_pattern = None,
        size = 128,
        transform = [None, None], # train ,val
        normalize: str = "percentile",
        img_suffix = ".tiff",
        max_pixel_value: int = 65535,
        val = 0.2, # val split ratio
        batch_size = 256
        ):
    # recursively find all TIFF images within the root directory
    image_paths = []
    for dirpath, dirnames, filenames in os.walk(data_dir, followlinks=False):
        for filename in filenames:
            if filename.endswith(img_suffix):
                image_paths.append(os.path.join(dirpath, filename))
    if include_pattern is not None:
        print(f"include files by {include_pattern}")
        image_paths = [f for f in image_paths if re.match(include_pattern, f)]
    if exclude_pattern is not None:
        print(f"exclude files by {exclude_pattern}")
        image_paths = [f for f in image_paths if not re.match(exclude_pattern, f)]
    if not image_paths:
        raise FileNotFoundError(f"No TIFF images found in directory: {data_dir}")
    
    # files split
    # get label for dataset split
    labels = [os.path.basename(os.path.dirname(f)) for f in image_paths]
    train_files, val_files, _, _ = train_test_split(
        image_paths, labels, test_size=val, random_state=42, stratify=labels)
    
    # dataset
    train_dataset = SingCell(
        image_paths=train_files, normalize=normalize, transform=transform[0])
    val_dataset = SingCell(
        image_paths=val_files, normalize=normalize, transform=transform[1])
    
    # loader
    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True, num_workers=20)
    val_loader = DataLoader(
        val_dataset, batch_size=batch_size, shuffle=False, num_workers=20)

    return train_loader, val_loader, train_dataset, val_dataset


train_loader, val_loader, train_dataset, val_dataset = prepare_classifier_data(
        data_dir='/media/hao/Data1/single_cell_dataset/single_channel_128/2025-01-06_BC_test', 
        include_pattern = None, exclude_pattern = '.*Mix.*',
        size = 128, transform = [train_transform, val_transform], # only for train
        normalize = None,
        img_suffix = ".tiff", max_pixel_value = 65535,
        val = 0.2, batch_size = 256)

# viewe samples
if True:
    test_samples = next(iter(train_loader))
    # x = test_samples[0][0,:,:,:].numpy()
    visualize_dataloader(test_samples[0], test_samples[1])
   

#%% model
class CellClassifier(pl.LightningModule):
    def __init__(self, model_name='resnet18', in_chans=1, num_classes=2, learning_rate=1e-3):
        super().__init__()
        self.save_hyperparameters() # Automatically save hyperparameters
        self.model = timm.create_model(
            model_name, pretrained=True, in_chans=in_chans, num_classes=num_classes)
        self.criterion = nn.CrossEntropyLoss()

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        images, labels, *_ = batch
        outputs = self(images)
        loss = self.criterion(outputs, labels)
        self.log('train_loss', loss, prog_bar=True, on_epoch=True, on_step=True)
        return loss

    def validation_step(self, batch, batch_idx):
        images, labels, *_ = batch
        outputs = self(images)
        loss = self.criterion(outputs, labels)
        acc = (outputs.argmax(dim=-1) == labels).float().mean()
        self.log('val_loss', loss, prog_bar=True, on_epoch=True, on_step=True)
        self.log('val_acc', acc, prog_bar=True, on_epoch=True, on_step=True)
        return loss

    def configure_optimizers(self):
        return optim.Adam(self.parameters(), lr=self.hparams.learning_rate)


in_chans = train_dataset[0][0].shape[0]
num_classes = len(train_dataset.idx_to_label)

model = CellClassifier(
    model_name='resnet34', 
    in_chans=in_chans, 
    num_classes=num_classes, 
    learning_rate=1e-3)


#%% train
trainer = pl.Trainer(max_epochs=100, precision='16')
trainer.fit(model, train_loader, val_loader)

v_num = 6
# save model
torch.save(model, f'lightning_logs/version_{v_num}/model.pth')

#%% extract features
# load model
model = torch.load('lightning_logs/version_6/model.pth')
extractor = nn.Sequential(*list(model.model.children())[:-1])
# features = feature_extractor(test_samples[0]).detach().numpy()

# dataloader
loader, *_ = prepare_classifier_data(
    data_dir='/media/hao/Data1/single_cell_dataset/single_channel_128/2025-01-06_BC_test', 
    include_pattern = None, exclude_pattern = '.*Mix.*',
    size = 128, transform = [train_transform, val_transform], # only for train
    normalize = None, img_suffix = ".tiff", max_pixel_value = 65535,
    val = 0.0, batch_size = 256)

# get features
features = extract_feature(extractor, loader, 
                           output_umap_enbed=True, 
                           plot=True,
                           output_encode_img=True)
# save features
fname = f'lightning_logs/version_3/features.pkl'
with open(fname, 'wb') as file:
    pickle.dump(features, file)
    

# check feature embeding
# python '/home/hao/Documents/GitHub/myScript/python_functions/vis_features_app.py' 






