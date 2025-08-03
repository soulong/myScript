#!/usr/bin/env python3
"""
Semantic segmentation model for cell spot detection.
Uses U-Net architecture with ResNet34 encoder.
"""

#%%
import os
import glob
import random
import numpy as np
import torch
# import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import albumentations as A
from albumentations.pytorch import ToTensorV2
import cv2
import tifffile
import segmentation_models_pytorch as smp
from torch.utils.tensorboard import SummaryWriter
import matplotlib.pyplot as plt



class Config:
    def __init__(self, 
                 data_root, 
                 checkpoint_dir=None, 
                 log_dir=None, 
                 image_size=512, 
                 in_channels=1, 
                 classes=1, 
                 val_split=0.2,                  
                 batch_size=64, 
                 num_epochs=500, 
                 learning_rate=1e-3, 
                 weight_decay=1e-4, 
                 encoder='resnet34', 
                 encoder_weights='imagenet', 
                 device=torch.device("cuda" if torch.cuda.is_available() else "cpu"), 
                 num_workers=10):
        
        # Data paths
        self.data_root = data_root
        self.checkpoint_dir = checkpoint_dir if checkpoint_dir else data_root + '/logs'
        self.log_dir = log_dir if log_dir else checkpoint_dir

        # Training params
        self.batch_size = batch_size
        self.num_epochs = num_epochs
        self.image_size = image_size
        self.val_split = val_split
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay

        # Model params
        self.encoder = encoder
        self.encoder_weights = encoder_weights
        self.in_channels = in_channels
        self.classes = classes

        # Hardware
        self.device = device
        self.num_workers = num_workers


class SegmentationDataset(Dataset):
    def __init__(self, image_paths, mask_paths, transform=None):
        self.image_paths = image_paths
        self.mask_paths = mask_paths
        self.transform = transform

    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # Load image and normalize to [0,1]
        image = tifffile.imread(self.image_paths[idx])
        assert image.ndim == 2, f"Image has {image.ndim} dimensions, expected 2"
        image = (image / 65535).astype(np.float32)
        
        # Load mask
        mask = cv2.imread(self.mask_paths[idx], cv2.IMREAD_UNCHANGED).astype(np.uint16)
        # it's somehow specifical here that the mask is 2 for background
        mask[mask == 2] = 0
    
        # Apply augmentations
        if self.transform:
            augmented = self.transform(image=image, mask=mask)
            image, mask = augmented["image"], augmented["mask"]
        
        return image, mask, self.image_paths[idx]


def get_transforms(image_size):
    """Get training and validation transforms"""
    train_transform = A.Compose([
        A.PadIfNeeded(image_size, image_size),
        A.RandomCrop(image_size, image_size),
        # A.Resize(image_size, image_size),
        A.Rotate(),
        A.Affine(scale=(0.8, 1.2), translate_percent=0.2),
        A.RandomBrightnessContrast(brightness_limit=0.2, contrast_limit=0.2, p=0.2),
        # A.ElasticTransform(alpha=120, sigma=120 * 0.05, alpha_affine=120 * 0.03),
        ToTensorV2()
    ])

    val_transform = A.Compose([
        A.PadIfNeeded(image_size, image_size),
        A.RandomCrop(image_size, image_size),
        # A.Resize(image_size, image_size),
        ToTensorV2()
    ])
    
    return train_transform, val_transform


def prepare_data(config):
    """Prepare train and validation dataloaders"""
    # Find all image/mask pairs
    mask_paths = sorted(glob.glob(os.path.join(config.data_root, "**", "*.png"), recursive=True))
    image_paths = [p.replace("_Simple Segmentation.png", ".tiff") for p in mask_paths]
    
    # Split data
    dataset_size = len(image_paths)
    val_size = int(dataset_size * config.val_split)
    train_size = dataset_size - val_size
    indices = list(range(dataset_size))
    random.shuffle(indices)
    train_indices, val_indices = indices[:train_size], indices[train_size:]

    # Get transforms
    train_transform, val_transform = get_transforms(config.image_size)

    # Create datasets
    train_dataset = SegmentationDataset(
        [image_paths[i] for i in train_indices],
        [mask_paths[i] for i in train_indices],
        transform=train_transform
    )
    val_dataset = SegmentationDataset(
        [image_paths[i] for i in val_indices],
        [mask_paths[i] for i in val_indices],
        transform=val_transform
    )

    # Create dataloaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=config.batch_size,
        shuffle=True,
        num_workers=config.num_workers
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=config.batch_size,
        shuffle=False,
        num_workers=config.num_workers
    )

    return train_loader, val_loader, train_transform, val_transform


def prepare_model(config):
    """Prepare model, criterion and optimizer"""
    # Initialize model
    model = smp.Unet(
        encoder_name=config.encoder,
        encoder_weights=config.encoder_weights,
        in_channels=config.in_channels,
        classes=config.classes
    ).to(config.device)
    
    # Loss and optimizer
    criterion = smp.losses.DiceLoss(mode='binary')
    # criterion = smp.losses.FocalLoss(mode='binary')
    optimizer = optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay
    )
    
    return model, criterion, optimizer


def visualize_augmentations(dataset, num_samples=4, cols=2):
    """Visualize augmented training samples"""
    rows = int(np.ceil(num_samples / cols))
    fig, axes = plt.subplots(rows, cols * 2, figsize=(10, 5 * rows))
    axes = axes.ravel()
    
    for idx in range(num_samples):
        image, mask, _ = dataset[idx]
        if isinstance(image, torch.Tensor):
            image = image.numpy().transpose(1, 2, 0)
        if isinstance(mask, torch.Tensor):
            mask = mask.numpy()
            
        ax_idx = idx * 2
        axes[ax_idx].imshow(image, cmap='gray')
        axes[ax_idx].set_title(f'Augmented Image {idx+1}')
        axes[ax_idx].axis('off')
        
        axes[ax_idx + 1].imshow(mask, cmap='jet', alpha=0.5)
        axes[ax_idx + 1].set_title(f'Augmented Mask {idx+1}')
        axes[ax_idx + 1].axis('off')
    
    plt.tight_layout()
    return fig


class SegmentationTrainer:
    def __init__(self, config):
        self.config = config
        self.best_val_loss = float("inf")
        self.writer = SummaryWriter(config.log_dir)
        
        # Create directories
        os.makedirs(config.checkpoint_dir, exist_ok=True)
        os.makedirs(config.log_dir, exist_ok=True)
        
        # Setup data and model
        self.train_loader, self.val_loader, self.train_transform, self.val_transform = prepare_data(config)
        self.model, self.criterion, self.optimizer = prepare_model(config)

    def save_checkpoint(self, epoch, val_loss):
        if val_loss < self.best_val_loss:
            self.best_val_loss = val_loss
            checkpoint = {
                'epoch': epoch + 1,
                'model': self.model,
                # 'optimizer_state_dict': self.optimizer.state_dict(),
                'val_loss': val_loss,
                'train_transform': self.train_transform,
                'val_transform': self.val_transform,
                'config': self.config.__dict__,
                'criterion': self.criterion
            }
            checkpoint_path = os.path.join(self.config.checkpoint_dir, "model.pth")
            torch.save(checkpoint, checkpoint_path)
            print(f"Checkpoint saved: {checkpoint_path}")

    def visualize_results(self, epoch, images, masks, outputs):
        N = min(8, len(images))
        images_np = images.cpu().numpy()[:N, 0]
        masks_np = masks.cpu().numpy()[:N, 0]
        outputs_np = outputs.cpu().numpy()[:N, 0]
        
        fig, axes = plt.subplots(N, 3, figsize=(10, 3 * N))
        for i in range(N):
            axes[i, 0].imshow(images_np[i], cmap="gray")
            axes[i, 0].set_title("Image")
            axes[i, 1].imshow(masks_np[i], cmap="jet", alpha=0.5)
            axes[i, 1].set_title("Ground Truth")
            axes[i, 2].imshow(outputs_np[i], cmap="jet", alpha=0.5)
            axes[i, 2].set_title("Prediction")
            for ax in axes[i]:
                ax.axis("off")
        
        self.writer.add_figure("Validation Results", fig, epoch)
        plt.close(fig)

    def train_epoch(self):
        self.model.train()
        train_loss = 0.0
        for images, masks, _ in self.train_loader:
            images = images.to(self.config.device)
            masks = masks.to(self.config.device).float().unsqueeze(1)
            self.optimizer.zero_grad()
            outputs = self.model(images)
            loss = self.criterion(outputs, masks)
            loss.backward()
            self.optimizer.step()
            train_loss += loss.item()
            
        return train_loss / len(self.train_loader)

    def validate_epoch(self):
        self.model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for images, masks, _ in self.val_loader:
                images = images.to(self.config.device)
                masks = masks.to(self.config.device).float().unsqueeze(1)
                outputs = self.model(images)
                loss = self.criterion(outputs, masks)
                val_loss += loss.item()
                
        return val_loss / len(self.val_loader)

    def train(self):
        print(f"Training on device: {self.config.device}")
        print(f"Total epochs: {self.config.num_epochs}")
        print(f"checkpoint_dir: {self.config.checkpoint_dir}")
        print(f"log_dir: {self.config.log_dir}")
        
        for epoch in range(self.config.num_epochs):
            train_loss = self.train_epoch()
            val_loss = self.validate_epoch()
            
            # Logging
            self.writer.add_scalar("Loss/Train", train_loss, epoch)
            self.writer.add_scalar("Loss/Validation", val_loss, epoch)
            
            # Save checkpoint
            self.save_checkpoint(epoch, val_loss)
            
            # Visualize results periodically
            if epoch % 10 == 0:
                images, masks, _ = next(iter(self.val_loader))
                images = images.to(self.config.device)
                masks = masks.to(self.config.device).float().unsqueeze(1)
                outputs = torch.sigmoid(self.model(images))
                outputs = (outputs > 0.5).float()
                self.visualize_results(epoch, images, masks, outputs)
            
            print(f"Epoch {epoch+1}/{self.config.num_epochs}, "
                  f"Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}")
        
        self.writer.close()
        torch.cuda.empty_cache()


def main(**kwargs):
    config = Config(**kwargs)
    trainer = SegmentationTrainer(config)
    
    # Visualize some augmented training samples before training
    train_dataset = trainer.train_loader.dataset
    aug_fig = visualize_augmentations(train_dataset)
    plt.show()
    
    trainer.train()






#%%
if __name__ == '__main__':

    config = Config(
        data_root = '/media/hao/Data1/spot_segmentation_traindata',
        checkpoint_dir = '/media/hao/Data1/spot_segmentation_traindata/logs_v3',
        log_dir = None,
        image_size = 512,
        in_channels = 1,
        classes = 1,
        val_split = 0.2,
        batch_size = 16,
        num_epochs = 500,
        learning_rate = 1e-3,
        weight_decay = 1e-4,
        encoder = 'resnet34',
        encoder_weights = 'imagenet',
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu"),
        num_workers = 10)
    
    trainer = SegmentationTrainer(config)

    train_dataset = trainer.train_loader.dataset
    next_batch = next(iter(trainer.train_loader))
    print(next_batch[0].shape)
    print(next_batch[1].shape)

    aug_fig = visualize_augmentations(train_dataset)
    plt.show()
    
    trainer.train()

