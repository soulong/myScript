

import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
import imageio.v3 as iio
import torch
import torch.nn as nn
import torchvision
import pytorch_lightning as pl
from torch.utils.data import Dataset, DataLoader, random_split, Subset
from torchvision.models.detection import maskrcnn_resnet50_fpn, MaskRCNN_ResNet50_FPN_Weights
from pytorch_lightning.loggers import TensorBoardLogger
from torchvision.transforms import v2 as T

# --- Data Transforms ---
def get_transform(train=True):
    transforms = []
    transforms.append(T.Resize((256, 256), antialias=True))
    
    if train:
        transforms.extend([
            T.RandomHorizontalFlip(p=0.5),
            T.RandomVerticalFlip(p=0.5),
            T.RandomRotation(30),
            T.RandomAdjustSharpness(2, p=0.3),
            T.GaussianBlur(3, p=0.3)
        ])
    
    transforms.extend([
        T.ConvertImageDtype(torch.float32),
        T.Normalize(mean=[0], std=[65535.0])
    ])
    
    return T.Compose(transforms)

# --- Dataset Class ---
class InstanceSegmentationDataset(Dataset):
    def __init__(self, root_dir, train=True):
        self.root_dir = root_dir
        self.train = train
        self.transform = get_transform(train)  # Pass train flag to get appropriate transforms
        self.image_paths = []
        self.mask_paths = []

        # Recursively find all image/mask pairs in subdirectories
        for dirpath, _, filenames in os.walk(root_dir):
            for filename in filenames:
                if filename.lower().endswith(('.tiff', '.tif')):
                    image_path = os.path.join(dirpath, filename)
                    mask_name = os.path.splitext(filename)[0] + "_spot.png"
                    mask_path = os.path.join(dirpath, mask_name)
                    
                    if os.path.exists(mask_path):
                        self.image_paths.append(image_path)
                        self.mask_paths.append(mask_path)

    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # Load image and mask
        img = iio.imread(self.image_paths[idx])
        img = img[:,:,None] # Add channel dimension
        img = torch.from_numpy(img).permute(2, 0, 1)  # Convert to CxHxW format
        
        mask = iio.imread(self.mask_paths[idx])
        mask = torch.from_numpy(mask)

        # Get instance objects (excluding background)
        objs = torch.unique(mask)
        objs = objs[objs != 0]
        
        # Create binary masks and bounding boxes for each instance
        masks = [(mask == obj).to(torch.uint8) for obj in objs]
        bboxes = []
        for msk in masks:
            msk_np = msk.numpy()
            x, y, w, h = cv2.boundingRect(cv2.findNonZero(msk_np))
            bboxes.append([x, y, x + w, y + h])
            
        # All instances are same class (1)
        labels = torch.ones(len(objs), dtype=torch.int64)
        
        target = {
            "masks": torch.stack(masks) if masks else torch.zeros((0, mask.shape[0], mask.shape[1]), dtype=torch.uint8),
            "boxes": torch.tensor(bboxes, dtype=torch.float32) if bboxes else torch.zeros((0, 4), dtype=torch.float32),
            "labels": labels
        }

        # Apply transforms
        if self.transform:
            img = self.transform(img)

        return img, target



# --- Model Architecture ---
def get_model_instance_segmentation(in_chans, num_classes):
    # Load pre-trained model
    model = maskrcnn_resnet50_fpn(weights=MaskRCNN_ResNet50_FPN_Weights.DEFAULT)
    
    # Modify first conv layer for single channel input if needed
    in_channels = model.backbone.body.conv1.in_channels
    if in_channels != in_chans:
        model.backbone.body.conv1 = nn.Conv2d(in_chans, 64, kernel_size=(7, 7), 
                                            stride=(2, 2), padding=(3, 3), bias=False)

    # Modify anchor sizes and aspect ratios for varying object sizes
    model.rpn.anchor_generator.sizes = ((8,), (16,), (32,), (64,), (128,)) # Wider range of anchor sizes
    model.rpn.anchor_generator.aspect_ratios = ((0.25, 0.5, 1.0, 2.0, 4.0),) * len(model.rpn.anchor_generator.sizes)

    # Increase number of proposals for better detection of varying sizes
    model.rpn.pre_nms_top_n_train = 2000  # Default is 1000
    model.rpn.post_nms_top_n_train = 1000 # Default is 500
    model.rpn.pre_nms_top_n_test = 1000   # Default is 1000
    model.rpn.post_nms_top_n_test = 500   # Default is 500

    # Replace box predictor with larger hidden layer
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    model.roi_heads.box_predictor = torchvision.models.detection.faster_rcnn.FastRCNNPredictor(
        in_features, num_classes)

    # Replace mask predictor with larger hidden layer
    in_features_mask = model.roi_heads.mask_predictor.conv5_mask.in_channels
    hidden_layer = 512  # Increased from 256
    model.roi_heads.mask_predictor = torchvision.models.detection.mask_rcnn.MaskRCNNPredictor(
        in_features_mask, hidden_layer, num_classes)

    return model

# --- Lightning Module ---
class InstanceSegmentationModule(pl.LightningModule):
    def __init__(self, model, lr=0.005):
        super().__init__()
        self.model = model
        self.lr = lr

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        images, targets = batch
        images = list(image for image in images)
        targets = [{k: v for k, v in t.items()} for t in targets]
        
        loss_dict = self.model(images, targets)
        losses = sum(loss for loss in loss_dict.values())
        
        # Log losses
        for k, v in loss_dict.items():
            self.log(f'train_{k}', v, on_step=True, on_epoch=True)
        self.log('train_loss', losses, on_step=True, on_epoch=True)
        
        return losses

    def validation_step(self, batch, batch_idx):
        images, targets = batch
        images = list(image for image in images)
        targets = [{k: v for k, v in t.items()} for t in targets]
        
        loss_dict = self.model(images, targets)
        losses = sum(loss for loss in loss_dict.values())
        
        # Log losses
        for k, v in loss_dict.items():
            self.log(f'val_{k}', v, on_epoch=True)
        self.log('val_loss', losses, on_epoch=True)
        
        # Log example predictions at end of epoch
        if batch_idx == 0:
            self.log_predictions(images[0], targets[0])
            
        return losses
    
    def log_predictions(self, image, target):
        # Get model predictions
        self.model.eval()
        with torch.no_grad():
            prediction = self.model([image])[0]
        
        # Convert tensors to numpy
        image_np = image.cpu().numpy().transpose(1,2,0)
        
        fig, ax = plt.subplots(figsize=(8,8))
        ax.imshow(image_np, cmap='gray')
        
        # Plot ground truth boxes and masks
        boxes = target['boxes'].cpu().numpy()
        masks = target['masks'].cpu().numpy()
        
        for box, mask in zip(boxes, masks):
            x1, y1, x2, y2 = box
            rect = plt.Rectangle((x1,y1), x2-x1, y2-y1, 
                               fill=False, color='green', linewidth=1)
            ax.add_patch(rect)
            ax.imshow(np.ma.masked_where(mask==0, mask),
                     alpha=0.3, cmap='Greens')
            
        # Plot predictions
        pred_boxes = prediction['boxes'].cpu().numpy()
        pred_masks = prediction['masks'].cpu().numpy()
        scores = prediction['scores'].cpu().numpy()
        
        for box, mask, score in zip(pred_boxes, pred_masks, scores):
            if score > 0.5:
                x1, y1, x2, y2 = box
                rect = plt.Rectangle((x1,y1), x2-x1, y2-y1,
                                   fill=False, color='red', linewidth=1)
                ax.add_patch(rect)
                ax.imshow(np.ma.masked_where(mask[0]==0, mask[0]),
                         alpha=0.3, cmap='Reds')
                
        plt.axis('off')
        self.logger.experiment.add_figure('predictions', fig, self.current_epoch)
        plt.close()

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.model.parameters(), lr=self.lr)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=3, gamma=0.1)
        return {
            'optimizer': optimizer,
            'lr_scheduler': scheduler
        }

# --- Training Functions ---
def split_dataset(dataset_root, val_ratio=0.2, seed=42):
    # Create train dataset with augmentations
    train_dataset = InstanceSegmentationDataset(dataset_root, train=True)
    
    # Create validation dataset with validation transforms
    val_dataset = InstanceSegmentationDataset(dataset_root, train=False)
    
    # Split indices
    total_size = len(train_dataset)
    val_size = int(total_size * val_ratio)
    train_size = total_size - val_size
    
    # Split indices using random_split
    train_indices, val_indices = random_split(
        range(total_size),
        [train_size, val_size],
        generator=torch.Generator().manual_seed(seed)
    )
    
    # Create subset datasets
    train_dataset = Subset(train_dataset, train_indices)
    val_dataset = Subset(val_dataset, val_indices)
    
    return train_dataset, val_dataset

def train_model(model, train_dataset, val_dataset, max_epochs=10):
    # Create data loaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=8,
        shuffle=True,
        num_workers=4,
        collate_fn=lambda x: tuple(zip(*x))
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=8,
        shuffle=False,
        num_workers=4,
        collate_fn=lambda x: tuple(zip(*x))
    )
    
    # Initialize Lightning module
    lightning_model = InstanceSegmentationModule(model)
    
    # Setup logger
    logger = TensorBoardLogger('lightning_logs', name='instance_segmentation')
    
    # Initialize trainer
    trainer = pl.Trainer(
        max_epochs=max_epochs,
        logger=logger,
        accelerator='auto'
    )
    
    # Train model
    trainer.fit(
        lightning_model,
        train_loader,
        val_loader
    )
    
    return lightning_model


def generate_test_data(num_samples=100, image_size=(256, 256)):
    # Create random test dataset
    test_dataset = []
    
    for _ in range(num_samples):
        # Generate random image
        image = torch.randint(0, 65535, (1, *image_size), dtype=torch.float32)
        
        # Generate random number of objects (1-5)
        num_objects = torch.randint(1, 6, (1,)).item()
        
        # Generate random boxes and masks
        boxes = []
        masks = []
        labels = []
        
        for _ in range(num_objects):
            # Random box dimensions
            x1 = torch.randint(0, image_size[0]-60, (1,)).item()
            y1 = torch.randint(0, image_size[1]-60, (1,)).item()
            x2 = x1 + torch.randint(30, 60, (1,)).item()
            y2 = y1 + torch.randint(30, 60, (1,)).item()
            
            boxes.append([x1, y1, x2, y2])
            
            # Create random mask
            mask = torch.zeros(image_size)
            mask[y1:y2, x1:x2] = 1
            masks.append(mask)
            
            # Add label (1 for object)
            labels.append(1)
            
        # Convert to tensors
        boxes = torch.tensor(boxes, dtype=torch.float32)
        masks = torch.stack(masks)
        labels = torch.tensor(labels, dtype=torch.int64)
        
        # Create target dictionary
        target = {
            'boxes': boxes,
            'labels': labels,
            'masks': masks
        }
        
        test_dataset.append((image, target))
    
    return test_dataset



test_dataset = generate_test_data()

# Create model
model = get_model_instance_segmentation(in_chans=1, num_classes=2)  # Use 91 classes as expected by model

# Split dataset
train_data, val_data = split_dataset(test_dataset)

# Train model
trained_model = train_model(
    model=model,
    train_dataset=train_data,
    val_dataset=val_data,
    max_epochs=5
)



if __name__ == "__main__":
    # Generate test dataset
    test_dataset = generate_test_data()
    
    # Create model
    model = get_model_instance_segmentation(in_chans=1, num_classes=2)  # Use 91 classes as expected by model
    
    # Split dataset
    train_data, val_data = split_dataset(test_dataset)
    
    # Train model
    trained_model = train_model(
        model=model,
        train_dataset=train_data,
        val_dataset=val_data,
        max_epochs=5
    )

