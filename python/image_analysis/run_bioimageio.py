#!/usr/bin/env python3
"""
EM Image Segmentation using bioimage.io predict_many

This script segments electron microscopy images using bioimage.io models.
It processes multiple images and saves the segmentation results.
"""
#%%
import os
import glob
from pathlib import Path
import numpy as np
from typing import List, Optional, Union
import logging
import argparse

try:
    from bioimageio.core import load_resource
    from bioimageio.core.prediction import predict_many
    import imageio.v3 as imageio
    import tifffile
except ImportError as e:
    print(f"Missing dependencies: {e}")
    print("Install with: pip install bioimageio.core imageio tifffile")
    exit(1)

def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Segment EM images using bioimage.io models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default cellpose model
  python em_segmentation.py -i /path/to/images/ -o results/
  
  # Use specific model
  python em_segmentation.py -i image.tif -o results/ -m "deepimagej/cellpose_nuclei_model"
  
  # Process single image with custom threshold
  python em_segmentation.py -i single_image.tif -o results/ -t 0.7
  
  # Batch process with verbose output
  python em_segmentation.py -i /path/to/images/ -o results/ -v
  
Popular models:
  - deepimagej/cellpose_model_bacteria_omnipose (bacteria)
  - ilastik/cellpose_model (general cells)
  - deepimagej/cellpose_nuclei_model (nuclei)
  
Find more models at: https://bioimage.io/
        """
    )
    
    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Input image file or directory containing images"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output directory for segmentation results"
    )
    
    parser.add_argument(
        "-m", "--model",
        type=str,
        default="deepimagej/cellpose_model_bacteria_omnipose",
        help="Bioimage.io model name (default: deepimagej/cellpose_model_bacteria_omnipose)"
    )
    
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.5,
        help="Probability threshold for binary segmentation (default: 0.5)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    parser.add_argument(
        "--list-formats",
        action="store_true",
        help="List supported image formats and exit"
    )
    
    parser.add_argument(
        "--prefix",
        type=str,
        default="segmented",
        help="Prefix for output filenames (default: segmented)"
    )
    
    return parser.parse_args()

def list_supported_formats():
    """List supported image formats"""
    formats = [
        "TIFF (.tif, .tiff)",
        "PNG (.png)",
        "JPEG (.jpg, .jpeg)",
        "BMP (.bmp)",
        "GIF (.gif)"
    ]
    
    print("Supported image formats:")
    for fmt in formats:
        print(f"  - {fmt}")
    print("\nNote: TIFF format is recommended for scientific imaging data.")

def setup_logging(verbose: bool):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s' if verbose else '%(levelname)s: %(message)s'
    
    logging.basicConfig(
        level=level,
        format=format_str,
        datefmt='%Y-%m-%d %H:%M:%S'
    )

# Initialize logger
logger = logging.getLogger(__name__)

class EMSegmenter:
    """
    A class to handle EM image segmentation using bioimage.io models
    """
    
    def __init__(self, model_name: str, threshold: float = 0.5):
        """
        Initialize the segmenter with a bioimage.io model
        
        Args:
            model_name: Name/ID of the bioimage.io model (will be downloaded automatically)
            threshold: Probability threshold for binary segmentation
        """
        self.model_name = model_name
        self.threshold = threshold
        self.model = None
        self._load_model()
    
    def _load_model(self):
        """Load the bioimage.io model (downloads automatically if needed)"""
        try:
            logger.info(f"Loading/downloading model: {self.model_name}")
            
            # Load resource description (this will download the model if needed)
            self.model = load_resource(self.model_name)
            
            logger.info("Model loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load model '{self.model_name}': {e}")
            logger.info("Available models can be found at: https://bioimage.io/")
            raise
    
    def load_images(self, input_path: Union[str, List[str]]) -> List[np.ndarray]:
        """
        Load images from file paths
        
        Args:
            input_path: Single file path, directory path, or list of file paths
            
        Returns:
            List of loaded images as numpy arrays
        """
        images = []
        
        if isinstance(input_path, str):
            if os.path.isdir(input_path):
                # Load all images from directory
                image_extensions = ['*.tif', '*.tiff', '*.png', '*.jpg', '*.jpeg']
                file_paths = []
                for ext in image_extensions:
                    file_paths.extend(glob.glob(os.path.join(input_path, ext)))
                    file_paths.extend(glob.glob(os.path.join(input_path, ext.upper())))
            else:
                # Single file
                file_paths = [input_path]
        else:
            # List of file paths
            file_paths = input_path
        
        logger.info(f"Found {len(file_paths)} images to process")
        
        for file_path in sorted(file_paths):
            try:
                if file_path.lower().endswith(('.tif', '.tiff')):
                    image = tifffile.imread(file_path)
                else:
                    image = imageio.imread(file_path)
                
                # Ensure image is 2D for EM data
                if image.ndim > 2:
                    if image.shape[-1] == 3:  # RGB image
                        image = np.mean(image, axis=-1)  # Convert to grayscale
                    else:
                        image = np.squeeze(image)
                
                images.append(image.astype(np.float32))
                logger.info(f"Loaded image: {file_path} - Shape: {image.shape}")
                
            except Exception as e:
                logger.error(f"Failed to load image {file_path}: {e}")
                continue
        
        return images
    
    def preprocess_images(self, images: List[np.ndarray]) -> List[np.ndarray]:
        """
        Preprocess images for the model
        
        Args:
            images: List of input images
            
        Returns:
            List of preprocessed images
        """
        preprocessed = []
        
        for img in images:
            # Normalize to 0-1 range
            img_norm = (img - img.min()) / (img.max() - img.min())
            
            # Add batch and channel dimensions if needed
            if img_norm.ndim == 2:
                img_norm = img_norm[np.newaxis, np.newaxis, ...]  # (1, 1, H, W)
            elif img_norm.ndim == 3:
                img_norm = img_norm[np.newaxis, ...]  # (1, C, H, W)
            
            preprocessed.append(img_norm)
        
        return preprocessed
    
    def segment_images(self, images: List[np.ndarray]) -> List[np.ndarray]:
        """
        Segment images using predict_many
        
        Args:
            images: List of preprocessed images
            
        Returns:
            List of segmentation results
        """
        try:
            logger.info(f"Starting segmentation of {len(images)} images...")
            
            # Use predict_many for batch processing
            predictions = predict_many(
                model=self.model,
                inputs=images
            )
            
            logger.info("Segmentation completed successfully")
            return predictions
            
        except Exception as e:
            logger.error(f"Segmentation failed: {e}")
            raise
    
    def postprocess_results(self, predictions: List[np.ndarray]) -> List[np.ndarray]:
        """
        Postprocess segmentation results
        
        Args:
            predictions: Raw model predictions
            
        Returns:
            List of processed segmentation masks
        """
        processed = []
        
        for pred in predictions:
            # Remove batch dimension and get the segmentation mask
            if pred.ndim > 2:
                mask = np.squeeze(pred)
            else:
                mask = pred
            
            # Convert probabilities to binary mask if needed
            if mask.max() <= 1.0:  # Probability map
                mask = (mask > self.threshold).astype(np.uint8)
            else:  # Already segmentation labels
                mask = mask.astype(np.uint8)
            
            processed.append(mask)
        
        return processed
    
    def save_results(self, masks: List[np.ndarray], output_dir: str, 
                    input_paths: Optional[List[str]] = None, prefix: str = "segmented"):
        """
        Save segmentation results
        
        Args:
            masks: List of segmentation masks
            output_dir: Directory to save results
            input_paths: Original input file paths for naming
            prefix: Prefix for output filenames
        """
        os.makedirs(output_dir, exist_ok=True)
        
        for i, mask in enumerate(masks):
            if input_paths and i < len(input_paths):
                # Use original filename
                base_name = Path(input_paths[i]).stem
                output_path = os.path.join(output_dir, f"{base_name}_{prefix}.tif")
            else:
                output_path = os.path.join(output_dir, f"{prefix}_{i:03d}.tif")
            
            tifffile.imwrite(output_path, mask)
            logger.info(f"Saved segmentation result: {output_path}")

def main():
    """
    Main function to run EM image segmentation
    """
    # Parse command line arguments
    args = parse_arguments()
    
    # Handle special commands
    if args.list_formats:
        list_supported_formats()
        return
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    # Validate inputs
    if not os.path.exists(args.input):
        logger.error(f"Input path does not exist: {args.input}")
        return
    
    logger.info(f"Starting EM image segmentation")
    logger.info(f"Model: {args.model}")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Threshold: {args.threshold}")
    
    try:
        # Initialize segmenter (model will be downloaded automatically)
        logger.info("Initializing segmenter...")
        segmenter = EMSegmenter(args.model, args.threshold)
        
        # Load images
        logger.info("Loading images...")
        images = segmenter.load_images(args.input)
        if not images:
            logger.error("No images loaded. Check input path and supported formats.")
            return
        
        # Preprocess images
        logger.info("Preprocessing images...")
        preprocessed_images = segmenter.preprocess_images(images)
        
        # Segment images
        logger.info("Running segmentation...")
        predictions = segmenter.segment_images(preprocessed_images)
        
        # Postprocess results
        logger.info("Postprocessing results...")
        segmentation_masks = segmenter.postprocess_results(predictions)
        
        # Prepare input paths for naming
        input_paths = None
        if isinstance(args.input, str):
            if os.path.isdir(args.input):
                # Directory input
                image_extensions = ['*.tif', '*.tiff', '*.png', '*.jpg', '*.jpeg', '*.bmp', '*.gif']
                input_paths = []
                for ext in image_extensions:
                    input_paths.extend(glob.glob(os.path.join(args.input, ext)))
                    input_paths.extend(glob.glob(os.path.join(args.input, ext.upper())))
                input_paths = sorted(input_paths)
            elif os.path.isfile(args.input):
                # Single file input
                input_paths = [args.input]
        
        # Save results
        logger.info("Saving results...")
        segmenter.save_results(segmentation_masks, args.output, input_paths, args.prefix)
        
        logger.info("EM image segmentation completed successfully!")
        logger.info(f"Results saved in: {args.output}")
        logger.info(f"Processed {len(segmentation_masks)} images")
        
    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
    except Exception as e:
        logger.error(f"Segmentation pipeline failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()

#%%
if __name__ == "__main__":
    main()
# %%
