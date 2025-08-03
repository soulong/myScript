// FIJI/ImageJ Script for Batch RGB Processing
// Author: Assistant
// Description: Split RGB images, convert to grayscale, set scale, save green channel

// Get input and output directories
inputDir = getDirectory("Choose input directory with RGB images");
outputDir = getDirectory("Choose output directory for processed images");
//inputDir = "C:/Users/haohe/Desktop/20250523数据处理/20250523";
//outputDir = "C:/Users/haohe/Desktop/20250523数据处理/green";

// Get list of image files in input directory
fileList = getFileList(inputDir);

// Filter for common image formats
imageFiles = newArray();
for (i = 0; i < fileList.length; i++) {
    if (endsWith(fileList[i], ".tif") || endsWith(fileList[i], ".tiff") || 
        endsWith(fileList[i], ".jpg") || endsWith(fileList[i], ".jpeg") || 
        endsWith(fileList[i], ".png") || endsWith(fileList[i], ".bmp")) {
        imageFiles = Array.concat(imageFiles, fileList[i]);
    }
}

// Check if any image files were found
if (imageFiles.length == 0) {
    showMessage("No image files found in the selected directory!");
    exit();
}

// Set scale parameters (modify these values according to your microscope setup)
pixel_distance = 30.7692;
known_distance = 1;
unit = "micron";

// Set output file parameters
outputSuffix = "_green";  // Will create files like "image_green.tif"
                         // Examples: "_G", "_green_channel", "_processed", etc.
outputFormat = "Tiff";    // Output format: "Tiff", "Jpeg", "PNG", "BMP"
outputExtension = ".tif"; // File extension: ".tif", ".jpg", ".png", ".bmp"
                         // Make sure format and extension match!

// Channel selection for multi-channel images
targetChannel = 2;        // Which channel to keep (1=first, 2=second, 3=third, etc.)
                         // For RGB: 1=red, 2=green, 3=blue
                         // For multi-channel stacks: depends on your specific setup


// Process each image file
setBatchMode(false); // Speed up processing by hiding image windows


for (i = 0; i < imageFiles.length; i++) {
    // Show progress
    showProgress(i, imageFiles.length);
    print("Processing file " + (i+1) + "/" + imageFiles.length + ": " + imageFiles[i]);
    
    // Open the image
    open(inputDir + imageFiles[i]);
    originalTitle = getTitle();
    // Get filename without extension for saving
    baseName = substring(originalTitle, 0, lastIndexOf(originalTitle, "."));
    
    // Determine image type and process accordingly
    isRGB = (bitDepth() == 24);
    isMultiChannel = (nSlices > 1) || (indexOf(originalTitle, "C1-") >= 0) || (indexOf(originalTitle, "C2-") >= 0);
    print("isRGB: " + isRGB);
    print("isMultiChannel: " + isMultiChannel);
    
    if (isRGB) {
        // Handle RGB images
        print("Processing RGB image: " + originalTitle);
        
        // 1. Split RGB channels
        run("Split Channels");
        
        // Close unwanted channels, keep only target channel
        channelNames = newArray(baseName + ".tif (red)", baseName + ".tif (green)", baseName + ".tif (blue)");
        targetChannelName = channelNames[targetChannel - 1];
        
        // Close other channels
        for (j = 0; j < channelNames.length; j++) {
            if (j != (targetChannel - 1)) {
                if (isOpen(channelNames[j])) {
                    selectWindow(channelNames[j]);
                    close();
                }
            }
        }
        
        // Select target channel
        selectWindow(targetChannelName);
        
    } else if (isMultiChannel || nSlices > 1) {
        // Handle multi-channel stack images
        print("Processing multi-channel stack: " + originalTitle);
        
        // Split channels
        run("Split Channels");
        
        // Find available channel windows
        availableChannels = newArray();
        channelTitles = newArray();
        
        // Check for C1-, C2-, C3- etc. format windows
        for (c = 1; c <= 10; c++) {  // Check up to 10 channels
            channelTitle = "C" + c + "-" + originalTitle;
            if (isOpen(channelTitle)) {
                availableChannels = Array.concat(availableChannels, c);
                channelTitles = Array.concat(channelTitles, channelTitle);
            }
        }
        
        if (availableChannels.length == 0) {
            print("Warning: No channel windows found for " + originalTitle + ". Skipping...");
            close();
            continue;
        }
        
        print("Found " + availableChannels.length + " channels: " + String.join(availableChannels, ", "));
        
        // Check if target channel exists
        targetExists = false;
        targetWindowTitle = "";
        for (j = 0; j < availableChannels.length; j++) {
            if (availableChannels[j] == targetChannel) {
                targetExists = true;
                targetWindowTitle = channelTitles[j];
                break;
            }
        }
        
        if (!targetExists) {
            print("Warning: Target channel " + targetChannel + " not found in " + originalTitle + ". Using channel " + availableChannels[0] + " instead.");
            targetChannel = availableChannels[0];
            targetWindowTitle = channelTitles[0];
        }
        
        // Close unwanted channels
        for (j = 0; j < channelTitles.length; j++) {
            if (channelTitles[j] != targetWindowTitle) {
                selectWindow(channelTitles[j]);
                close();
            }
        }
        
        // Select target channel
        selectWindow(targetWindowTitle);
        
    } else {
        // Handle single channel images
        print("Processing single channel image: " + originalTitle);
        // Image is already selected, no splitting needed
    }
    
    
    // 2. Convert to 8-bit grayscale (if not already)
    if (bitDepth() != 8) {
        run("8-bit");
    }
    
    // 3. Set scale to microns
    run("Set Scale...", "distance=" + pixel_distance + " known=" + known_distance + " unit=" + unit);
    
    // 4. Save green channel to separate file
    saveAs(outputFormat, outputDir + baseName + outputSuffix + outputExtension);
    
    // Close the green channel
    close();
    
    // Clean up any remaining windows
    while (nImages > 0) {
        close();
    }
}


setBatchMode(false); // Re-enable image display

// Final message
print("\\Update:Processing complete!");
print("Processed " + imageFiles.length + " images.");
print("Green channel images saved to: " + outputDir);
showMessage("Batch processing complete!\n\nProcessed " + imageFiles.length + " images.\nGreen channel images saved to output directory.");

// Optional: Clear the log window
//run("Clear Results");