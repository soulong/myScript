/**
 * ImageJ Macro Language (.ijm) Script for Batch IHC Staining Image Analysis
 * 1. Asks the user to select an input folder containing TIF images.
 * 2. Processes each TIF file:
 * a. Opens the image.
 * b. Converts it to RGB if necessary.
 * c. Performs Colour Deconvolution (using H DAB vector).
 * d. Measures the mean intensity for the three resulting channels (DAB, Hema, Residual).
 * e. Closes temporary images.
 * 3. Saves the results to a CSV file.
 * Note: Assumes a typical Hematoxylin/DAB staining. Adjust the 'vector' variable if needed.
 */

// Define the name of the results file
resultsFileName = "IHC_COL1.csv";

// Define the Colour Deconvolution vector for Hematoxylin (H) and DAB (D)
vector = "H PAS"; 

// Speed up processing by hiding image windows
setBatchMode(false); 


// --- Main Script Logic ---

// 1. Get the input folder from the user
inputDir = getDirectory("Select Input Folder Containing TIF Images");
if (inputDir == "") exit();

// 2. Define the output file path (e.g., in the parent directory of the input folder)
parentDir = File.getParent(inputDir);
resultsPath = parentDir + File.separator + resultsFileName;

// 3. Clear existing results and set measurements
run("Clear Results");
// Measure Mean Gray Value only
run("Set Measurements...", "mean decimal=3");

// Open and prepare the results file
f = File.open(resultsPath);
print(f, "Filename,Mean_C1,Mean_C2,Mean_C3");

// 4. Get the list of TIF files
list = getFileList(inputDir);

// Check if any image files were found
if (list.length == 0) {
    showMessage("No image files found in the selected directory!");
    exit();
}


// 5. Process each file
for (i = 0; i < list.length; i++) {
    fileName = list[i];
    print(i + ": " + fileName);
    
    // Only process TIF files
    if (endsWith(fileName, ".tif") || endsWith(fileName, ".TIF")) {
        filePath = inputDir + fileName;

        // Open the image
        open(filePath);
        title = getTitle();
        getDimensions(width, height, channels, slices, frames);
        //print(channels);
        // Ensure image is RGB (Colour Deconvolution requires 3 channels)
        /*if (channels != 3) {
            print("Warning: Image " + fileName + " is not RGB. Converting to RGB.");
            run("RGB Color");
        }*/
        
        // Perform Colour Deconvolution
        // This generates three new images: C1, C2, C3 (DAB, Hema, Residual)
        run("Colour Deconvolution", "vectors=[" + vector + "] hide");
        /*openedimagelist = getList("image.titles");
        for (i = 0; i < openedimagelist.length; i++) {
        	print("File " + (i + 1) + ": " + openedimagelist[i]);
        }*/
        
        // Get the resultant images' names (must match the original title)
        C1_title = title + "-(Colour_1)";
        C2_title = title + "-(Colour_2)";
        C3_title = title + "-(Colour_3)";
        
        // Initialize intensity variables
        meanC1 = 0;
        meanC2 = 0;
        meanC3 = 0;

        // --- Measure Mean Intensity for DAB channel (C1) ---
        selectWindow(C1_title);
        run("Measure");
        meanC1 = getResult("Mean", nResults-1);
        close();

        // --- Measure Mean Intensity for Hema channel (C2) ---
        selectWindow(C2_title);
        run("Measure");
        meanC2 = getResult("Mean", nResults-1);
        close();

        // --- Measure Mean Intensity for Residual channel (C3) ---
        selectWindow(C3_title);
        run("Measure");
        meanC3 = getResult("Mean", nResults-1);
        close();

        // Close the original image
        selectWindow(title);
        close();

        // Write results to the CSV file
        print(f, fileName + "," + meanC1 + "," + meanC2 + "," + meanC3);
    }
}

// 6. Finalize and report
File.close(f);
print("--- Analysis Complete ---");
print("Results saved to: " + resultsPath);
run("Clear Results");