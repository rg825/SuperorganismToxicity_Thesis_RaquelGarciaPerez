// AUTOMATIC MEASUREMENT OF BLOBS DETECTED IN R BY AN OPTIMIZED BLOB DETECTOR

// AIM: Measure the area of all rounded structures detected by the blob detector from R and export 
// their measurements as a single CSV file.

// APPROACH: The macro processes each image individually by importing the blob coordinates
// generated in R and creating circular ROIs centred on each detected blob. Once the ROIs are created, 
// Fiji measures the selected parameters and exports all measurements into a single results table.


// 1. INITIAL SETTINGS --------------------------------------------------------------------------------------------------------------------------------

// Folder containing the colony images:

imageDir = "/Users/raquelgarciaperez/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB/Annotated_Images/";

// Folder containing the blob coordinates exported from R:

blobDir = "/Users/raquelgarciaperez/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB/Detected_ROIs/";

// Output folder for the automated measurements:

outputDir = "/Users/raquelgarciaperez/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB/";

// Output folder to store all ROI sets in ZIP format for future analysis:

roiDir = "/Users/raquelgarciaperez/Documents/VIDA ACADÉMICA/IMPERIAL/TFM/DATA/GITHUB/Detected_ROIs/";


// Measurements to record in the macro:

run("Clear Results");  // Clearing the results in case the macro was run before

run("Set Measurements...",
    "area mean standard min max centroid shape display redirect=None decimal=3");


// 2. IMPORTING ROIs ON EACH IMAGE --------------------------------------------------------------------------------------------------------------------

// Listing all images that have been processed by the optimized blob detector:
imageList = getFileList(imageDir);

for(i=0; i<lengthOf(imageList); i++){
	
	imageName = imageList[i];

    if (!endsWith(imageName, ".jpeg")) 
        continue;
        
    print("");
    print("Processing: " + imageName);

    open(imageDir + imageName);

    width = getWidth();
    height = getHeight();

    roiManager("Reset");
    
    blobNumber = 0;

    csvName = "Detected_ROIs_" +
    replace(imageName, ".jpeg", ".csv");

    csvPath = blobDir + csvName;

    if (!File.exists(csvPath)){print("CSV not found.");
        close();
        continue;}

    // Reading the CSV file:
    text = File.openAsString(csvPath);   
    
    lines = split(text, "\n");    // Splitting the file into lines

    nBlobs = lengthOf(lines) - 1;   // Extracting the number of detected blobs

    print("Blobs detections in CSV: " + nBlobs);

// Process every detected blob
    for (j = 1; j < lengthOf(lines); j++) {

         line = trim(lines[j]);

          if (line == "")
              continue;

          values = split(line, ",");

          x = parseFloat(values[0]);
          y = parseFloat(values[1]);
          sigma = parseFloat(values[2]);
          
          // The DoH detector returns the Gaussian sigma at which the blob is detected.
          // Assuming that the blobs are approximately circular, I calculate the radius of
          // the blob using sigma: radius = sigma × sqrt(2):
          radius = sigma * sqrt(2);   

           // Removing any blob coinciding with the border. All colony images have the nest
           // centred, with no possibility of cocoons coinciding with any border, so we 
           // exclude any false positive in this position: 
           if (x-radius < 0) continue;
           if (y-radius < 0) continue;
           if (x+radius > width) continue;
           if (y+radius > height) continue;

           blobNumber = blobNumber + 1;    // Assigning a unique ID to blobs kept for analysis.
    
           makeOval(x-radius,
                    y-radius,
                    radius*2,
                    radius*2);

           roiManager("Add");
           roiManager("Select", roiManager("count")-1);
           roiManager("Rename", "Blob_" + blobNumber);     // Renaming ROIs using the blob numbers.
           
           run("Measure");   // Measuring parameters of interest.
           
           // We track the row where the measurements have been stored and we add metadata 
           // to the latest measurement so each ROI can be linked back to its image and 
           // the corresponding detected blob.
           
           row = nResults - 1;   
           setResult("Image", row, imageName);
           setResult("Blob_number", row, blobNumber);
           setResult("DoH_Sigma", row, sigma);
           setResult("ROI_radius", row, radius);
           updateResults();
}
    
    //Saving measurements:
    roiZip = roiDir + 
    "Detected_ROIs_" +
    replace(imageName, ".jpeg", ".zip");
         
    roiManager("Save", roiZip);
    
print(imageName + ": " + blobNumber + " blobs measured.");

close(); // Closing image before moving on to the next one.

}

// Saving the complete results table as a CSV file.
saveAs("Results", outputDir + "automated_measurements.csv");


