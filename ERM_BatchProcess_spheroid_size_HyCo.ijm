// Macro to measure Area, Intensity, Perimeter, and Shape of directory of images

run("Clear Results"); // clear the results table of any previous measurements

// The next line prevents ImageJ from showing the processing steps during 
// processing of a large number of images, speeding up the macro
setBatchMode(true); 

// Show the user a dialog to select a directory of images
inputDirectory = getDirectory("Choose a Directory of Images");

// Get the list of files from that directory
// NOTE: if there are non-image files in this directory, it may cause the macro to crash
fileList = getFileList(inputDirectory);

for (i = 0; i < fileList.length; i++)
{
    processImage(fileList[i]);
}

setBatchMode(false); // Now disable BatchMode since we are finished
updateResults();  // Update the results table so it shows the filenames

//SAVING THE DIAMETERS IN AN EXCEL SHEET
path=inputDirectory;
Dialog.create("Name Excel");
Dialog.addString("Name of the excel file", "Enter name");
Dialog.show();
name=Dialog.getString(); 
saveAs("results", path + name+".xls");


function processImage(imageFile)
{
    // Store the number of results before executing the commands,
    // so we can add the filename just to the new results
    prevNumResults = nResults;  
    
    open(imageFile);
    // Get the filename from the title of the image that's open for adding to the results table
    // We do this instead of using the imageFile parameter so that the
    // directory path is not included on the table
    filename = getTitle();
		    
		//FLATTENING THE TIFF FILE
		run("Flatten");
		
		
		//MAKING A COPY OF THE IMAGE TO APPLY THE CHANGES
		run("Duplicate...", " ");
		image_duplicate = getImageID();
		selectImage(image_duplicate);
		rename("duplicate.tif");

		//set scale
		run("Set Scale...", "distance=0.19 known=1 pixel=1 unit=um");
		
		//converting to gray scale
		run("8-bit");
	
		//Applying median filter
		run("Median...", "radius=" + 15);
		
		//Setting a threshold
		setAutoThreshold("Default");
		run("Threshold...");
		setThreshold(0, 45);
		
		//Applying threshold to make image binary
		run("Convert to Mask");
		
		//Dilating 8 times
		n_dilation=10;
		for (i = 0; i < n_dilation; i++) 
			{
		        run("Dilate");
			}
		
		//Analyzing particle
		run("Analyze Particles...","size=300000-600000 circularity=0.00-1.00 show=Outlines display include");

    // Now loop through each of the new results, and add the filename to the "Filename" column
    for (row = prevNumResults; row < nResults; row++)
    {
        setResult("Filename", row, filename);
    }

    close("*");  // Closes all images
}