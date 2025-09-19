// Fiji macro to calculate the area in each channel and the overlap of green in red ROI
// REMEMBER TO ADD MEDIAN TO MEASURE BEFORE STARTING
//Load the RGB image


run("Set Scale...", "distance=3.1008 known=1 unit=µm global");
filename = getTitle() ;


// channel splitting and only keeping the red channel
run("Split Channels");
close (filename + " (blue)");
close (filename + " (green)");

//background substraction to normalize image intensity

selectImage(filename + " (red)");
run("Subtract Background...", "rolling=50");

//selecting the spheroid to remove unwanted tissue
setTool("freehand");
waitForUser("Waiting for user to circle the spheroid. Press Okay to continue....");

run("Make Inverse");
setBackgroundColor(0, 0, 0);
run("Clear", "slice");

// measure of MFI for the whole spheroid
run("Make Inverse");
run("Measure");

// to comment when working with small spheroid 



run("Duplicate...", "ignore");
run("Select None");




// measure of MFI of hypoxic region
setAutoThreshold("Default dark no-reset");
setTool("freehand");
waitForUser("Waiting for user to circle the hypoxic region. Press Okay to continue....");

run("Measure");
run("Clear", "slice");
run("Select All");

resetThreshold;

// generating a mask for the normoxic donut
setMinAndMax(0, 0);
run("Apply LUT");

run("Mean...", "radius=10");


setThreshold(10, 255);
setOption("BlackBackground", true);
run("Convert to Mask");

run("Create Selection");

 //measuring MFI of the normoxic donut

selectImage(filename + " (red)");
run("Restore Selection");
run("Measure"); 

