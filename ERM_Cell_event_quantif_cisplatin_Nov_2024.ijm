run("Set Scale...", "distance=3.1008 known=1 unit=µm global");
filename = getTitle() ;

run("Split Channels");

// In the blue channel, nuclei area measurment
selectWindow(filename + " (blue)");

run("Duplicate...", " ");

copyDAPI = getTitle();
rename(copyDAPI + "TotNucArea");

//threshold blue channel to retrieve nuclei area
setAutoThreshold("Default dark");
run("Threshold...");
run("Smooth");

//Choose appropriate threshold
//STS117
//N1[30-255]
setThreshold(30, 255);
setOption("BlackBackground", true);
run("Convert to Mask");

// to select the spheroid
setTool("polygon");
waitForUser("Waiting for user to draw a polygon. Press Okay to continue....");


// to remove all the unwanted tissue
run("Make Inverse");
setBackgroundColor(0, 0, 0);
run("Clear", "slice");
run("Create Selection");
run("Measure");

run("Duplicate...", " ignore");
copyDAPI_CAIX = getTitle();
rename(copyDAPI_CAIX + "TotNucAreaCAIXpos");

// In the red channel, CAIX positive nuclei area measurment

selectWindow(filename + " (red)");

run("Duplicate...", " ");
copyCAIX = getTitle();
rename(copyCAIX + "CAIXPosArea");

selectWindow(copyCAIX + "CAIXPosArea");

run("Subtract Background...", "rolling=50 sliding");
setAutoThreshold("Default dark");
run("Threshold...");
setTool("polygon");
waitForUser("Waiting for user to draw a polygon. Press Okay to continue....");


// to retrieve CAIX positive nuclei area
selectWindow(copyCAIX + "CAIXPosArea");
selectWindow(copyDAPI_CAIX + "TotNucAreaCAIXpos");

run("Restore Selection");
run("Make Inverse");
run("Clear", "slice");
run("Create Selection");
run("Measure");


// In the green channel positive CellEvent 3/7 Tot nuclei area

selectWindow(filename + " (green)");

run("Duplicate...", " ");
copyyH2AX = getTitle();
rename(copyyH2AX + "ActivatedCaspase3/7")

selectWindow(copyyH2AX + "ActivatedCaspase3/7");

run("Subtract Background...", "rolling=50 sliding");

/*selectWindow(copyDAPI + "TotNucArea");
selectWindow(copyyH2AX + "ActivatedCaspase3/7");

run("Restore Selection");
run("Make Inverse");
run("Clear", "slice");
run("Select None");

setAutoThreshold("Default dark");
run("Threshold...");
run("Smooth");

//Choose appropriate threshold
//STS117
//N1[30-255]
setThreshold(25, 255);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Create Selection");
run("Measure");

// In the green channel positive CellEvent 3/7 CAIX pos nuclei area
selectWindow(copyyH2AX + "ActivatedCaspase3/7");
run("Duplicate...", " ");
copyCAIX = getTitle();
rename(copyCAIX + "ActivatedCaspase3/7_CAIXPosArea");

selectWindow(copyDAPI_CAIX + "TotNucAreaCAIXpos");
selectWindow(copyCAIX + "ActivatedCaspase3/7_CAIXPosArea");
run("Restore Selection");
run("Make Inverse");
run("Clear", "slice");
run("Select None");
run("Create Selection");
run("Measure");

// save CAIX pos area

selectWindow(copyCAIX + "CAIXPosArea");
run("Save");

// close everything

//close("*");
