//run("Set Scale...", "distance=0.2949 known=1 unit=µm global"); // N1
run("Set Scale...", "distance=0.3550 known=1 unit=µm global"); // N2; N3

// Core
filename = getTitle() ;
selectWindow(filename);

run("8-bit");

// To crop scale and time point
//makeRectangle(0, 0, 1280, 860); // N : 1
//makeRectangle(0, 0, 1536, 1048); // N : 2
//run("Crop");

run("Duplicate...", " ");

copyCore = getTitle();
rename(copyCore + "_Core");
///////////////////////////////////////////////

// Intermediate
selectWindow(filename);

run("Duplicate...", " ");

copyIntermediate = getTitle();
rename(copyIntermediate + "_Intermediate");
//////////////////////////////////////////////

// Edge

selectWindow(filename);

run("Duplicate...", " ");

copyEdge = getTitle();
rename(copyEdge + "_Edge");
//////////////////////////////////////////////

// Core Thresholding - N1 : 0-55; N2 : 0-65; N3 : 0-65; N3-E8 day 0 to day 3 not the same brighness, 
// Difference of 40 in pixel value for the same region
// N3-E8 : 105
// N4 : 0-70
selectWindow(copyCore + "_Core");
setAutoThreshold("Default dark");
run("Threshold...");
//setThreshold(0, 55);
setThreshold(0, 65);
//setThreshold(0, 70);
//setThreshold(0, 105);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");

//Analyzing particle
run("Analyze Particles...","size=15000-900000 circularity=0.00-1.00 show=Outlines display include");

///////////////////////////////////////////

// Intermediate Thresholding - N1 : 0-75; N2 : 0-75; N3: 0-75; N4: 0-80
selectWindow(copyIntermediate + "_Intermediate");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 75);
//setThreshold(0,115);
//setThreshold(0, 80);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");

//Analyzing particle
run("Analyze Particles...","size=30000-900000 circularity=0.00-1.00 show=Outlines display include");

///////////////////////////////////////////

// Edge Thresholding - N1 : 0-150; N2 : 0-150; N3: 0-150; N4: 0-150
selectWindow(copyEdge + "_Edge");
setAutoThreshold("Default dark");
run("Threshold...");
//setThreshold(0, 140);
setThreshold(0, 150);
//setThreshold(0, 195);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");


//Analyzing particle
run("Analyze Particles...","size=90000-5900000 circularity=0.00-1.00 show=Outlines display include");

//selectWindow(copyCore + "_Core");
//close();
//selectWindow(copyIntermediate + "_Intermediate");
//close();
//selectWindow(copyEdge + "_Edge");
//close();
//selectWindow(filename);
//close();
