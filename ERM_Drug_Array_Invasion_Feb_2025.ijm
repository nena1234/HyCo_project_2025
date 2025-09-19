//run("Set Scale...", "distance=0.292 known=1 unit=µm global"); // if 1000µm
run("Set Scale...", "distance=0.3500 known=1 unit=µm global"); // if 800µm

// Core
filename = getTitle() ;
selectWindow(filename);

run("8-bit");

// To crop scale and time point
//makeRectangle(0, 0, 1280, 860);
makeRectangle(0, 0, 1536, 1048); 
run("Crop");

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

// Core Thresholding - N1 : 0-55; N2-3 : 0-65; N4 : 0-75
selectWindow(copyCore + "_Core");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 55);
//setThreshold(0, 75);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");

//Analyzing particle
run("Analyze Particles...","size=3000-1200000 circularity=0.00-1.00 show=Outlines display include");

///////////////////////////////////////////

// Intermediate Thresholding - N1 : 0-70; N2-N3: 0-80; N4: 0-85
selectWindow(copyIntermediate + "_Intermediate");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 70);
//setThreshold(0, 85);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");

//Analyzing particle
run("Analyze Particles...","size=30000-1200000 circularity=0.00-1.00 show=Outlines display include");

///////////////////////////////////////////

// Edge Thresholding - N1-N4 : 0-120; N2 : 0-140 N4: 0-130
selectWindow(copyEdge + "_Edge");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 120);
//setThreshold(0, 140);
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
