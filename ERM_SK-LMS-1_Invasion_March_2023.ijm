//run("Set Scale...", "distance=0.2949 known=1 unit=µm "); // if 1000 µm scale
run("Set Scale...", "distance=0.3550 known=1 unit=µm "); // if 800 µm scale

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

// Core Thresholding - N1, N2 : 0-45;
selectWindow(copyCore + "_Core");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 45);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");

//Analyzing particle
run("Analyze Particles...","size=3000-1200000 circularity=0.00-1.00 show=Outlines display include");

///////////////////////////////////////////

// Intermediate Thresholding - N1, N2 : 0-60;
selectWindow(copyIntermediate + "_Intermediate");
setAutoThreshold("Default dark");
run("Threshold...");

setThreshold(0, 60);
run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");

//Analyzing particle
run("Analyze Particles...","size=30000-1200000 circularity=0.00-1.00 show=Outlines display include");

///////////////////////////////////////////

// Edge Thresholding - N1, N2 : 125
selectWindow(copyEdge + "_Edge");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 125);

run("Smooth");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Fill Holes");


//Analyzing particle
run("Analyze Particles...","size=90000-5900000 circularity=0.00-1.00 show=Outlines display include");

/*selectWindow(copyCore + "_Core");
close();
selectWindow(copyIntermediate + "_Intermediate");
close();
selectWindow(copyEdge + "_Edge");
close();
selectWindow(filename);
close();*/
run("Close All");