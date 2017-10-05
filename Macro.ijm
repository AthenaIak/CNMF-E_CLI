//path = "/home/athina/Data/";
//nam = "recording_20170724_102219";
//nam = 'test_TRP';

// Get arguments from command line (in the format path+recordingID+#files)
input = getArgument();
options = split(input, "+");
path = options[0]+"/";
nam = options[1];
numFiles = parseInt(options[2]);

if (path=="") exit ("No path argument!");
if (nam=="") exit ("No nam argument!");

// Figure out the names of the files used and created
ppMovie1 = "pp_" + "recording_" + nam + "-1.tif";
ppMovie2 = "pp_" + "recording_" + nam + "-2.tif";
mcMovie1 = "mc_" + "recording_" + nam + "-1.tif";
mcMovie2 = "mc_" + "recording_" + nam + "-2.tif";
refImg = "pp_" + "recording_" + nam + "-ref.tif";

width = 1440; height = 1080;


// Create the reference image
open(path+ppMovie1);
run("Duplicate...", "title="+refImg);
saveAs("Tiff",path+refImg);

// motion correct the first movie
turboreg(ppMovie1, refImg, path+mcMovie1);
selectWindow(mcMovie1);
//saveAs("Tiff", path+mcMovie1); //not needed because java does that now

// temporarily calculate min intensity (to figure out boundaries created by motion correction)
run("Z Project...", "projection=[Min Intensity]");
selectWindow(ppMovie1);
close();
selectWindow(mcMovie1);
close();

// check that a second file exists for the same movie
if (numFiles == 2){
	open(path+ppMovie2);

	// motion correct the second movie
	turboreg(ppMovie2, refImg, path+mcMovie2);
	selectWindow(ppMovie2);
	close();
	selectWindow(refImg);
	close();

	selectWindow(mcMovie2);
	run("Z Project...", "projection=[Min Intensity]");
	selectWindow(mcMovie2);
	close();

	run("Concatenate...", "  title=[Concatenated Stacks] image1=MIN_"+mcMovie1+" image2=MIN_"+mcMovie2+" image3=[-- None --]");
	run("Z Project...", "projection=[Min Intensity]");
}
else {
	selectWindow("MIN_"+mcMovie1);
}


saveAs("Tiff",path+"mc_minIntensity_"+nam+".tif");
close();

if (numFiles == 2){
	selectWindow("Concatenated Stacks");
	close();
}

run("Quit");


function turboreg(source, reference, outputFilename)
{
	run("TurboReg FullCLI",
		"-align "
		+ "-window " + source + " "// Source (window reference).
		+ "0 0 " + (width-1) + " " + (height-1) + " " // No cropping.
		+ "-window " + reference + " "// Target (file reference).
		+ "0 0 " + (width-1) + " " + (height-1) + " " // No cropping.
		+ "-rigidBody " // This corresponds to rotation and translation.
		+ (width / 2) + " " + (height / 2) + " " // Source translation landmark.
		+ (width / 2) + " " + (height / 2) + " " // Target translation landmark.
		+ "0 " + (height / 2) + " " // Source first rotation landmark.
		+ "0 " + (height / 2) + " " // Target first rotation landmark.
		+ (width - 1) + " " + (height / 2) + " " // Source second rotation landmark.
		+ (width - 1) + " " + (height / 2) + " " // Target second rotation landmark.
		+ "-batch "
		+ "-hideOutput "
		+ outputFilename);
}
