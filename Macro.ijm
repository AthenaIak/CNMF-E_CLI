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

// Figure out the name of the first file used and created
ppMovieCurr = "pp_" + "recording_" + nam + "-1.tif";
mcMovieCurr = "mc_recording_"+nam+"-1.tif"
refImg = "pp_" + "recording_" + nam + "-ref.tif";

width = 1440; height = 1080;

// Create the reference image
open(path+ppMovieCurr);
run("Duplicate...", "title="+refImg);
saveAs("Tiff",path+refImg);

// motion correct the first movie
turboreg(ppMovieCurr, refImg, path+mcMovieCurr);
//mc movie is automatically saved by plugin
selectWindow(mcMovieCurr);

// temporarily calculate min intensity (to figure out boundaries created by motion correction)
run("Z Project...", "projection=[Min Intensity]");
rename("MIN_" + "mc_" + "recording_" + nam + ".tif");
selectWindow(ppMovieCurr);
close();
selectWindow(mcMovieCurr);
close();

for (i=2; i-1<numFiles; i++) {
	//figure out the next file names (and remember the previous file)
	ppMovieCurr = "pp_" + "recording_" + nam + "-" + i + ".tif";
	mcMovieCurr = "mc_" + "recording_" + nam + "-" + i + ".tif";
	
	open(path+ppMovieCurr);
	
	turboreg(ppMovieCurr, refImg, path+mcMovieCurr);
	//mc movie is automatically saved by plugin
	selectWindow(mcMovieCurr);
	
	// temporarily calculate min intensity (to figure out boundaries created by motion correction)
	run("Z Project...", "projection=[Min Intensity]");
	selectWindow(ppMovieCurr);
	close();
	selectWindow(mcMovieCurr);
	close();
	
	run("Concatenate...", "  title=[Concatenated Stacks] image1=MIN_mc_recording_" + nam + ".tif" + " image2=MIN_" + mcMovieCurr + " image3=[-- None --]");
	selectWindow("Concatenated Stacks");
	run("Z Project...", "projection=[Min Intensity]");
	
	selectWindow("MIN_Concatenated Stacks");
	rename("MIN_" + "mc_" + "recording_" + nam + ".tif");
	selectWindow("Concatenated Stacks");
	close();
}

selectWindow("MIN_" + "mc_" + "recording_" + nam + ".tif");
saveAs("Tiff",path+"mc_MINintensity_"+nam+".tif");
close();
selectWindow(refImg);
close();

function turboreg(source, reference, outputFilename)
{
	run("TurboReg FullCLI",
		"-align "
		+ "-window " + source + " "// Source (window reference).
		+ "0 0 " + (width) + " " + (height) + " " // No cropping.
		+ "-window " + reference + " "// Target (file reference).
		+ "0 0 " + (width) + " " + (height) + " " // No cropping.
		+ "-rigidBody " // This corresponds to rotation and translation.
		+ (width / 2) + " " + (height / 2) + " " // Source translation landmark.
		+ (width / 2) + " " + (height / 2) + " " // Target translation landmark.
		+ "0 " + (height / 2) + " " // Source first rotation landmark.
		+ "0 " + (height / 2) + " " // Target first rotation landmark.
		+ (width) + " " + (height / 2) + " " // Source second rotation landmark.
		+ (width) + " " + (height / 2) + " " // Target second rotation landmark.
		+ "-batch "
		+ "-hideOutput "
		+ outputFilename);
}
