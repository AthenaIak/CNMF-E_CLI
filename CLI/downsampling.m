function [] = downsampling(movieFiles, down_factor)
%% downsampling
% Takes a list of movie files, concatenates them into a movie and 
% downsamples by the specified factor. Saves the resulting movie in the 
% same directory, to one or more files of a 4 GB maximum size.
% input
%   movieFiles  :   A cell containing the filenames (path included)
%   down_factor :   The downsampling factor.

disp('--- PREPROCESSING ---');

TIFF_MAX_FRAMES = 1349;
%movieFiles={'/home/athina/Data/recording_20170711_131010.tif'};preprocessing(movieFiles,4);
[inDir,rec_nam,~] = fileparts(movieFiles{1});
run_setup;
clear CNMF_dir;

% load all files to Y
numFiles = length(movieFiles); curr_id = 1;
for f = 1:numFiles
	%load movie file
	fprintf('Loading file %d/%d.\n', f, numFiles); tic;
    tmpY = bigread2( movieFiles{f});
    fprintf('Time cost of loading images: %.2f seconds\n', toc);
    
    % concatenate all to Y
	numFrames = size(tmpY,3);
    Y(:,:,curr_id:curr_id+numFrames-1) = tmpY;
    clear tmpY;
    curr_id = curr_id + numFrames;
end
clear curr_id in_nam;

%temporally downsample
downY = Y(:,:,1:down_factor:end);
clear Y;

%save as tiff
disp('Saving movie...');
curr_id=1; numFrames = size(downY,3);
for img=1:numFrames
	if mod(img,TIFF_MAX_FRAMES)==1
    	outputFileName = fullfile(inDir,sprintf('pp_%s-%d.tif',rec_nam,curr_id));
        curr_id = curr_id+1;
	end
    imwrite(downY(:, :, img), outputFileName, 'WriteMode', 'append', 'Compression','none');
        
	if mod(img,250)==0
    	fprintf('%d/%d frames saved\n', img, numFrames);
	end
end
    
disp('Done saving movie.');    
   

disp('--- END OF PREPROCESSING ---');

clear numMovies

