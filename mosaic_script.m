clc; clear; close all;

srcDir = '~/tp/athina/Mosaic-MATLAB/src';

set_parameters='~/Data/32364/parameters_an003_mosaic';
run (set_parameters);

%parameters that lead to extensive motion detected but seemingly good
%results
%crop = [100 210 1124 722];
%roi=[302 226 522 410];

prefsFile = '/home/athina/Data/inscopix/examplePrefs.mat';
workDir = '/home/athina/Data/inscopix/temp';
memoryQuota = 28;
driveQuota = 60;

% numMovies = length(movieFiles);
% let's start with one
numMovies = length(movieFiles);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

addpath(srcDir);

prefs = mosaic.Preferences(prefsFile, workDir, memoryQuota, driveQuota);
mosaic.initialize('preferences', prefs);

%% motion correct each movie
refName = fullfile(outDir,'referenceFrame.tif');

for m = 1:numMovies
    % load movie
    movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}));
    
    % crop movie
    croppedMovie = mosaic.cropMovie(movie, ...
        crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');
    clear movie;
    
    % preprocess movie
    ppMovie = mosaic.preprocessMovie(croppedMovie, ...
       'fixDefectivePixels', fixDefectivePixels, ...
       'fixRowNoise', fixRowNoise, ...
       'fixDroppedFrames', fixDroppedFrames, ...
       'spatialDownsampleFactor', spatialDownsampleFactor);
   clear croppedMovie;
    
    % motion correct movie
    % - define the reference frame
    if m==1
        referenceFrame = mosaic.extractFrame(ppMovie, 'frame', 1);
    else
        referenceFrame = mosaic.loadImage(refName);
    end
    
    mcRoi = mosaic.RectangleRoi(mosaic.Point(roi(1), roi(2)), mosaic.Point(roi(3), roi(4)));
    %mcRoi.view(ppMovie);
    
    [mcMovie, translations] = mosaic.motionCorrectMovie(ppMovie, ... 
    'referenceImage', referenceFrame, ...
    'motionType', 'Translation', ... 
    'roi', mcRoi, ... 
    'speedWeight', speedWeight, ... 
    'parallelProcess', parallelProcess, ... 
    'invertImage', invertImage, ... 
    'normalizeImage', normalizeImage, ... 
    'subtractSpatialMean', subtractSpatialMean, ... %true, ... 
    'subtractSpatialMeanPixels', subtractSpatialMeanPixels, ... 
    'applySpatialMean', applySpatialMean, ... %true, ... 
    'applySpatialMeanPixels', applySpatialMeanPixels, ... 
    'minimumValue', minimumValue, ... 
    'maximumValue', maximumValue, ... 
    'autoCrop', false ... % true ...
    );
    %clear ppMovie;
    
    % save motion corrected movie
    movName = fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}));
    mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'None');
    
    % save the last frame as a reference for next part of the movie
    xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
    xtranslast = round(xtrans(end));
    ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();
    % find(ytrans==max(ytrans))
    ytranslast = round(ytrans(end));
    idx = mcMovie.getTimingInfo().getNumTimes();
    lastFrame = mosaic.extractFrame(mcMovie, 'frame', idx);
    % something goes wrong here
    % check 
    lastFrame2 = mosaic.extractFrame(ppMovie, 'frame', idx);
    lastFrame2.view();
    % the coloring is off, which is probably the cause of bad motion
    % correction
    % solution: do the translation ourselves from the ppMovie. 
    % maybe mosaic has a function that facilitates this?
    mosaic.saveImageTiff(lastFrame, refName,'compression','None');
    
    clear mcMovie;
end

clear fixDefectivePixels fixRowNoise fixDroppedFrames spatialDownsampleFactor;
clear speedWeight parallelProcess invertImage normalizeImage;
clear subtractSpatialMean subtractSpatialMeanPixels applySpatialMean;
clear applySpatialMeanPixels minimumValue maximumValue;
    
clear referenceFrame idx translations;
%% load all motion corrected movies and concatenate movie them
% load tiff movies
movies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    movie = mosaic.loadMovieTiff(fullfile(outDir, sprintf('mcorr_%s', movieFiles{m})));
    movies.add(movie);
end
clear movie;

concatMovie = mosaic.concatenateMovies(movies, ...
    'gapType', 'Add one time period between movies');
clear ppMovies;

%% Temporally downsample and save as .mat for CNMF-E
% actually, CNMF-E doesn't need 5 Hz input. Let's try with normal input and
% compare.

if time_down
desiredStep = 0.2; % CNMF wants 5 Hz data input
downsamplingFactor = round(desiredStep / currentStep);

downsMovie = mosaic.resampleMovie(concatMovie, ...
    'spatialReduction', 1, ...
    'temporalReduction', downsamplingFactor);

Y = downsMovie.getData();
currentStep = desiredStep;
else
    Y = concatMovie.getData();
end
%Y = Y(:,:,1:downsamplingFactor:end);

% remove unwanted rows (caused by motion correction)
% 1. sum over columns if sum=0, then the row is blank and it should not be
%used (the motion correction algorithm made it blank).
% 2. find the minimum sub over columns for each row, over all frames (even
% if only one frame does not contain data for this row, we cannot use the
% whole row.
% rule: min(sum(Y,2),[],3)==0;
% remove unwanted columns, rule: min(sum(Y,1),[],3)==0;

% keep only rows and columns that never had sum=0.
Y = Y(find(min(sum(Y,2),[],3)~=0),find(min(sum(Y,1),[],3)~=0),:);

% save some meta-data
Ysiz = size(Y)';
Yfs = 1/currentStep;

nam_mat = fullfile(outDir, outFileName);
save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');

clear Y Ysiz Yfs;
clear currentStep desiredStep downsamplingFactor;

%% Terminate the Mosaic script session.
mosaic.terminate();