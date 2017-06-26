clc; clear; close all;

srcDir = '~/tp/athina/Mosaic-MATLAB/src';
inDir = '~/tp/athina/data/4294/';
outDir = '~/tp/athina/data/4294/output/20170626_of_20160531';

outFileName = 'ready_20160531.mat';
movieFiles = {...
    'recording_20160531_141257.tif', ...
    'recording_20160531_141257-001.tif', ...
    'recording_20160531_141257-002.tif', ...
    'recording_20160531_141829.tif', ...
    'recording_20160531_141829-001.tif', ...
    'recording_20160531_141829-002.tif', ...
    'recording_20160531_141829-003.tif', ...
    'recording_20160531_141829-004.tif', ...
    'recording_20160531_142630.tif', ...
    'recording_20160531_142630-001.tif', ...
    'recording_20160531_142630-002.tif', ...
    'recording_20160531_142630-003.tif', ...
    'recording_20160531_142630-004.tif', ...
    'recording_20160531_143408.tif', ...
    'recording_20160531_143408-001.tif', ...
    'recording_20160531_143408-002.tif', ...
    'recording_20160531_143408-003.tif', ...
    'recording_20160531_143408-004.tif', ...
    'recording_20160531_144143.tif', ...
    'recording_20160531_144143-001.tif', ...
    'recording_20160531_144143-002.tif'
};

% crop parameters
crop = [100 210 1196 796];

% RoI used for motion correction
roi = [253 644 576 830];

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
for m = 1:numMovies
    % load movie
    movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}));
    
    % crop movie
    croppedMovie = mosaic.cropMovie(movie, ...
        crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');
    clear movie;
    
    % preprocess movie
    ppMovie = mosaic.preprocessMovie(croppedMovie, ...
       'fixDefectivePixels', true, ...
       'fixRowNoise', false, ...
       'fixDroppedFrames', false, ...
       'spatialDownsampleFactor', 1); % do not downsample
    clear croppedMovie;
    
    % motion correct movie
    % - define the reference frame
    refName = fullfile(outDir,'referenceFrame.tif');
    if m==1
        referenceFrame = mosaic.extractFrame(ppMovie, 'frame', 1);
    else
        referenceFrame = mosaic.loadImage(refName);
    end
    
    mcRoi = mosaic.RectangleRoi(mosaic.Point(roi(1), roi(2)), mosaic.Point(roi(3), roi(4)));
    %mcRoi.view(concatMovie);
    
    [mcMovie, translations] = mosaic.motionCorrectMovie(ppMovie, ... 
    'referenceImage', referenceFrame, ...
    'motionType', 'Translation', ... 
    'roi', mcRoi, ... 
    'speedWeight', 0.1, ... 
    'parallelProcess', true, ... 
    'invertImage', true, ... 
    'normalizeImage', false, ... 
    'subtractSpatialMean', true, ... %true, ... 
    'subtractSpatialMeanPixels', 20, ... 
    'applySpatialMean', true, ... %true, ... 
    'applySpatialMeanPixels', 5, ... 
    'minimumValue', -102, ... 
    'maximumValue', 108, ... 
    'autoCrop', false ... % true ...
    );
    clear ppMovie;
    
    % save motion corrected movie
    movName = fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}));
    mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'None');
    
    % save the last frame as a reference for next part of the movie
    idx = mcMovie.getTimingInfo().getNumTimes();
    lastFrame = mosaic.extractFrame(mcMovie, 'frame', idx);
    mosaic.saveImageTiff(lastFrame, refName,'compression','None');
    
    clear mcMovie;
end

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
currentStep = concatMovie.getTimingInfo().getStep();
desiredStep = 0.2; % CNMF wants 5 Hz data input
downsamplingFactor = round(desiredStep / currentStep);

downsMovie = mosaic.resampleMovie(concatMovie, ...
    'spatialReduction', 1, ...
    'temporalReduction', downsamplingFactor);

Y = downsMovie.getData();
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
Yfs = 1/desiredStep;

nam_mat = fullfile(outDir, outFileName);
save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');

clear Y Ysiz Yfs;
clear currentStep desiredStep downsamplingFactor;

%% Terminate the Mosaic script session.
mosaic.terminate();