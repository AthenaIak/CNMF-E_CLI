clc; clear; close all;

srcDir = '~/tp/athina/Mosaic-MATLAB/src';
inDir = '~/tp/athina/data/4294/';
outDir = '~/tp/athina/data/4294/output/20170623_of_20160531_143408';

% movieFiles = {...
%     %'recording_20160118_140521.xml', ...
%     'recording_20160125_114209.xml'...
% };

movieFiles = {...
    'recording_20160531_141257.tif', ...
    %'recording_20160531_141257-001.tif', ...
    %'recording_20160531_141257-002.tif', ...
    %failed to correct motion:
    %'recording_20160531_143408.tif',
    %'recording_20160531_143408-001.tif',
    %'recording_20160531_143408-002.tif',
    %'recording_20160531_143408-003.tif',
    %'recording_20160531_143408-004.tif'
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

%% load tiff movies

movies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}));
    movies.add(movie);
end
clear movie;

%% crop movies

croppedMovies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    croppedMovie = mosaic.cropMovie(movies.get(m), ...
        crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');
    croppedMovies.add(croppedMovie);
end
clear movies croppedMovie;

%% Preprocess the data to reduce data volume and remove artifacts.
ppMovies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    ppMovie = mosaic.preprocessMovie(croppedMovies.get(m), ...
       'fixDefectivePixels', true, ...
       'fixRowNoise', false, ...
       'fixDroppedFrames', false, ...
       'spatialDownsampleFactor', 1); % do not downsample
    ppMovies.add(ppMovie);
    %movName = fullfile(outDir, sprintf('preproc_%s.tif',croppedMovies.get(m).getName()));
    %mosaic.saveMovieTiff(ppMovie, movName, 'compression', 'None');
end
clear ppMovie movName croppedMovies;

%% concatenate movie files to one file
concatMovie = mosaic.concatenateMovies(ppMovies, ...
    'gapType', 'Add one time period between movies');
clear ppMovies;

%% Correct motion in each movie.

%mcMovies = mosaic.List('mosaic.Movie');
referenceFrame = mosaic.extractFrame(concatMovie, 'frame', 1);
mcRoi = mosaic.RectangleRoi(mosaic.Point(roi(1), roi(2)), mosaic.Point(roi(3), roi(4)));
%mcRoi.view(concatMovie);

[mcMovie, translations] = mosaic.motionCorrectMovie(concatMovie, ... 
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
    'autoCrop', true);
clear concatMovie;
    
movName = fullfile(outDir, sprintf('mcorr_%s.tif',movieFiles{1}));
mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'None');

clear concatMovie movName;

%% Temporally downsample and save as .mat for CNMF-E
currentStep = mcMovie.getTimingInfo().getStep();
desiredStep = 0.2; % CNMF wants 5 Hz data input
downsamplingFactor = round(desiredStep / currentStep);
Y = mcMovie.getData();
Y = Y(:,:,1:downsamplingFactor:end);

Ysiz = size(Y)';
Yfs = 1/desiredStep;

nam_mat = fullfile(outDir, sprintf('ready_%s.mat',movieFiles{1}));
save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');

clear Y Ysiz Yfs;
clear currentStep desiredStep downsamplingFactor;

%% Normalize the movie by the mean frame (df/f).
dffMovie = mosaic.normalizeMovie(mcMovie, 'method', '(f-f0)/f0');
clear mcMovie;
%dffMovie.view()

dffMovieFile = fullfile(outDir, 'processedMovie.mat');
%mosaic.saveOneObject(dffMovie, dffMovieFile);

%% Identify cells using PCA-ICA.

ics = mosaic.pcaIca(dffMovie, 'unmix', 'traces', ...
    'numPCs', 350, 'numICs', 270);
% save object to view later
icsFile = fullfile(outDir, 'processedMovie-ICs.mat');
mosaic.saveOneObject(ics, icsFile);

% view traces etc
icTraces = ics.getList('types', {'mosaic.Trace'});
% icTraces.get(1).view()
icMixingImages = ics.getMixingImages();
icUnmixingImages = ics.getUnmixingImages();
% icMixingImages.get(1).view()

%% Detect events in IC traces.
eventTraces = mosaic.List('mosaic.Trace');
for i = 1:icTraces.getSize()
    eventTrace = mosaic.detectEvents(icTraces.get(i));
    eventTraces.add(eventTrace);
end
% eventTraces.get(1).view()

%% Save the events.
mosaic.saveList(eventTraces, eventsFile);

%% Load a previous list of traces
eventsFile = fullfile(outDir, 'processedMovie-ICs-events.mat');
eventTraces = mosaic.loadObjects(eventsFile).getList();
%% Terminate the Mosaic script session.
mosaic.terminate();