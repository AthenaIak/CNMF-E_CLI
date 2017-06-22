clc; clear; close all;

srcDir = '~/tp/athina/Mosaic-MATLAB/src';
inDir = '~/tp/athina/data/4294';
outDir = '~/tp/athina/data/4294/output20160531';

% movieFiles = {...
%     %'recording_20160118_140521.xml', ...
%     'recording_20160125_114209.xml'...
% };

movieFiles = {...
    'recording_20160531_141257.tif', ...
    %'recording_20160531_141257-001.tif', ...
    %'recording_20160531_141257-002.tif', ...
    %'recording_20160531_141829.tif', ...
    %'recording_20160531_141829-001.tif', ...
    %'recording_20160531_141829-002.tif', ...
    %'recording_20160531_141829-003.tif', ...
    %'recording_20160531_141829-004.tif', ...
    %'recording_20160531_142630.tif', ...
    %'recording_20160531_142630-001.tif', ...
    %'recording_20160531_142630-002.tif', ...
    %'recording_20160531_142630-003.tif', ...
    %'recording_20160531_142630-004.tif', ...
    %'recording_20160531_143408.tif', ...
    %'recording_20160531_143408-001.tif', ...
    %'recording_20160531_143408-002.tif', ...
    %'recording_20160531_143408-003.tif', ...
    %'recording_20160531_143408-004.tif', ...
    %'recording_20160531_144143.tif', ...
    %'recording_20160531_144143-001.tif', ...
    %'recording_20160531_144143-002.tif'
};

% RoI used for motion correction
roi = [390 527 986 800];

% crop parameters
crop = [91 216 1162 614];

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


% reference frame used for motion correction
refName = fullfile(outDir,'referenceFrame.tif');
referenceFrame = mosaic.loadImage(refName);

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
    ppMovies.add(croppedMovies.get(m));
    %movName = fullfile(outDir, sprintf('preproc_%s.tif',croppedMovies.get(m).getName()));
    %mosaic.saveMovieTiff(ppMovie, movName, 'compression', 'None');
end
clear ppMovie movName movies;


%% Correct motion for each movie using a common reference frame

mcMovies = mosaic.List('mosaic.Movie');
referenceFrame = mosaic.extractFrame(ppMovies.get(1), 'frame', 1);
%mosaic.saveImageTiff(referenceFrame, refName,'compression','None');
mcRoi = mosaic.RectangleRoi(mosaic.Point(roi(1), roi(2)), mosaic.Point(roi(3), roi(4)));
%mcRoi.view(ppMovies.get(1));

for i=1:numMovies 
    [mcMovie, translations] = mosaic.motionCorrectMovie(ppMovies.get(i), ... 
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
    mcMovies.add(mcMovie);
end
clear mcMovie ppMovies;

%% concatenate movie files to one file
concatMovie = mosaic.concatenateMovies(mcMovies, ...
    'gapType', 'Add one time period between movies');
clear mcMovies;

movName = fullfile(outDir, sprintf('mcorr_%s.tif',movieFiles{1}));
mosaic.saveMovieTiff(concatMovie, movName, 'compression', 'None');
clear movName;

%% Temporally downsample and save as .mat for CNMF-E
currentStep = concatMovie.getTimingInfo().getStep();
desiredStep = 0.2; % CNMF wants 5 Hz data input
downsamplingFactor = round(desiredStep / currentStep);
Y = concatMovie.getData();
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