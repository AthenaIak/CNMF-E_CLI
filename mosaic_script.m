clc; clear; close all;

srcDir = '/home/athina/Software/Mosaic-MATLAB/src';
inDir = '/home/athina/Data/iHPC5 raw';
outDir = '/home/athina/Data/iHPC5 raw/output';

% movieFiles = {...
%     %'recording_20160118_140521.xml', ...
%     'recording_20160125_114209.xml'...
% };

movieFiles = {...
    'recording_20160118_140521.tif', ...
    %'recording_20160125_114209.tif'...
};

% RoI for each recording (used for motion correction)
ROIs = {
    [390 370 ; 1310 930]
};

prefsFile = '/home/athina/Data/inscopix/examplePrefs.mat';
workDir = '/home/athina/Data/inscopix/temp';
memoryQuota = 16;
driveQuota = 40;

% numMovies = length(movieFiles);
% let's start with one
numMovies = 1;

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
clear movie

%% Preprocess the data to reduce data volume and remove artifacts.
ppMovies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    ppMovie = mosaic.preprocessMovie(movies.get(m), ...
       'fixDefectivePixels', true, ...
       'fixRowNoise', false, ...
       'fixDroppedFrames', false, ...
       'spatialDownsampleFactor', 1); % do not downsample
    ppMovies.add(ppMovie);
    ppMovies.add(movies.get(m));
    movName = fullfile(outDir, sprintf('preproc_%s.tif',movies.get(m).getName()));
    mosaic.saveMovieTiff(ppMovie, movName, 'compression', 'LZW');
end
clear ppMovie; clear movName;

%% Correct motion in each movie.

mcMovies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    referenceFrame = mosaic.extractFrame(ppMovies.get(m), 'frame', 1);
    mcRoi = mosaic.RectangleRoi(mosaic.Point(ROIs{m}(1,1), ROIs{m}(1,2)), mosaic.Point(ROIs{m}(2,1), ROIs{m}(2,2)));
    
    [mcMovie, translations] = mosaic.motionCorrectMovie(ppMovies.get(1), ...
    'motionType', 'Translation', ...
    'roi', mcRoi, ...
    'speedWeight', 0.1, ...
    'parallelProcess', true, ...
    'invertImage', true, ...
    'normalizeImage', false, ...
    'subtractSpatialMean', false, ... %true, ...
    'subtractSpatialMeanPixels', 20, ...
    'applySpatialMean', false, ... %true, ...
    'applySpatialMeanPixels', 5, ...
    'minimumValue', -102, ...
    'maximumValue', 108, ...
    'autoCrop', true);

    mcMovies.add(mcMovie);
    
    movName = fullfile(outDir, sprintf('mcorr_%s.tif',movies.get(m).getName()));
    mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'LZW');
end
clear ppMovies; clear mcMovie; clear movName;

%% Normalize the movie by the mean frame (df/f).
dffMovies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    dffMovie = mosaic.normalizeMovie(mcMovies.get(m), 'method', '(f-f0)/f0');
    dffMovies.add(dffMovie);
end
clear mcMovies; clear dffMovie;

% %% Save the processed movie for later.
% dffMovieFile = fullfile(outDir, 'processedMovie.mat');
% mosaic.saveOneObject(dffMovie, dffMovieFile);

%% Identify cells using PCA-ICA.
ICSs = mosaic.List('mosaic.Group');
for m=1:numMovies
    ics = mosaic.pcaIca(dffMovies.get(m), 'unmix', 'traces', ...
        'numPCs', 200, 'numICs', 150, 'icaMaxIterations', 800);
    %icTraces = ics.getList('types', {'mosaic.Trace'});
    %icMixingImages = ics.getMixingImages();
    %icUnmixingImages = ics.getUnmixingImages();
    
    ICSs.add(ics);
end

%% Save the ICs.
for m=1:numMovies
    icsFile = fullfile(outDir, 'ICSs_All.mat');
    mosaic.saveList(ICSs, icsFile);
end
%% 
icTracesAll = mosaic.List('mosaic.List');
icMixingImagesAll = mosaic.List('mosaic.List');
icUnmixingImagesAll = mosaic.List('mosaic.List');
for m=1:numMovies
    
    icTraces = ICSs.get(m).getList('types', {'mosaic.Trace'});
    icMixingImages = ICSs.get(m).getMixingImages();
    icUnmixingImages = ICSs.get(m).getUnmixingImages();
    
    icTracesAll.add(icTraces);
    icMixingImagesAll.add(icMixingImages);
    icUnmixingImagesAll.add(icUnmixingImages);
end
clear icTraces; clear icMixingImages; clear icUnmixingImages;

%% view the last three traces of the IC found
for m=1:numMovies
    numICs = icTracesAll.get(m).getSize();
    for n=numICs-2:numICs
        icTracesAll.get(m).get(n).view();
    end
end

%% Detect events in IC traces.
eventTracesAll = mosaic.List('mosaic.List');
for m=1:numMovies
    icTraces = icTracesAll.get(m);
    eventTraces = mosaic.List('mosaic.Trace');
    for i = 1:icTraces.getSize()
        eventTrace = mosaic.detectEvents(icTraces.get(i));
        eventTraces.add(eventTrace);
    end
    
    eventTracesAll.add(eventTraces);
end
clear eventTraces; clear icTraces;

%% Save the events.
for m=1:numMovies
    eventsFile = fullfile(outDir, sprintf('ICs-events_%s.mat',movies.get(m).getName()));
    mosaic.saveList(eventTracesAll.get(m), eventsFile);
end

%% Terminate the Mosaic script session.
mosaic.terminate();