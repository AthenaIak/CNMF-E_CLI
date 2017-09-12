clc; clear; close all;
pack;

srcDir = '~/tu/athina/Mosaic-MATLAB/src';

set_parameters='~/tu/athina/Data/analyzed/parameters/parameters_an001_mosaic';
run (set_parameters);

prefsFile = '/home/athina/Data/inscopix/examplePrefs.mat';
workDir = '/home/athina/Data/inscopix/temp';
memoryQuota = 28;
driveQuota = 60;

numMovies = length(movieFiles);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

addpath(srcDir);

prefs = mosaic.Preferences(prefsFile, workDir, memoryQuota, driveQuota);
mosaic.initialize('preferences', prefs);

%% rois to test

ROIx = [pw:xw:1440;xw-pw:xw:1440]';
ROIy = [pw:yw:1080;yw-pw:yw:1080]';

rois = [repmat(ROIx,size(ROIy,1),1) repelem(ROIy,size(ROIx,1),1)];
rois(:,[2,3]) = rois(:,[3,2]);
nComb = size(rois,1);

clear xw yw pw ROIx ROIy;

%% preprocess the session
%refName = fullfile(outDir,'referenceFrameInfo.mat');

%maxXtrans=0; minXtrans=0; maxYtrans=0; minYtrans=0;
m = 1; %no loop over sessions (for now)


numFiles = length(movieFiles{m});
% load movie files
movies = mosaic.List('mosaic.Movie');
for f = 1:numFiles
    movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}{f}));
    movies.add(movie);
end
clear movie f;

% concat session into one movie
session = mosaic.concatenateMovies(movies, ...
    'gapType', 'Add one time period between movies');
clear movies;

% preprocess movie
ppMovie = mosaic.preprocessMovie(session, ...
    'fixDefectivePixels', fixDefectivePixels, ...
    'fixRowNoise', fixRowNoise, ...
    'fixDroppedFrames', fixDroppedFrames, ...
    'spatialDownsampleFactor', spatialDownsampleFactor);
clear session;

% temporally downsample movie
if time_down
    ppMovie = mosaic.resampleMovie(ppMovie, ...
        'spatialReduction', 1, ...
        'temporalReduction', down_factor);
end

if m==1
    referenceFrame = mosaic.extractFrame(ppMovie, 'frame', 1);
end

% save movie that is ready for motion correction
movName = fullfile(outDir, sprintf('proc_%s', movieFiles{m}{1}));
mosaic.saveMovieTiff(ppMovie, movName, 'compression', 'None');
    
%% correct motion using different ROIs
transPerROI = zeros(nComb, ppMovie.getTimingInfo().getNumTimes()*2);
for r=1:nComb
    mcRoi = mosaic.RectangleRoi(mosaic.Point(rois(r,1), rois(r,2)), ...
        mosaic.Point(rois(r,3), rois(r, 4)));
    %mcRoi.view(ppMovie);
    [~, translations] = mosaic.motionCorrectMovie(ppMovie, ...
        'referenceImage', referenceFrame, ...
        'motionType', motionType, ...
        'roi', mcRoi, ...
        'speedWeight', speedWeight, ...
        'parallelProcess', parallelProcess, ...
        'invertImage', invertImage, ...
        'normalizeImage', normalizeImage, ...
        'subtractSpatialMean', subtractSpatialMean, ...
        'subtractSpatialMeanPixels', subtractSpatialMeanPixels, ...
        'applySpatialMean', applySpatialMean, ...
        'applySpatialMeanPixels', applySpatialMeanPixels, ...
        'minimumValue', minimumValue, ...
        'maximumValue', maximumValue, ...
        'autoCrop', false ... % true ...
        );
    transPerROI(r,:) = [translations.getList('types', {'mosaic.Trace'}).get(1).getData();translations.getList('types', {'mosaic.Trace'}).get(1).getData()]';
    
    xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
    ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();
    
    disp(sprintf('ROI %d\nvariation x = %f | range x = %f\nvariation y = %f | range y = %f', ...
        r, var(xtrans), range(xtrans), ...
        var(ytrans), range(ytrans)));
end

% save the last frame as a reference for next session
clear ppMovie mcRoi xtrans ytrans;

save(fullfile(outDir,'traslationsPerROI.mat'),'transPerROI','-v7.3');
clear croppedMovie;

clear translations xtrans xtranslast ytrans ytranslast;
clear last2frames mcFrames lastFrame;
clear mcMovie;

clear fixDefectivePixels fixRowNoise fixDroppedFrames spatialDownsampleFactor;
clear speedWeight parallelProcess invertImage normalizeImage;
clear subtractSpatialMean subtractSpatialMeanPixels applySpatialMean;
clear applySpatialMeanPixels minimumValue maximumValue;

clear referenceFrame idx translations;

%% Terminate the Mosaic script session.
mosaic.terminate();