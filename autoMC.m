clc; clear; close all;
pack;

srcDir = '~/tu/athina/Mosaic-MATLAB/src';

set_parameters='~/tu/athina/Data/32366/parameters_an009_mosaic';
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
xw = 180; % window width, x axis
yw = 180; % widnow width, y axis
pw = 20;  % padding width

ROIx = [pw:xw:1440;xw-pw:xw:1440]';
ROIy = [pw:yw:1080;yw-pw:yw:1080]';

rois = [repmat(ROIx,size(ROIy,1),1) repelem(ROIy,size(ROIx,1),1)];
rois(:,[2,3]) = rois(:,[3,2]);
nComb = size(rois,1);

%% motion correct each movie
refName = fullfile(outDir,'referenceFrameInfo.mat');

maxXtrans=0; minXtrans=0; maxYtrans=0; minYtrans=0;
m = 1; %no loop over sessions (for now)


numFiles = length(movieFiles{m});
% load movie files
movies = mosaic.List('mosaic.Movie');
for f = 1:numFiles
    movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}{f}));
    movies.add(movie);
end
clear movie;

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

transPerROI = zeros(nComb, ppMovie.getTimingInfo().getNumTimes()*2);
for r=1:nComb
    mcRoi = mosaic.RectangleRoi(mosaic.Point(rois(r,1), rois(r,2)), ...
        mosaic.Point(rois(r,3), rois(r, 4)));
    %mcRoi.view(ppMovie);
    [mcMovie, translations] = mosaic.motionCorrectMovie(ppMovie, ...
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

save(fullfile(outDir,'traslationsPerROI.mat'),'transPerROI','-v7.3');
clear croppedMovie;

% save motion corrected movie
movName = fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}{1}));
mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'None');

% save the last frame as a reference for next part of the movie

clear translations xtrans xtranslast ytrans ytranslast;
clear last2frames mcFrames lastFrame;
clear mcMovie;

clear fixDefectivePixels fixRowNoise fixDroppedFrames spatialDownsampleFactor;
clear speedWeight parallelProcess invertImage normalizeImage;
clear subtractSpatialMean subtractSpatialMeanPixels applySpatialMean;
clear applySpatialMeanPixels minimumValue maximumValue;

clear referenceFrame idx translations;

%% remove empty space caused by motion correction
for m=1:numMovies
    movie = mosaic.loadMovieTiff(fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}{1})));
    
    datasize = movie.getDataSize(); % returns Y X T ?
    
    % crop movie
    croppedMovie = mosaic.cropMovie(movie, ...
        floor(maxXtrans)+1, floor(maxYtrans)+1, ...
        datasize(2)+floor(minXtrans)-1, ...
        datasize(1)+floor(minYtrans)-1, 'coordinateSystem', 'pixels');
    
    movName = fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}{1}));
    mosaic.saveMovieTiff(croppedMovie, movName, 'compression', 'None');
end

%% remove border of video


%% clear matlab's memory
% by this time, matlab has allocated a lot of space that it didn't free. we
% need to free matlab memory
pack;
mosaic.initialize('preferences', prefs);
% if this fails, then run the first section (starts with clc; clear; close
% all;) move forward, skipping the 2nd section (motion correction).

%% load all motion corrected movies and concatenate movie them
% load tiff movies
movies = mosaic.List('mosaic.Movie');
for m = 1:numMovies
    movie = mosaic.loadMovieTiff(fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}{1})));
    movies.add(movie);
end
clear movie;

concatMovie = mosaic.concatenateMovies(movies, ...
    'gapType', 'Add one time period between movies');
clear ppMovies;

%% Save as .mat for CNMF-E

currentStep = concatMovie.getTimingInfo().getStep();
Y = concatMovie.getData();

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
Yfs = (1/currentStep)/down_factor;

nam_mat = fullfile(outDir, outFileName);
save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');

clear Y Ysiz Yfs;
clear currentStep desiredStep downsamplingFactor;

%% Terminate the Mosaic script session.
mosaic.terminate();