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
refName = fullfile(outDir,'referenceFrameInfo.mat');

for m = 1:numMovies
    % load movie
    movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}));
    
    % preprocess movie
    ppMovie = mosaic.preprocessMovie(movie, ...
       'fixDefectivePixels', fixDefectivePixels, ...
       'fixRowNoise', fixRowNoise, ...
       'fixDroppedFrames', fixDroppedFrames, ...
       'spatialDownsampleFactor', spatialDownsampleFactor);
   clear movie;
   
   % extract the last 2 frames (mosaic doesn't provide a function to crop
   % an image, so we need a movie of minimum 2 frames) The point is to
   % extract the last frame after we perform motion correction.
   idx = ppMovie.getTimingInfo().getNumTimes();
   last2frames = mosaic.trimMovie(ppMovie, idx-1, idx);
   
    % crop movie
    croppedMovie = mosaic.cropMovie(ppMovie, ...
        crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');
    clear ppMovie;
    
    % motion correct movie
    % - define the reference frame
    if m==1
        referenceFrame = mosaic.extractFrame(croppedMovie, 'frame', 1);
        
    else
       % referenceFrame = mosaic.loadImage(refName);
        load(refName);
        refFrames = mosaic.cropMovie(last2frames_prev, ...
        crop(1)+xtranslast, crop(2)+ytranslast, ...
        crop(3)+xtranslast, crop(4)+ytranslast, 'coordinateSystem', 'pixels');
    
        referenceFrame = mosaic.extractFrame(refFrames, 'frame', 2);
        
        moviePart = mosaic.trimMovie(croppedMovie,1,100);
    
        [~, translations] = mosaic.motionCorrectMovie(moviePart, ... 
            'referenceImage', referenceFrame, ...
            'motionType', 'Translation', ... 
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
    
        xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
        ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();
        xtransmedian = round(median(xtrans));
        ytransmedian = round(median(ytrans));
    
        refFrames = mosaic.cropMovie(last2frames_prev, ...
        crop(1)+xtransmedian, crop(2)+ytransmedian, ...
        crop(3)+xtransmedian, crop(4)+ytransmedian, 'coordinateSystem', 'pixels');
    
        referenceFrame = mosaic.extractFrame(refFrames, 'frame', 2);
        
        clear last2frames_prev moviePart translations xtras ytrans;
        clear xtransmedian ytransmedian;
    end
    
    mcRoi = mosaic.RectangleRoi(mosaic.Point(roi(1), roi(2)), mosaic.Point(roi(3), roi(4)));
    %mcRoi.view(croppedMovie);
    
    [mcMovie, translations] = mosaic.motionCorrectMovie(croppedMovie, ... 
    'referenceImage', referenceFrame, ...
    'motionType', 'Translation', ... 
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
    clear croppedMovie;
    
    % save motion corrected movie
    movName = fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}));
    mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'None');
    
    % save the last frame as a reference for next part of the movie
    xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
    ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();
    xtranslast = round(xtrans(end));
    ytranslast = round(ytrans(end));
    
    %disp(min(xtrans));disp( max(xtrans));
    disp(var(xtrans));disp(range(xtrans));
    %disp(min(ytrans));disp( max(ytrans));
    disp(var(ytrans));disp(range(ytrans));
    last2frames_prev = last2frames;
    save(refName, 'last2frames_prev', 'xtranslast', 'ytranslast','-v7.3');
    
    clear translations xtrans xtranslast ytrans ytranslast;
    clear last2frames mcFrames lastFrame;
    clear mcMovie; 
end

clear fixDefectivePixels fixRowNoise fixDroppedFrames spatialDownsampleFactor;
clear speedWeight parallelProcess invertImage normalizeImage;
clear subtractSpatialMean subtractSpatialMeanPixels applySpatialMean;
clear applySpatialMeanPixels minimumValue maximumValue;
    
clear referenceFrame idx translations;

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
currentStep = concatMovie.getTimingInfo().getStep();
downsamplingFactor = round(desiredStep / currentStep);

downsMovie = mosaic.resampleMovie(concatMovie, ...
    'spatialReduction', 1, ...
    'temporalReduction', downsamplingFactor);

Y = downsMovie.getData();
currentStep = desiredStep;
else
    Y = concatMovie.getData();
end
clear concatMovie;
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