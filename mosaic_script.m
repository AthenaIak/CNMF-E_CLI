clc; clear; close all; 
pack;

srcDir = '~/tu/athina/Mosaic-MATLAB/src';

set_parameters='~/tu/athina/Data/32366/parameters_an006_mosaic';
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
mcRoi = mosaic.RectangleRoi(mosaic.Point(roi(1), roi(2)), mosaic.Point(roi(3), roi(4)));

maxXtrans=0; minXtrans=0; maxYtrans=0; minYtrans=0;
for m = 1:numMovies
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
   
   % extract the last 2 frames (mosaic doesn't provide a function to crop
   % an image, so we need a movie of minimum 2 frames) The point is to
   % extract the last frame after we perform motion correction.
   idx = ppMovie.getTimingInfo().getNumTimes();
   last2frames = mosaic.trimMovie(ppMovie, idx-1, idx);
   
    % crop movie
    croppedMovie = mosaic.cropMovie(ppMovie, ...
        crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');
    
    % motion correct movie
    % - define the reference frame
    if m==1
        referenceFrame = mosaic.extractFrame(croppedMovie, 'frame', 1);
        clear ppMovie;
    else
       % referenceFrame = mosaic.loadImage(refName);
        load(refName);
        refFrames = mosaic.cropMovie(last2frames_prev, ...
        crop(1)+xtranslast, crop(2)+ytranslast, ...
        crop(3)+xtranslast, crop(4)+ytranslast, 'coordinateSystem', 'pixels');
    
        referenceFrame = mosaic.extractFrame(refFrames, 'frame', 2);
        moviePart = mosaic.trimMovie(croppedMovie, 1, 100);
        
        [~, translations] = mosaic.motionCorrectMovie(moviePart, ... 
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
    
        xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
        ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();
        xtransmedian = round(median(xtrans));
        ytransmedian = round(median(ytrans));
        
        tmp = max(xtrans);
         if maxXtrans < tmp
            maxXtrans = tmp;
        end
        tmp = min(xtrans);
        if minXtrans > tmp
           minXtrans = tmp;
        end
        tmp = max(ytrans);
        if maxYtrans < tmp
           maxYtrans = tmp;
        end
        tmp = min(ytrans);
        if minYtrans > tmp
           minYtrans = tmp;
        end

        clear croppedMovie;
        croppedMovie = mosaic.cropMovie(ppMovie, ...
        crop(1)+xtransmedian, crop(2)+ytransmedian, ... 
        crop(3)+xtransmedian, crop(4)+ytransmedian, 'coordinateSystem', 'pixels');
        clear ppMovie;
    
        clear last2frames_prev moviePart translations xtras ytrans;
        clear xtransmedian ytransmedian;
    end
    
    %mcRoi.view(croppedMovie);
    
    [mcMovie, translations] = mosaic.motionCorrectMovie(croppedMovie, ... 
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
    clear croppedMovie;
    
    % save motion corrected movie
    movName = fullfile(outDir, sprintf('mcorr_%s', movieFiles{m}{1}));
    mosaic.saveMovieTiff(mcMovie, movName, 'compression', 'None');
    
    % save the last frame as a reference for next part of the movie
    xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
    ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();
    xtranslast = round(xtrans(end));
    ytranslast = round(ytrans(end));
    
    disp(sprintf('File %s\nvariation x = %f | range x = %f\nvariation y = %f | range y = %f', ...
        movieFiles{m}{1}, var(xtrans), range(xtrans), ...
        var(ytrans), range(ytrans)));
    
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