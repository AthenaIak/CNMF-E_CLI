clc; clear; close all;
pack;

srcDir = '~/tu/athina/Mosaic-MATLAB/src';

set_parameters='~/tu/athina/Data/analyzed/parameters/parameters_an013_mosaic';
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

%% find best ROI

% load translations
load(fullfile(outDir,'traslationsPerROI.mat'));

minError = inf; bestIdxA = 0; bestIdxB = 0;
for i=1:size(transPerROI,1)
    error = sqrt(sum(bsxfun(@minus,transPerROI([1:i-1,i+1:end],:),transPerROI(i,:)).^2,2));
    [val,idx] = min(error);
    if idx >= i
        idx = idx +1;
    end
    disp(sprintf('Pair: %d-%d, error: %f', i, idx, val));
    
    if val<minError
        minError = val;
        bestIdxA = i;
        bestIdxB = idx;
    end
end

disp(sprintf('Best pair: %d-%d, error: %f', bestIdxA, bestIdxB, minError));

%% rois

ROIx = [pw:xw:1440;xw-pw:xw:1440]';
ROIy = [pw:yw:1080;yw-pw:yw:1080]';

rois = [repmat(ROIx,size(ROIy,1),1) repelem(ROIy,size(ROIx,1),1)];
rois(:,[2,3]) = rois(:,[3,2]);
nComb = size(rois,1);

clear xw yw pw ROIx ROIy;

%% calculate best motion correction

movies = mosaic.List('mosaic.Movie');
for m = 1:size(procSessionFiles,2)
    movie = mosaic.loadMovieTiff(fullfile(outDir, procSessionFiles{m}));
    movies.add(movie);
end
clear movie;

ppMovie = mosaic.concatenateMovies(movies, ...
    'gapType', 'Add one time period between movies');
clear movies;

% calculate again
mcRoi = mosaic.RectangleRoi(mosaic.Point(rois(bestIdxA,1), rois(bestIdxA,2)), ...
        mosaic.Point(rois(bestIdxA,3), rois(bestIdxA, 4)));
%mcRoi.view(ppMovie);
referenceFrame = mosaic.extractFrame(ppMovie, 'frame', 1);
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
clear ppMovie;

xtrans = translations.getList('types', {'mosaic.Trace'}).get(1).getData();
ytrans = translations.getList('types', {'mosaic.Trace'}).get(2).getData();

disp(sprintf('ROI %d\nvariation x = %f | range x = %f\nvariation y = %f | range y = %f', ...
    bestIdxA, var(xtrans), range(xtrans), ...
    var(ytrans), range(ytrans)));

%crop = [20 20 1140 840];

% crop movie
endMovie = mosaic.cropMovie(mcMovie, ...
    crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');

%% output for CNMF

% endMovie.view();
% move to parameters script
badFrames = [117 118 119 120];
%
    
currentStep = endMovie.getTimingInfo().getStep();
Y = endMovie.getData();
Y(:,:,badFrames)=0;
Ysiz = size(Y)';
Yfs = 20/down_factor;

nam_mat = fullfile(outDir, outFileName);
save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');

disp(sprintf('Ouput: see file %s', nam_mat));

%% Terminate the Mosaic script session.
mosaic.terminate();