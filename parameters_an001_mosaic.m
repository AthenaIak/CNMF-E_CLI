outFileName = 'ready_20160531-001.mat';  % unique analysis ID 001

inDir = '~/tp/athina/data/4294/';
outDir = '~/tp/athina/data/4294/output/20170626_of_20160531/';

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


% preprocess movie
fixDefectivePixels = true;
fixRowNoise = false;
fixDroppedFrames = false;
spatialDownsampleFactor = 1; % 1: do not downsample

% motion correction
motionType = 'Translation';
speedWeight = 0.1;
parallelProcess = true;
invertImage = true;
normalizeImage = false;
subtractSpatialMean = true; 
subtractSpatialMeanPixels = 20;
applySpatialMean = true; 
applySpatialMeanPixels = 5;
minimumValue = -102;
maximumValue = 108; 

% temporally downsample to 5Hz
time_down = true;

