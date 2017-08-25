function [ ] = run_analysis( set_parameters )

srcDir = '~/tu/athina/Mosaic-MATLAB/src';

%set_parameters='~/tu/athina/Data/analyzed/parameters/parameters_preproc_32363_day1';
run (set_parameters);

prefsFile = '/home/athina/Data/inscopix/examplePrefs.mat';
workDir = '/home/athina/Data/inscopix/temp';

numMovies = length(movieFiles);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

addpath(srcDir);

prefs = mosaic.Preferences(prefsFile, workDir, memoryQuota, driveQuota);
mosaic.initialize('preferences', prefs);

%% preprocess the session
for m=1:numMovies
    
    disp(sprintf('TRIAL #%d',m));

    numFiles = length(movieFiles{m});
    
    % load movie files
    disp('Loading files...');
    movies = mosaic.List('mosaic.Movie');
    for f = 1:numFiles
        movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}{f}));
        movies.add(movie);
        disp(sprintf('Loaded %s', movieFiles{m}{f}));
    end
    clear movie f;
    disp('Done loading');

    % concat trial into one movie
    disp('Concatenating files...');
    trial = mosaic.concatenateMovies(movies, ...
        'gapType', 'Add one time period between movies');
    clear movies;
    disp('Done concatenating');

    % preprocess trial
    disp('Preprocessing trial...');
    ppMovie = mosaic.preprocessMovie(trial, ...
        'fixDefectivePixels', fixDefectivePixels, ...
        'fixRowNoise', fixRowNoise, ...
        'fixDroppedFrames', fixDroppedFrames, ...
        'spatialDownsampleFactor', spatialDownsampleFactor);
    clear trial;
    disp('Done preprocessing');

    % temporally downsample movie
    disp('Temporally downsampling...');
    if time_down
        ppMovie = mosaic.resampleMovie(ppMovie, ...
            'spatialReduction', 1, ...
            'temporalReduction', down_factor);
    end
    disp('Done downsampling');

    % save trial that is ready for motion correction
    disp('Saving preprocessed movie...');
    movName = fullfile(outDir, sprintf('pp_%s', movieFiles{m}{1}));
    mosaic.saveMovieTiff(ppMovie, movName, 'compression', 'None');
    disp(sprintf('Saved as %s', movName));
end

%% Terminate the Mosaic script session.
mosaic.terminate();