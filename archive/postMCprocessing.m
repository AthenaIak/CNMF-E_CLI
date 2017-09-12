srcDir = '~/tu/athina/Mosaic-MATLAB/src';

set_parameters='~/tu/athina/Data/analyzed/parameters/parameters_preproc2_32363_day1';
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

%% load trial

for m=1:numMovies
numFiles = length(movieFiles{m});
movies = mosaic.List('mosaic.Movie');
    for f = 1:numFiles
        movie = mosaic.loadMovieTiff(fullfile(inDir, movieFiles{m}{f}));
        movies.add(movie);
    end
    clear movie;
    
    % concat session into one movie
    trial = mosaic.concatenateMovies(movies, ...
        'gapType', 'Add one time period between movies');
    clear movies;
    
    croppedMovie = mosaic.cropMovie(trial, ...
        crop(1), crop(2), crop(3), crop(4), 'coordinateSystem', 'pixels');
    
    Y = croppedMovie.getData();

    for badSet=1:size(badFrames{m},1)
        prevFrame = badFrames{m}{badSet}(1) - 1;
        nextFrame = badFrames{m}{badSet}(end) + 1;
        numBad = size(badFrames{m}{badSet},2);

        weights = linspace(1,0,numBad+2);
        weights = weights(2:end-1);
        for bf = 1:numBad
            Y(:,:,prevFrame+bf) = Y(:,:,prevFrame)*weights(bf)+...
            Y(:,:,nextFrame)*(1-weights(bf));
        end
    end
    clear badSet weights bf;

    Ysiz = size(Y)';
    Yfs = fs;

    nam_mat = fullfile(outDir, outFileNames{m});
    save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');
        
end

%% Terminate the Mosaic script session.
mosaic.terminate();