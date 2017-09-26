function[moviefiles] = getMovieFiles(movieID, numFiles)
%% getMovieFiles 
% input:
%   movieID     The ID of the movie (e.g. recording_20170710_132114-001.tif
%               has ID 20170710_132114)
%   numFiles    The number of .tif files that were generated from the .raw
%               file.
% output:
%   moviefiles  A cell(1,numFiles) containing the names of the .tif files.


moviefiles = cell(1,numFiles);
moviefiles{1} = sprintf('recording_%s.tif',movieID);
for i=2:numFiles
    moviefiles{i} = sprintf('recording_%s-00%d.tif',movieID,i-1);
end
