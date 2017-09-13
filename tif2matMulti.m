function nam_mat = tif2matMulti(movieFiles, bd, nam_mat)
%% convert tiff files into mat files
% inputs:
%   movieFiles: file names
%   bd: boundaries of the movie to be cropped.
%   nam_mat: name of the *.mat file to be created (optional).
%   file saved with the indicated name.
% output:
%   nam_mat: name of the *.mat file

%% tiff file information 
if isempty(movieFiles)
    nam_mat = [];
    fprintf('empty tiff file');
    return;
end

if ~exist('nam_mat','var')
    [tmp_dir, tmp_file, ~] = fileparts(movieFiles{1});
    nam_mat = sprintf('%s%s%s-all.mat', tmp_dir, filesep, tmp_file);
end

if ~exist('bd','var')
    bd = 0;
end

tic;

info = imfinfo(movieFiles{1});
d1 = info.Height;   % height of the image 
d2 = info.Width;    % width of the image 
d1 = d1 -2*bd;
d2 = d2 -2*bd;
T = length(info);   % number of frames 
for m=2:length(movieFiles)
    info = imfinfo(movieFiles{m});
    T= T + length(info);    
end
Ysiz = [d1, d2, T]'; 

fprintf('CNMF_E is converting TIFF file to *.mat file\n'); 
curr_idx = 1;
numFiles = length(movieFiles);
for m=1:numFiles
    fprintf('Loading file %d/%d...\n', m,numFiles);
    Y = bigread2(movieFiles{m}); 
    Y = Y(bd+1:end-bd,bd+1:end-bd,:);
    num_frames = size(Y,3);
    
    disp('Saving file...');
    if m==1
        % create a mat file 
        save(nam_mat, 'Y', 'Ysiz', '-v7.3'); 
        data = matfile(nam_mat, 'Writable', true); 
    else
        data.Y(:,:,(1:num_frames)+curr_idx) = Y;
    end
    
    curr_idx = curr_idx + num_frames;
end

fprintf('Time cost of converting tiffs to mat: %.1f sec.\n', toc);    
