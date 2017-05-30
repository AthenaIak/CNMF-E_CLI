function [ data ] = loadRawData( filename )
%LOADRAWDATA Loads the motion-corrected data.
%   Loads data from .tif or .mat file and makes sure that the data are
%   clean, meaning that there are no rows or columns filled with zero. If
%   the data are not clean, it cleans them and saves the results for future
%   use.

% load tif or mat file
if ~exist('filename', 'var') || isempty(filename)
    filename = 'Z:\athina\data\icdc.mat';
end

[path,name,ext] = fileparts(filename);

if strcmp(ext,'.tif')
    nam_mat = tif2mat(filename);
elseif strcmp(ext,'.mat')
    nam_mat = filename;
else
    disp('Unsupported format: only .tif and .mat files allowed');
    exit(1);
end

% temporally downsample the data
data = matfile(nam_mat);
try
    Yfs = data.Yfs;
catch
    warning('Frequency not specified. Assumed frequency is 20Hz.');
    Yfs = 20;
end

fsChanged = false;
dfc = floor(Yfs/5);
if dfc ~= 1
    disp(sprintf('Temporally downsampling data from %0.1fHz by a factor of %d.\n',Yfs,dfc));
    old_num_frames = size(data.Y,3);
    Y = data.Y(:,:,1:dfc:end);
    new_num_frames = size(Y,3);
    Yfs = Yfs * new_num_frames / old_num_frames;
    fsChanged = true;
else
    Y = data.Y;
end

% remove unwanted rows
% 1. sum over columns if sum=0, then the row is blank and it should not be
%used (the motion correction algorithm made it blank).
% 2. find the minimum sub over columns for each row, over all frames (even
% if only one frame does not contain data for this row, we cannot use the
% whole row.
% rule: min(sum(Y,2),[],3)==0;
% remove unwanted columns, rule: min(sum(Y,1),[],3)==0;

% keep only rows and columns that never had sum=0.
Y = Y(find(min(sum(Y,2),[],3)~=0),find(min(sum(Y,1),[],3)~=0),:);
Ysiz = size(Y)';
if max(data.Ysiz - Ysiz) ~= 0
    disp('Data cleaned.');
end

% save the data if any changes were made
if or(fsChanged, max(data.Ysiz - Ysiz) ~= 0)
    disp('Saving data...');
    nam_mat = sprintf('%s%s%s_proc.mat',path,filesep,name);
    save(nam_mat, 'Y', 'Ysiz', 'Yfs', '-v7.3');
    disp('Data saved.');
    
    disp('Loading the changed data.');
    data = matfile(nam_mat);
    disp('Data loaded');
end

end

