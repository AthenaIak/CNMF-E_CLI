function [ ans ] = run_analysis( filename )
%RUN_ANALYSIS runs CNMF-E analysis on the data
%   Filename: the path+name of the file containing the data to be analyzed.

% add paths
[CNMF_dir,~,~]= fileparts(which('run_analysis'));
addpath(sprintf('%s%sca_source_extraction%s', CNMF_dir, filesep, filesep));
addpath(sprintf('%s%sca_source_extraction%sutilities%s', CNMF_dir, filesep, filesep, filesep));
addpath(sprintf('%s%sca_source_extraction%sendoscope', CNMF_dir, filesep, filesep, filesep));
addpath(sprintf('%s%sCLI', CNMF_dir, filesep));
addpath(sprintf('%s%soasis', CNMF_dir, filesep));
addpath(sprintf('%s%scnmfe_scripts', CNMF_dir, filesep));

addpath(sprintf('%s%sGUI', CNMF_dir, filesep));
addpath(sprintf('%s%sGUI%sgui_callbacks', CNMF_dir, filesep, filesep));
addpath(sprintf('%s%sGUI%smodules', CNMF_dir, filesep, filesep));
clear CNMF_dir;

data = loadRawData(filename);




