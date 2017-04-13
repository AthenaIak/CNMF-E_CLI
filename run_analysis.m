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

% filename = 'D:\Διπλωματική\CNMF_E\demos\data_endoscope.tif'
data = loadRawData(filename);
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames
fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));

% create helper structure neuron_raw
neuron_raw = neuronForData(d1,d2);
[Y, neuron, neuron_raw] = downsample(neuron_raw, data);

% compute correlation image and peak-to-noise ratio image.
% this step is not necessary, but it can give you some ideas of how your
% data look like
[Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, numFrame, 1000))));
% save to be viewed where GUI is available
[path,name,ext] = fileparts(filename);
disp('Saving correlation image and peak-to-noise ratio image...');
nam_mat = sprintf('%s%s%s%sf01-cn&pnr.mat',path,filesep,name,filesep);
mkdir(path,name);
save(nam_mat, 'Cn', 'pnr');
disp(sprintf('Saved as %s', nam_mat));

%% initialization of A, C
tic;
debug_on = false; %true; 
save_avi = false; 
neuron.options.min_corr = 0.85;  % min correlation
neuron.options.min_pnr = 10;  % min peak-to-noise ratio
patch_par = [2,2]; %1;  % divide the optical field into m X n patches and do initialization patch by patch
K = 300; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
neuron.options.bd = 1; % boundaries to be removed due to motion correction
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 

disp('Saving correlation image of initialized neuron...');
nam_mat = sprintf('%s%s%s%sf02-cn_after_init.mat',path,filesep,name,filesep);
save(nam_mat, 'Cn', 'center');
disp(sprintf('Saved as %s', nam_mat));

[~, srt] = sort(max(neuron.C, [], 2)./get_noise_fft(neuron.C), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();

%% merge neurons, order neurons and delete some low quality neurons (cell 0, before running iterative udpates)
% only parts implemented

neuron_bk = neuron.copy();
merge_thr = [0.1, 0.7, 0];     % thresholds for merging neurons corresponding to
%{sptial overlaps, temporal correlation of C, temporal correlation of S}
[merged_ROI, newIDs] = neuron.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
% A: spatial shapes; S: spike counts; C: calcium traces 

% sort neurons
[Cpnr, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
neuron.orderROIs(srt);
[Ain, Cin] = neuron.snapshot();   % keep the initialization results

disp('Saving contours of neurons...');
nam_mat = sprintf('%s%s%s%sf03-contours_after_init.mat',path,filesep,name,filesep);
save(nam_mat, 'Cn', 'neuron');
disp(sprintf('Saved as %s', nam_mat));


%% udpate background (cell 1, the following three blocks can be run iteratively)
% determine nonzero pixels for each neuron
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    sn = neuron.estNoise(Y);
else
    sn = neuron.P.sn; 
end
thresh = 5;     % threshold for detecting large cellular activity in each pixel. (mean + thresh*sn)
% start approximating theb background
tic;
clear Ysignal;
Ybg = Y-neuron.A*neuron.C;
ssub = 3;   % downsample the data to improve the speed
rr = neuron.options.gSiz*1;  % average neuron size, it will determine the neighbors for regressing each pixel's trace
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
Ybg = neuron.localBG(Ybg, ssub, rr, active_px, sn, 5); % estiamte local background.
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);

% subtract the background from the raw data.
Ysignal = Y - Ybg;
neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

for i=1:10
%% update spatial components (cell 2), we can iteratively run cell 2& 3 for few times and then run cell 1
% use HALS to update the spatial components
neuron.options.dist = 5;
IND = determine_search_location(neuron.A, 'ellipse', neuron.options);
neuron.A = HALS_spatial(Ysignal, neuron.A, neuron.C, IND, 10);
%neuron.post_process_spatial(); % uncomment this line to postprocess
ind = find(sum(neuron.A, 1)<=neuron.options.min_pixel);
neuron.delete(ind);
clear IND;
fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);

%% update C  (cell 3)
% update temporal components. 
tic;
smin = 5;       % thresholding the amplitude of the spike counts as smin*noise level
neuron.options.maxIter = 4;   % iterations to update C 
neuron.updateTemporal_endoscope(Ysignal, smin); 
fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);
end

%% udpate background (cell 1, the following three blocks can be run iteratively)
% determine nonzero pixels for each neuron
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    sn = neuron.estNoise(Y);
else
    sn = neuron.P.sn; 
end
thresh = 5;     % threshold for detecting large cellular activity in each pixel. (mean + thresh*sn)
% start approximating theb background
tic;
clear Ysignal;
Ybg = Y-neuron.A*neuron.C;
ssub = 3;   % downsample the data to improve the speed
rr = neuron.options.gSiz*1;  % average neuron size, it will determine the neighbors for regressing each pixel's trace
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
Ybg = neuron.localBG(Ybg, ssub, rr, active_px, sn, 5); % estiamte local background.
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);

% subtract the background from the raw data.
Ysignal = Y - Ybg;
neuron.playMovie(Ysignal); % play the video data after subtracting the background components.


%% pick neurons from the residual (cell 4). It's not always necessary
Yres = Ysignal - neuron.A*neuron.C;
neuron.options.min_corr = 0.9;
neuron.options.min_pnr = 10;
patch_par = [2, 2];
[center_new, Cn_res, pnr_res] = neuron.pickNeurons(Yres, patch_par, 'auto'); % method can be either 'auto' or 'manual'
% [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Yres, patch_par, 'manual'); % method can be either 'auto' or 'manual'

%% save results
result_nm = [dir_nm, file_nm, '_results.mat'];
neuron.save_results(result_nm); %save variable 'neuron' only.
% neuron.save_results(result_nm, neuron.reshape(Ybg, 2)); % save background as well



%% display neurons
dir_neurons = sprintf('%s%s_neurons%s', dir_nm, file_nm, filesep);
if exist('dir_neurons', 'dir')
    temp = cd();
    cd(dir_neurons);
    delete *;
    cd(temp);
else
    mkdir(dir_neurons);
end
neuron.viewNeurons([], neuron.C_raw, dir_neurons);
%% display contours of the neurons
figure;
% neuron.viewContours(Cn, 0.9, 0);  % correlation image computed with
% spatially filtered data
% [Cn, pnr] = neuron.correlation_pnr(Ysignal); % very slow 
Cnn = correlation_image(Ysignal(:, 1:5:end), 4, d1, d2);
neuron.viewContours(Cnn, 0.9, 0); % correlation image computed with background-subtracted data
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');
% plot contours with IDs
figure;
plot_contours(neuron.A, Cn, 0.9, 0, [], neuron.Coor);
