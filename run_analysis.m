function [ ans ] = run_analysis( filename, crop, patch_par )
%RUN_ANALYSIS runs CNMF-E analysis on the data
%   filename: the path+name of the file containing the data to be analyzed.
%   crop :  crops the border of the data. Array with 4 values (top, right,
%   bottom left). Default=[0 0 0 0].
%   patch_par: defines how the raw data will be split. E.g. [2 2] splits
%   the image to 2x2 images and then performs calculations. Default=1. If
%   you can avoid splitting the image (enough RAM), avoid it.

if ~exist('crop', 'var'), crop = [0, 0 0 0]; end
if ~exist('patch_par', 'var'), patch_par = 1; end

% add paths
run_setup;
clear CNMF_dir;

% filename = 'D:\�����������\CNMF_E\demos\data_endoscope.tif'
% filename = '~/tp/athina/data/4294/output/ready_recording_20160531_144143.mat';
% filename = '~/tp/athina/data/4294/output/small_recording_20160531_144143.mat'
% filename = '~/tp/data/iHPC5 raw/recording_20160125_114832/output/mcorr_mosaic-recording_20160125_114832.tif';
% filename = '~/tp/data/iHPC5 raw/recording_20160125_114832/output/mcorr_mosaic_128-recording_20160125_114832.tif';
% filename = '~/tp/data/iHPC5 raw/recording_20160125_114832/output/mcorr_moco-recording_20160125_114832.tif';
% filename = '~/tp/data/iHPC5 raw/recording_20160125_114832/output/mcorr_moco_128-recording_20160125_114832.tif';
% filename = '~/Data/iHPC5 raw/output/mcorr_128-2_recording_20160118_140521.tif'
% filename = '/home/athina/Data/demo_endoscope/data_endoscope.tif'
% filename = '~/Data/iHPC5 raw/output/mcorr_128_recording_20160118_140521.tif';
% filename = 'D:\�����������\data\mcorr_128_recording_20160118_140521.tif'
data = loadRawData(filename);

%% cut the borders of the image to reduce size (should not normally do that)
Y = data.Y(crop(1)+1:end-crop(3),crop(4)+1:end-crop(2),:);
Yfs = data.Yfs;
clear data;
Ysiz = size(Y)';
filename = '~/tp/athina/data/4294/output/gSig56.mat';
save(filename, 'Y', 'Ysiz','Yfs','-v7.3');
clear Y Ysiz Yfs crop;
%
data = loadRawData(filename);
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames
%Yfs = data.Yfs;

% create container structure neuron
neuron_raw = neuronForData(d1,d2);
[Y, neuron, ~] = downsample(neuron_raw, data);
clear data neuron_raw;

% compute correlation image and peak-to-noise ratio image.
% this step is not necessary, but it can give you some ideas of how your
% data look like
neuron.options.gSiz = 49; % Average size of neuron (default = 15)
%neuron.options.nb = 1; % Number of background elements (default = 1)
%neuron.options.min_corr = 0.3; % Minimum correlation for separating neurons (default = 0.3)
[Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, numFrame, 1000)))); % calls correlation_image_endoscope
% save to be viewed where GUI is available
[path,name,~] = fileparts(filename);
output_dir = sprintf('%s%s%s',path,filesep,name);
if exist(output_dir, 'dir')
    temp = cd();
    cd(output_dir);
    delete *;
    cd(temp);
else
    mkdir(output_dir);
end
clear output_dir;
disp('Saving correlation image and peak-to-noise ratio image...');
nam_mat = sprintf('%s%s%s%sf01-cn&pnr.mat',path,filesep,name,filesep);
save(nam_mat, 'Cn', 'pnr', '-v7.3'); % specify version 7.3 to allow partial loading
clear Cn pnr;
disp(sprintf('Saved as %s', nam_mat));

%% initialization of A, C
tic;
debug_on = false; %true; 
save_avi = false; 
%neuron.options.nk = 5; % number of knots for creating spline basis
neuron.options.min_corr = 0.3;  % min correlation (default = 0.3)
neuron.options.min_pnr = 10;  % min peak-to-noise ratio % (default = 10)
neuron.options.merge_thr = .7; % merge threshold (the higher, the more seed pixels are detected and not merged into 1nu) (default = 0.7)
neuron.options.gSig = 14; % width of the gaussian used for spatial filtering (default = 4)
%patch_par = [2 2]; %1;  % divide the optical field into m X n patches and do initialization patch by patch
K = 300; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
neuron.options.bd = 1; % boundaries to be removed due to motion correction
[center, Cn, ~] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 

disp('Saving correlation image of initialized neuron...');
nam_mat = sprintf('%s%s%s%sf02-cn_after_init.mat',path,filesep,name,filesep);
save(nam_mat, 'Cn', 'center');
clear center;
disp(sprintf('Saved as %s', nam_mat));

% Terminate the script if no neurons were detected
neurDetected = size(neuron.A, 2);
if neurDetected == 0
    disp(sprintf('\nATTENTION: No neurons detected! Termination.\n'));
    exit(1);
end
disp(sprintf('\n%d neurons detected.\n', neurDetected));
clear neurDetected;

[~, srt] = sort(max(neuron.C, [], 2)./get_noise_fft(neuron.C), 'descend'); % can crush if no neurons detected
neuron.orderROIs(srt);
%neuron_init = neuron.copy();

%% merge neurons, order neurons and delete some low quality neurons (cell 0, before running iterative udpates)
% only parts implemented

%neuron_bk = neuron.copy();
merge_thr = [0.1, 0.7, 0];     % thresholds for merging neurons corresponding to
%{sptial overlaps, temporal correlation of C, temporal correlation of S}
%[merged_ROI, newIDs] = neuron.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
neuron.quickMerge(merge_thr); % merge neurons based on the correlation computed with {'A', 'S', 'C'}
% A: spatial shapes; S: spike counts; C: calcium traces 
clear merge_thr;

% sort neurons
%[Cpnr, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
[~, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
neuron.orderROIs(srt);
clear srt;
%[Ain, Cin] = neuron.snapshot();   % keep the initialization results

disp('Saving contours of neurons...');
nam_mat = sprintf('%s%s%s%sf03-contours_after_init.mat',path,filesep,name,filesep);
save(nam_mat, 'Cn', 'neuron','-v7.3');
disp(sprintf('Saved as %s', nam_mat));


%% update background (cell 1, the following three blocks can be run iteratively)
% determine nonzero pixels for each neuron
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    sn = neuron.estNoise(Y);
else
    sn = neuron.P.sn; 
end
%thresh = 5;     % threshold for detecting large cellular activity in each pixel. (mean + thresh*sn)
% start approximating theb background
tic;
Ybg = Y-neuron.A*neuron.C;
ssub = 3;   % downsample the data to improve the speed
rr = neuron.options.gSiz*1;  % average neuron size, it will determine the neighbors for regressing each pixel's trace
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
Ybg = neuron.localBG(Ybg, ssub, rr, active_px, sn, 5); % estiamte local background.
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);

% subtract the background from the raw data.
Ysignal = Y - Ybg;
clear Ybg ssub  rr active_px;
disp('Saving video data after subtracting the background...');
nam_mat = sprintf('%s%s%s%sf04-Ysignal_background_subtracted.mat',path,filesep,name,filesep);
save(nam_mat, 'Ysignal', 'neuron', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% HERE
for i=1:10
%% update spatial components (cell 2), we can iteratively run cell 2& 3 for few times and then run cell 1
% use HALS to update the spatial components
neuron.options.dist = 5;
IND = determine_search_location(neuron.A, 'ellipse', neuron.options);
neuron.A = HALS_spatial(Ysignal, neuron.A, neuron.C, IND, 10);
%neuron.post_process_spatial(); % uncomment this line to postprocess
ind = find(sum(neuron.A, 1)<=neuron.options.min_pixel);
neuron.delete(ind);
clear IND ind;
fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);

%% update C  (cell 3)
% update temporal components. 
tic;
smin = 5;       % thresholding the amplitude of the spike counts as smin*noise level
neuron.options.maxIter = 4;   % iterations to update C 
neuron.updateTemporal_endoscope(Ysignal, smin); 

fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);
end
clear Ysignal smin i;

%% update background (cell 1 again, after cell 2&3 are run iteratively)
% determine nonzero pixels for each neuron
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    sn = neuron.estNoise(Y);
else
    sn = neuron.P.sn; 
end
%thresh = 5;     % threshold for detecting large cellular activity in each pixel. (mean + thresh*sn)
% start approximating theb background
tic;
Ybg = Y-neuron.A*neuron.C;
ssub = 3;   % downsample the data to improve the speed
rr = neuron.options.gSiz*1;  % average neuron size, it will determine the neighbors for regressing each pixel's trace
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
Ybg = neuron.localBG(Ybg, ssub, rr, active_px, sn, 5); % estiamte local background.
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);

% subtract the background from the rawls data.
Ysignal = Y - Ybg;
clear Ybg ssub rr active_px sn;
disp('Saving video data after subtracting the background (end)...');
nam_mat = sprintf('%s%s%s%sf05-Ysignal_end.mat',path,filesep,name,filesep);
save(nam_mat, 'Ysignal', 'neuron', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));

%% pick neurons from the residual (cell 4). It's not always necessary
%{
% let's not do this for now. it doesn't seem to do anything + consumes too
much ram
Yres = Ysignal - neuron.A*neuron.C;

neuron.options.min_corr = 0.9;
neuron.options.min_pnr = 10;
%patch_par = [2, 2]; % 1;
neuron.pickNeurons(Yres, patch_par, 'auto');
%[center_new, Cn_res, pnr_res] = neuron.pickNeurons(Yres, patch_par, 'auto'); % method can be either 'auto' or 'manual'
% [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Yres, patch_par, 'manual'); % method can be either 'auto' or 'manual'
clear Yres;
%}

%% save results
result_nm = [path, 'results.mat'];
neuron.save_results(result_nm, '-v7.3'); %save variable 'neuron' only.
% neuron.save_results(result_nm, neuron.reshape(Ybg, 2)); % save background as well
clear result_nm;

%% save neurons for display
dir_neurons = sprintf('%s%s%s%sneurons%s', path,filesep,name,filesep,filesep);
disp('Saving neuron and neurons dir...');
nam_mat = sprintf('%s%s%s%sf06-neurons.mat',path,filesep,name,filesep);
save(nam_mat, 'dir_neurons', 'neuron', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));
clear dir_neurons;

%% save neural contours for display
disp('Saving contours of neurons for display...');
nam_mat = sprintf('%s%s%s%sf07-contours.mat',path,filesep,name,filesep);
save(nam_mat, 'neuron', 'Ysignal', 'd1', 'd2', 'Cn', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));

clear Ysignal neuron d1 d2 Cn;

