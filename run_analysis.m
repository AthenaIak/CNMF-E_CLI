function [ ans ] = run_analysis( set_parameters )
%RUN_ANALYSIS runs CNMF-E analysis on the data
%   filename: the path+name of the file containing the data to be analyzed.
%   crop :  crops the border of the data. Array with 4 values (top, right,
%   bottom left). Default=[0 0 0 0].
%   patch_par: defines how the raw data will be split. E.g. [2 2] splits
%   the image to 2x2 images and then performs calculations. Default=1. If
%   you can avoid splitting the image (enough RAM), avoid it.

% add paths
run_setup;
clear CNMF_dir;

% set_parameters='~/tu/athina/Data/analyzed/parameters/parameters_an006_cnmf';
run (set_parameters);

% load the data
%data = loadRawData(filename);
[~,~,ext] = fileparts(filename);

if strcmp(ext,'.tif')
    nam_mat = tif2mat(filename);
elseif strcmp(ext,'.mat')
    nam_mat = filename;
else
    disp('Unsupported format: only .tif and .mat files allowed');
    exit(1);
end
data = matfile(nam_mat);
clear ext num_mat;

Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames
%Yfs = data.Yfs;

%% create Source2D class object for storing results and parameters
neuron_raw = Sources2D('d1',d1,'d2',d2);   % dimensions of datasets
neuron_raw.Fs = fs;         % frame rate
neuron_raw.updateParams('ssub', ssub,...  % spatial downsampling factor
    'tsub', tsub, ...  %temporal downsampling factor
    'gSig', gSig,... %width of the gaussian kernel, which can approximates the average neuron shape
    'gSiz', gSiz, ...% maximum diameter of neurons in the image plane. larger values are preferred. 
    'dist', ddist, ... % maximum size of the neuron: dist*gSiz
    'search_method', search_method, ... % searching method
    'merge_thr', merge_thr, ... % threshold for merging neurons
    'bas_nonneg', bas_nonneg);   % 1: positive baseline of each calcium traces; 0: any baseline

neuron_raw.kernel = kernel;

clear fs gSig gSiz ddist search_method merge_thr bas_nonneg;
clear ssub tsub tau_decay tau_rise nframe_decay bound_pars kernel;

%% downsample data for fast and better initialization
sframe=1;           % user input: first frame to read (optional, default:1)
num2read= numFrame; % user input: how many frames to read   (optional, default: until the end)

tic;
if and(neuron_raw.options.ssub==1, neuron_raw.options.tsub==1)
    neuron = neuron_raw;
    Y = double(data.Y(:, :, sframe+(1:num2read)-1));
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else
    [Y, neuron] = neuron_raw.load_data(nam_mat, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
end
Y = neuron.reshape(Y, 1);
%neuron_raw.P.p = 2;      %order of AR model

fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);
clear sframe num2read d1s d2s T data neuron_raw;

% compute correlation image and peak-to-noise ratio image.
% this step is not necessary, but it can give you some ideas of how your
% data look like

neuron.options.nb = nb; % Number of background elements (default = 1)
[Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, numFrame, 1000)))); % calls correlation_image_endoscope
% save to be viewed where GUI is available
[path,name,~] = fileparts(filename);
output_dir = sprintf('%s%s%s-%s',path,filesep,name,tag);
if exist(output_dir, 'dir')
    temp = cd();
    cd(output_dir);
    delete *;
    cd(temp);
else
    mkdir(output_dir);
end
clear nb output_dir;
disp('Saving correlation image and peak-to-noise ratio image...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f01-cn&pnr.mat');
save(nam_mat, 'Cn', 'pnr', '-v7.3'); % specify version 7.3 to allow partial loading
clear Cn pnr;
disp(sprintf('Saved as %s', nam_mat));

%% initialization of A, C
tic;

neuron.options.nk = nk; % number of knots for creating spline basis
neuron.options.min_corr = min_corr;  % min correlation (default = 0.3)
neuron.options.min_pnr = min_pnr;  % min peak-to-noise ratio % (default = 10)
neuron.options.bd = bd; % boundaries to be removed due to motion correction
clear nk min_corr min_pnr bd;
[center, Cn, ~] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 

disp('Saving correlation image of initialized neuron...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f02-cn_after_init.mat');
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
%[merged_ROI, newIDs] = neuron.quickMerge(merge_thr_after);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
neuron.quickMerge(merge_thr_after); % merge neurons based on the correlation computed with {'A', 'S', 'C'}
% A: spatial shapes; S: spike counts; C: calcium traces 
clear merge_thr_after;

% sort neurons
%[Cpnr, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
[~, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
neuron.orderROIs(srt);
clear srt;
%[Ain, Cin] = neuron.snapshot();   % keep the initialization results

disp('Saving contours of neurons...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f03-contours_after_init.mat');
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
%disp('Saving video data after subtracting the background...');
%nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f04-Ysignal_background_subtracted.mat');
%save(nam_mat, 'Ysignal', 'neuron', '-v7.3');
%disp(sprintf('Saved as %s', nam_mat));

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
%disp('Saving video data after subtracting the background (end)...');
%nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f05-Ysignal_end.mat');
%save(nam_mat, 'Ysignal', 'neuron', '-v7.3');
%disp(sprintf('Saved as %s', nam_mat));

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
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'results.mat');
neuron.save_results(nam_mat, '-v7.3'); %save variable 'neuron' only.
% neuron.save_results(result_nm, neuron.reshape(Ybg, 2)); % save background as well
clear nam_mat;

%% save neurons for display
dir_neurons = sprintf('%s%s%s-%s%sneurons%s', path,filesep,name,tag,filesep,filesep);
disp('Saving neuron and neurons dir...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f06-neurons.mat');
save(nam_mat, 'dir_neurons', 'neuron', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));
clear dir_neurons;

%% save neural contours for display
disp('Saving contours of neurons for display...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f07-contours.mat');
save(nam_mat, 'neuron', 'Ysignal', 'd1', 'd2', 'Cn', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));

clear Ysignal neuron d1 d2 Cn;

