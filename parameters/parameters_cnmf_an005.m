
%crop = 20;
fs = 5; % frequency of the movie (frames per second)


% preprocessing and neuron expectations
patch_par=[2 2];% divide the optical field into m X n patches and do initialization patch by patch
ssub = 2; 	% spatial downsampling factor
tsub = 1; 	% temporal downsampling factor
gSig = 5; 	% width of the gaussian kernel, which can approximate the average neuron shape (default = 4)
gSiz = 25; 	% maximum diameter of neurons in the image plane. larger values are preferred. 
ddist = 2; 	% maximum size of the neuron: dist*gSiz
search_method = 'ellipse';	% searching method
merge_thr=0.7; 	% threshold for merging neurons (the higher, the more seed pixels are detected and not merged into 1nu) (default = 0.7)
bas_nonneg = 1; % 1: positive baseline of each calcium traces; 0: any baseline
nb = 1; 	% Number of background elements (default = 1)

% convolution kernel
tau_decay = 1;  % unit: second 
tau_rise = 0.1; 
nframe_decay = ceil(6*tau_decay*fs); 
bound_pars = false; 
kernel = create_kernel('exp2', [tau_decay, tau_rise]*fs, nframe_decay, [], [], bound_pars); 

% options for running deconvolution 
dectype = 'ar1'; 	% model of the calcium traces. {'ar1', 'ar2'}
decmethod = 'thresholded'; % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
optimize_pars = true;  	% optimize AR coefficients
optimize_b = true; 	% optimize the baseline
optimize_smin = false;  % optimize the threshold 

% options used during/after initialization of A, C
debug_on = false;
save_avi = false; 
K = 150; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
min_corr = 0.3; % Minimum correlation for separating neurons (default = 0.9)
min_pnr = 20;  % min peak-to-noise ratio (default = 20)
min_pixel = 15; % minimum number of (non-negative) pixels that describe the neuron
nk = 5; % number of knots for creating spline basis
bd = 20; % boundaries to be removed due to motion correction (does not work properly)

% merge options (used after initialization)
merge_thr_all = [1e-5, 0.70, .1];     % thresholds for merging neurons
% corresponding to {spatial overlaps, temporal correlation of C, temporal correlation of S}

% UPDATING
% background estimation
spatial_ds_factor = 3;      % spatial downsampling factor. it's for faster estimation when updating background
bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
maxIter_spatial = 5;       % number of iterations required

% parameters, estimate the temporal components
maxIter_temporal = 4;
smin = 5;       % thresholding the amplitude of the spike counts as smin*noise level

