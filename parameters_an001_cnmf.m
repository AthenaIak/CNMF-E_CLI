tag = '001'; 	% analysis unique identifier

% define the .mat file to be analysed
filename = '~/tp/athina/data/4294/output/20170626_of_20160531/ready_20160531.mat';
fs = 5; % frequency of the movie (frames per second)


% preprocessing and neuron expectations
patch_par=1; 	% divide the optical field into m X n patches and do initialization patch by patch
ssub = 1; 	% spatial downsampling factor
tsub = 1; 	% temporal downsampling factor
gSig = 14; 	% width of the gaussian kernel, which can approximate the average neuron shape (default = 4)
gSiz = 49; 	% maximum diameter of neurons in the image plane. larger values are preferred. 
dist = 2; 	% maximum size of the neuron: dist*gSiz
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

% options used during/after initialization of A, C
debug_on = false;
save_avi = false; 
K = 300; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
min_corr = 0.3; % Minimum correlation for separating neurons (default = 0.3)
min_pnr = 10;  % min peak-to-noise ratio (default = 10)
nk = 5; % number of knots for creating spline basis
bd = 1; % boundaries to be removed due to motion correction

% merge options (used after initialization)
merge_thr_after = [0.1, 0.7, 0]; % thresholds for merging neurons corresponding to
%{sptial overlaps, temporal correlation of C, temporal correlation of S}




