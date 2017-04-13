function [neuron_raw] = neuronForData(d1, d2)
%NEURONFORDATA initializes a neuron object from raw data
%	d1	The number of rows of each frame
%	d2	The number of columns of each frame

%% create Source2D class object for storing results and parameters
neuron_raw = Sources2D('d1',d1,'d2',d2);   % dimensions of datasets
neuron_raw.Fs = 5;         % frame rate
ssub = 1;           % spatial downsampling factor
tsub = 1;           % temporal downsampling factor
neuron_raw.updateParams('ssub', ssub,...  % spatial downsampling factor
    'tsub', tsub, ...  %temporal downsampling factor
    'gSig', 4,... %width of the gaussian kernel, which can approximates the average neuron shape
    'gSiz', 15, ...% maximum diameter of neurons in the image plane. larger values are preferred. 
    'dist', 2, ... % maximum size of the neuron: dist*gSiz
    'search_method', 'ellipse', ... % searching method
    'merge_thr', 0.7, ... % threshold for merging neurons
    'bas_nonneg', 1);   % 1: positive baseline of each calcium traces; 0: any baseline

% create convolution kernel to model the shape of calcium transients 
tau_decay = 1;  % unit: second 
tau_rise = 0.1; 
nframe_decay = ceil(6*tau_decay*neuron_raw.Fs); 

bound_pars = false; 
neuron_raw.kernel = create_kernel('exp2', [tau_decay, tau_rise]*neuron_raw.Fs, nframe_decay, [], [], bound_pars); 
