function [ ans ] = run_analysis( movieFiles, tag )
%RUN_ANALYSIS runs CNMF-E analysis on the data
%   set_parameters: a .m file containing all the parameters used by the
%   analysis algorithm.

% movieFiles={'/home/athina/Data/mc_recording_20170711_131010-1.tif'};
[inDir,nam,~] = fileparts(movieFiles{1});
recID = regexprep(nam, 'mc_recording_(\w*)-1','$1');
mat_nam = sprintf('%s.mat', recID);

% add paths
cnmfe_setup;
clear oasis_folder optimization_folder;

% import the parameter values
set_parameters = fullfile(CNMF_dir,'parameters',sprintf('parameters_cnmf_an%s',tag));
clear CNMF_dir;
run (set_parameters);
clear set_parameters;

% load the motion corrected movies
filename = fullfile(inDir, mat_nam);
if ~exist(filename,'file')
    tif2matMulti(movieFiles, filename);
end
data = matfile(filename);

Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

%% create Source2D class object for storing results and parameters
neuron_raw = Sources2D('d1',d1,'d2',d2,... % dimensions of datasets
    'ssub', ssub,...    % spatial downsampling factor
    'tsub', tsub, ...   %temporal downsampling factor
    'gSig', gSig,...    %width of the gaussian kernel, which can approximates the average neuron shape
    'gSiz', gSiz, ...   % maximum diameter of neurons in the image plane. larger values are preferred. 
    'dist', ddist, ...  % maximum size of the neuron: dist*gSiz
    'search_method', search_method, ... % searching method (neuron body vs. dendrites)
    'merge_thr', merge_thr, ... % threshold for merging neurons
    'bas_nonneg', bas_nonneg);  % 1: positive baseline of each calcium traces; 0: any baseline

neuron_raw.Fs = fs;         % frame rate
neuron_raw.kernel = kernel; % convolution kernel

clear fs gSig gSiz ddist search_method bas_nonneg merge_thr;
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
    [Y, neuron] = neuron_raw.load_data(filename, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
end
Y = neuron.reshape(Y, 1);
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);
clear sframe num2read d1s d2s T data neuron_raw;

%% compute correlation image and peak-to-noise ratio image.
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
fprintf('Saved as %s\n', nam_mat);

%% options for running deconvolution 
neuron.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ... % optimize the baseline
    'optimize_smin', false);  % optimize the threshold 

%% initialization of A, C
tic;

neuron.options.nk = nk; % number of knots for creating spline basis
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
clear nk min_corr min_pnr bd min_pixel;

[center, Cn, ~] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 

disp('Saving correlation image of initialized neuron...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f02-cn_after_init.mat');
save(nam_mat, 'Cn', 'center');
clear center;
fprintf('Saved as %s', nam_mat);

% Terminate the script if no neurons were detected
neurDetected = size(neuron.A, 2);
if neurDetected == 0
    fprintf('\nATTENTION: No neurons detected! Termination.\n');
    exit(1);
end
fprintf('\n%d neurons detected.\n', neurDetected);
clear neurDetected;

[~, srt] = sort(max(neuron.C, [], 2), 'descend');
neuron.orderROIs(srt);

%% iteratively update A, C and B

% parameters, estimate the background
if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
    %sn = neuron.estNoise(Y);
    neuron.P.sn = neuron.estNoise(Y);
%else
    %sn = neuron.P.sn; 
end

neuron.options.maxIter = maxIter_temporal;   % iterations to update C

% parameters for running iterations
nC = size(neuron.C, 1);    % number of neurons 

maxIter = 5;        % maximum number of iterations 
for miter=1:maxIter
    %% merge neurons, order neurons and delete some low quality neurons
    
    % merge neurons
    neuron.quickMerge(merge_thr_all);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
    % A: spatial shapes; S: spike counts; C: calcium traces

    % sort neurons
    [~, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
    neuron.orderROIs(srt);
    
    %% update background (cell 1, the following three blocks can be run iteratively)
    % estimate the background
    clear Ysignal;
    tic; 
    Ybg = Y-neuron.A*neuron.C;
    rr = ceil(neuron.options.gSiz * bg_neuron_ratio); 
    active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
    [Ybg, ~] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, neuron.P.sn, 5); % estiamte local background.

    % subtract the background from the raw data.
    Ysignal = Y - Ybg;
    fprintf('Time cost in estimating the background: %.2f seconds\n', toc);
    
    %% update spatial & temporal components
    tic;
    for m=1:5    
        %temporal
        neuron.updateTemporal_endoscope(Ysignal, smin);
        [merged_ROI, ~] = neuron.quickMerge(merge_thr_all); 
        
        % sort neurons
        [~, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
        neuron.orderROIs(srt);
        
        
        %spatial
        neuron.updateSpatial_endoscope(Ysignal, maxIter_spatial);
        if isempty(merged_ROI)
            break;
        end
    end
    fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
    
    %% pick neurons from the residual (cell 4).
    if miter==1
        neuron.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
        neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, 'auto'); % method can be either 'auto' or 'manual'
    end
    
    %% stop the iteration 
    temp = size(neuron.C, 1); 
    if nC~=temp 
        break;
    end
    
    nC = temp;
end

%% apply results to the full resolution
% not supported (is this needed?)

%% save results
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'results.mat');
neuron.save_results(nam_mat, '-v7.3'); %save variable 'neuron' only.
% neuron.save_results(result_nm, neuron.reshape(Ybg, 2)); % save background as well
clear nam_mat;

%% save neurons for display
dir_neurons = sprintf('%s%s%s-%s%sneurons%s', path,filesep,name,tag,filesep,filesep);
referenceImg = mean(reshape(Y(:,1:100),neuron.options.d1,neuron.options.d2,100),3);
disp('Saving neuron and neurons dir...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f06-neurons.mat');
save(nam_mat, 'dir_neurons', 'neuron', 'referenceImg', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));
clear dir_neurons;

%% save neural contours for display
disp('Saving contours of neurons for display...');
nam_mat = fullfile(path,sprintf('%s-%s',name,tag),'f07-contours.mat');
save(nam_mat, 'neuron', 'Ysignal', 'd1', 'd2', 'Cn', '-v7.3');
disp(sprintf('Saved as %s', nam_mat));

clear Ysignal d1 d2 Cn;

%% Save data into text files
saveDir = fullfile(path,sprintf('%s-%s',name,tag));
footprints = full(neuron.A);
csvwrite(fullfile(saveDir, 'footprints.txt'), footprints);

traces = neuron.C;
csvwrite(fullfile(saveDir, 'traces.txt'), traces);

spikes = neuron.S;
csvwrite(fullfile(saveDir, 'spikes.txt'), spikes);

% Save neuron object
data = neuron;
save(fullfile(saveDir, 'data.mat'), 'neuron');
clear saveDir;


