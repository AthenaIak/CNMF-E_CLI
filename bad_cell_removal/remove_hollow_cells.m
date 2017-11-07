function [ updated_neuron ] = remove_hollow_cells( data, out_dir )
%REMOVE_HOLLOW_CELLS Performs tests on the neurons' spatial footprints and
%removes those that have bad spatial footprints.
% Input:
%   data    :   A neuron object (class defined by CNMF-E)
%   out_dir :   Output directory. The new neuron object is saved there.
%   Depending on the options, also figures with the statistics per neuron
%   and other files will be saved in this directory. If not defined,
%   the new data will be saved in the same directory with the data.

% data = '/home/athina/Data/20170724_094625-014/data.mat'

% determine the output directory if not defined
if ~ exist('out_dir','var')
    path = fileparts(data);
end

% load data
run(which('cnmfe_setup'));  % we need to initialize cnmfe to load the neuron object
load(data);

% get parameters from script
parameters;

% retrieve information about the recording
[d1, d2] = size(neuron.Cn);     % number of pixels in the y and x axis
nframes = size(neuron.C, 2);    % number of frames
nNeurons = size(neuron.C,1);    % number of neurons detected by CNMF-E

% initialize vectors to save statistics per frame
mserrors = zeros(1,nNeurons);   % mean squared error when compared to an ideal shape
notgaus = zeros(1,nNeurons);    % gaussian hypothesis rejection
numHoles = zeros(1,nNeurons);   % number of wholes inside the spatial footprint
badCell = false(1,nNeurons);    % saves which neurons are badly shaped and will be deleted

if display
    figure; % create a figure that will host all the subplots
    
    % retrieve information about time
    [t, str_xlabel] = calculate_timeline(nframes,neuron.Fs);
end

for i=1:nNeurons % for each neuron
    % variables to store test results
    mse_pass = true; gaus_pass = true; holes_pass = true; gen_pass = true;
    
    % get the spatial footprints in 2-dimensional form
    neur2d = reshape(neuron.A(:,i), size(neuron.Cn));
    
    % calculate statistics for this spatial footprint
    clf; % clear everything drawn before
    if display, subplot('421'); plot_cropped_neuron(neur2d, i); end
    if display, subplot('422'); end
    numHoles(i) = count_holes(neur2d, display);
    if display, subplot('423'); end
    mserrors(i) = ideal_comparison(neur2d, display);
    if display, subplot('413'); end
    notgaus(i) = gaussian_test(neuron.A(:,i), display);
    if display, subplot('414'); 
    plot_calcium(neuron.C_raw(i,:), neuron.C(i,:), max(neuron.A(:,i)), t);
    xlabel(str_xlabel);
    end
    
    % compare all statistics with the user defined limits
    if numHoles(i) > ht, holes_pass = false; gen_pass = false; end
    if mserrors(i) > elim, mse_pass = false; gen_pass = false; end
    if rng
        if notgaus(i) == 1, gaus_pass = false; gen_pass = false; end
    else
        gaus_pass = NaN;
    end
    
    % show a summary of all test results
    if display
        subplot('424');
        text(0,0.5,sprintf('Summary\nHoles: \nShape: \nGaussian: ')); axis off
        if holes_pass, text(0.65,0.56,sprintf('Pass'),'color','green');
        else text(0.65,0.56,sprintf('Fail'),'color','red'); end
        if mse_pass text(0.65,0.45,sprintf('Pass'),'color','green');
        else text(0.65,0.45,sprintf('Fail'),'color','red'); end
        if isnan(gaus_pass) text(0.65,0.33,sprintf('---'),'color','blue');
        elseif gaus_pass text(0.65,0.33,sprintf('Pass'),'color','green');
        else text(0.65,0.33,sprintf('Fail'),'color','red'); end
    end
    
    % mark the neuron if it doesn't pass one of the tests
    if ~ gen_pass, badCell(i) = true; end
    
    %noiseLevel = mean(neur1d(find(neur1d>0)))-3*std(neur1d(find(neur1d>0)))
    
    if display
        pause(1);
    end
end

% create new, updated neuron object
updated_neuron = neuron.copy();
updated_neuron.A = neuron.A(:,~badCell);
updated_neuron.C = neuron.C(~badCell,:);
updated_neuron.S = neuron.S(~badCell,:);
updated_neuron.C_raw = neuron.C_raw(~badCell,:);




end
