function [ neuron ] = remove_bad_cells( data, out_dir )
%REMOVE_HOLLOW_CELLS Performs tests on the neurons' spatial footprints and
%removes those that have bad spatial footprints.
% Input:
%   data    :   A neuron object (class defined by CNMF-E)
%   out_dir :   Output directory. The new neuron object is saved there.
%   Depending on the options, also figures with the statistics per neuron
%   and other files will be saved in this directory. If not defined,
%   the new data will be saved in the same directory with the data.
% Output: 
%   neuron  :   The updated neuron object (bad cells removed).

% data = '/home/athina/Data/20170724_094625-014/data.mat'

% determine the output directory if not defined
if ~ exist('out_dir','var')
    out_dir = fileparts(data);
end

% load data
run(which('cnmfe_setup'));  % we need to initialize cnmfe to load the neuron object
load(data);                 % contains variable named neuron with all data

% get parameters from script
run('parameters');

% create directory to save neurons (if option is selected)
if saveFigures
    fig_out_dir = fullfile(out_dir, 'bc_stats');
    if ~ exist('fig_out_dir', 'dir'), mkdir(fig_out_dir); end
end

% retrieve information about the recording
[d1, d2] = size(neuron.Cn);     % number of pixels in the y and x axis
nframes = size(neuron.C, 2);    % number of frames
nNeurons = size(neuron.C,1);    % number of neurons detected by CNMF-E

% initialize vectors to save statistics per frame
mserrors = zeros(1,nNeurons);   % mean squared error when compared to an ideal shape
notgaus = zeros(1,nNeurons);    % gaussian hypothesis rejection
numHoles = zeros(1,nNeurons);   % number of wholes inside the spatial footprint
badCell = false(1,nNeurons);    % saves which neurons are badly shaped and will be deleted

if or(dispFig, saveFigures), doPlot = true; else doPlot = false; end

if doPlot
    if dispFig
        h = figure; % create a figure that will host all the subplots
    else
        h = figure('visible','off');
    end
    
    % retrieve information about time
    [t, str_xlabel] = calculate_timeline(nframes,neuron.Fs);
end

for i=1:nNeurons % for each neuron
    % variables to store test results
    mse_pass = true; gaus_pass = true; holes_pass = true; gen_pass = true;
    
    % get the spatial footprints in 2-dimensional form
    neur2d = reshape(neuron.A(:,i), size(neuron.Cn));
    
    % calculate statistics for this spatial footprint
    if dispFig, clf; 
        subplot('421'); plot_cropped_neuron(neur2d, i);
        subplot('422'); 
    end
    numHoles(i) = count_holes(neur2d, dispFig);
    if dispFig, subplot('423'); end
    mserrors(i) = ideal_comparison(neur2d, dispFig);
    if dispFig, subplot('413'); end
    notgaus(i) = gaussian_test(neuron.A(:,i), dispFig);
    if dispFig, subplot('414'); 
    plot_calcium(neuron.C_raw(i,:), neuron.C(i,:), max(neuron.A(:,i)), t);
    xlabel(str_xlabel);
    end
    
    % compare all statistics with the user defined limits
    if numHoles(i) > ht, holes_pass = false; gen_pass = false; end
    if mserrors(i) > elim, mse_pass = false; gen_pass = false; end
    if rejng
        if notgaus(i) == 1, gaus_pass = false; gen_pass = false; end
    else
        gaus_pass = NaN;
    end
    
    % show a summary of all test results
    if dispFig
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
    
    if dispFig
        pause(1);
    end
    
    
    if saveFigures
        if gen_pass, qual = 'good'; else qual = 'bad'; end
        figName = fullfile(fig_out_dir,sprintf('neuron_%04d_%s.png', i, qual));
        saveas(h,figName,'png')
    end
end

%% delete badly shaped cells
neuron.A = neuron.A(:,~badCell);
neuron.C = neuron.C(~badCell,:);
neuron.S = neuron.S(~badCell,:);
neuron.C_raw = neuron.C_raw(~badCell,:);

%% save results
if saveUpdNeuron
    % Save neuron object
    save(fullfile(out_dir, 'updated_data.mat'), 'neuron', '-v7.3');
end

if saveForCellReg
    % save spatial footprints as an object (compatible with CellReg)
    allFiltersMat = reshape(neuron.A,[],size(neuron.Cn,1),size(neuron.Cn,2));
    save(fullfile(out_dir, 'spatial_footprints.mat'), 'allFiltersMat', '-v7.3');
end

if saveAsText
    % save footprints, traces and spikes as text files for analysis
    csvwrite(fullfile(out_dir, 'footprints.txt'), neuron.A);
    csvwrite(fullfile(out_dir, 'traces.txt'), neuron.C);
    csvwrite(fullfile(out_dir, 'spikes.txt'), neuron.S);
end

end
