function [ neuron ] = bc_visualization( data, out_dir )
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

%% Initialization
% load data
run(which('cnmfe_setup'));  % we need to initialize cnmfe to load the neuron object
load(data);                 % contains variable named neuron with all data

% get parameters from script
run('parameters');

% determine the output directory if not defined
if ~ exist('out_dir','var')
    out_dir = fileparts(data);
end

% create directory to save neurons (if option is selected)
if saveFigures
    fig_out_dir = fullfile(out_dir, 'bc_stats');
    if ~ exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end
end

% retrieve information about the recording
[d1, d2] = size(neuron.Cn);     % number of pixels in the y and x axis
nframes = size(neuron.C, 2);    % number of frames
nNeurons = size(neuron.C,1);    % number of neurons detected by CNMF-E

% initialize vectors to save statistics per frame
mserrors = zeros(1,nNeurons);   % mean squared error when compared to an ideal shape
similarity = zeros(1,nNeurons); 
notgaus = zeros(1,nNeurons);    % gaussian hypothesis rejection
numHoles = zeros(1,nNeurons);   % number of holes inside the spatial footprint
numLocmins = zeros(1,nNeurons); % number of local minimums inside the spatial footprint
badCell = false(1,nNeurons);    % saves which neurons are badly shaped and will be deleted

% create plots if the user chose to display or save figures
if or(dispFig, saveFigures), doPlot = true; else, doPlot = false; end

% initialize plotting
if doPlot
    if dispFig
        h = figure; % create a figure that will host all the subplots
    else
        h = figure('visible','off');
    end
    
    % retrieve information about time
    [t, str_xlabel] = calculate_timeline(nframes,neuron.Fs);
end

%% perform tests to all detected neurons
for i=1:nNeurons % for each neuron
    if mod(i,100)== 0
        fprintf('Neuron %d of %d',i, nNeurons);
    end
    
    % variables to store test results
    holes_pass = true; locmin_pass = true; mse_pass = true; 
    sim_pass = true; gaus_pass = true; gen_pass = true;
    
    % get the spatial footprints in 2-dimensional form
    neur2d = reshape(neuron.A(:,i), size(neuron.Cn));
    [x,y] = ind2sub(size(neur2d),find(neur2d>0)); % find the indexes to ...
    neur2d = neur2d(min(x):max(x),min(y):max(y)); % ... crop the boundaries
    
    % calculate statistics for this spatial footprint 
    % |--- plot neuron and the ideal neuron, for its dimensions, next to it
    if doPlot, clf; 
        subplot('421'); imagesc(neur2d); title(sprintf('Neuron: %d', i));
        subplot('422'); 
    end
    [mserrors(i), similarity(i)] = ideal_comparison(neur2d, doPlot);
    
    % |--- plot the binary mask of the neuron
    if doPlot, subplot('423'); end
    numHoles(i) = count_holes(neur2d, doPlot);
    
    % |--- plot the local minimum mask of the neuron
    if doPlot, subplot('424'); end
    numLocmins(i) = count_local_minimums(neur2d, mn, doPlot);
    
    % |--- plot the distribution of the spatial footprint
    if doPlot, subplot(4,3,7:8); end
    notgaus(i) = gaussian_test(neuron.A(:,i), pval, doPlot);
    
    % |--- plot the calcium concentration of the neuron, useful as context
    if doPlot, subplot('414'); 
    plot_calcium(neuron.C_raw(i,:), neuron.C(i,:), max(neuron.A(:,i)), t);
    xlabel(str_xlabel);
    end
    
    % compare all test results with the user defined limits
    % if any test fails, then the neuron fails the text in general (get_pass)
    if numHoles(i) > ht, holes_pass = false; gen_pass = false; end
    if numLocmins(i) > mt, locmin_pass = false; gen_pass = false; end
    if mserrors(i) > elim, mse_pass = false; gen_pass = false; end
    if similarity(i) < simthres, sim_pass = false; gen_pass = false; end
    if rejng
        if notgaus(i) == 1, gaus_pass = false; gen_pass = false; end
    else
        gaus_pass = NaN;
    end
    
    % display a summary of all test results next to the ideal neuron shape
    if doPlot
        subplot('439'); axis off;
        if holes_pass, result1='\color{green}Pass'; 
        else result1='\color{red}Fail'; end
        if locmin_pass, result2='\color{green}Pass'; 
        else result2='\color{red}Fail'; end
        if mse_pass, result3='\color{green}Pass'; 
        else result3='\color{red}Fail'; end
        if sim_pass, result4='\color{green}Pass'; 
        else result4='\color{red}Fail'; end
        if gaus_pass, result5='\color{green}Pass'; 
        else result5='\color{red}Fail'; end
        text(0,0.5,{'Summary';'Holes: ';'Mins: ';'Shape: ';'Simil: ';'Gaussian: '}); 
        text(0.85,0.5,{'',result1,result2,result3,result4,result5});
    end
    
    % mark the neuron as bad if it doesn't pass one of the tests
    if ~ gen_pass, badCell(i) = true; end
    
    if dispFig
        % pause so the user can have a chance to look at the figure
        pause(interval);
    end
    
    if saveFigures
        % construct the file's name and save as a .png image
        if gen_pass, qual = 'good'; else, qual = 'bad'; end
        figName = fullfile(fig_out_dir,sprintf('neuron_%04d_%s.png', i, qual));
        saveas(h,figName);
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

if saveStatistics 
    % export statistics
    csvwrite(fullfile(out_dir, 'rbc_meansquarederror.txt'), mserrors);
    csvwrite(fullfile(out_dir, 'rbc_similarity.txt'), similarity);
    csvwrite(fullfile(out_dir, 'rbc_isnotgaussian.txt'), notgaus);
    csvwrite(fullfile(out_dir, 'rbc_numberholes.txt'), numHoles);
    csvwrite(fullfile(out_dir, 'rbc_numbermins.txt'), numLocmins);
    csvwrite(fullfile(out_dir, 'rbc_isbadcell.txt'), badCell);
end

if saveAsText
    % save footprints, traces and spikes as text files for analysis
    csvwrite(fullfile(out_dir, 'footprints.txt'), neuron.A);
    csvwrite(fullfile(out_dir, 'traces.txt'), neuron.C);
    csvwrite(fullfile(out_dir, 'spikes.txt'), neuron.S);
end

end
