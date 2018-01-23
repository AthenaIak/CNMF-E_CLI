%% cell rejection parameters
% maximum number of allowed holes in the spatial footprint
ht = 0;                 % hole tolerance
mt = 2;                 % local minimum tolerance
mn = 8;                 % number of neighbours checked for the local 
                        % minimum calculation (4 or 8)

elim = 0.035;           % gaussian error limit
simthres = 0.45;       % similarity threshold

rejng = true;           % reject non-gaussian
pval = 0.025;           % minimum significance level for rejecting the 
                        % normal distribution hypothesis 
                        % (the default p value is 0.05)

%% other options
dispFig=false;          % displays a figure with all information if true
interval=1;             % time interval in between figures (in seconds)

% save options
saveFigures = false;    % saves a figure for each neuron (as .png images)
saveUpdNeuron = true;   % save the updated neuron object (as updated_data.mat)
saveAsText = true;      % saves footprints, traces and spikes of kept neurons 
saveForCellReg = true;  % saves the footprints in a form compatible with CellReg
saveStatistics = true;  % saves all statistics for this neuron
