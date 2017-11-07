%% cell rejection parameters
% hole tolerance (maximum number of allowed holes in the spatial footprint)
ht = 0;

% gaussian error limit
elim = 6;

% remove neurons for whose spatial footprints the normality hypothesis is
% rejected
rejng = true;           % reject non-gaussian

%% other options
dispFig=false;          % displays a figure with all information if true

% save options
saveFigures = true;     % saves a figure for each neuron
figFormat = 'png';      % only used if saveFigures is true (see function saveas)
saveUpdNeuron = true;   % save the updated neuron object (as updated_data.mat)
saveAsText = true;      % saves footprints, traces and spikes of kept neurons 
saveForCellReg = true;  % saves the footprints in a form compatible with CellReg
