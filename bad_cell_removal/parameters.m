%% cell rejection parameters
% hole tolerance (maximum number of allowed holes in the spatial footprint)
ht = 0;

% gaussian error limit
elim = 6;

% remove neurons for whose spatial footprints the normality hypothesis is
% rejected
rng = true;     % reject non-gaussian

%% other options
display=false;      % displays a figure with all information if true
saveFigure = true;  % saves a figure for each neuron

