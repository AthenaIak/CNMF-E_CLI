%% Initialization
% load data
data='~/Data/20170724_102838-014/data.mat';
load(data);                 % contains variable named neuron with all data

% get parameters from script
%run('parameters');

% determine the output directory if not defined
out_dir = fileparts(data);
nNeurons = size(neuron.C,1);    % number of neurons detected by CNMF-E

% initialize a matrix to save statistics per frame
results = zeros(nNeurons,18);   % mean squared error when compared to an ideal shape
% mean squared error, similarity, 
% number of holes, accummulative error, normalized accummulative error
% number of local mins, accummulative error, normalized accummulative error for 4 connections.
% number of local mins, accummulative error, normalized accummulative error for 8 connections.
% sum of holes+local mins, error, normalized error for 4 connections
% sum of holes+local mins, error, normalized error for 8 connections
% pvalue of gaussian test

%% perform tests to all detected neurons
for i=1:nNeurons % for each neuron
    
    % get the spatial footprints in 2-dimensional form
    neur2d = reshape(neuron.A(:,i), size(neuron.Cn));
    [x,y] = ind2sub(size(neur2d),find(neur2d>0)); % find the indexes to ...
    neur2d = neur2d(min(x):max(x),min(y):max(y)); % ... crop the boundaries
    
    % calculate statistics for this spatial footprint 
    [results(i,1), results(i,2)] = ideal_comparison(neur2d, false);
    [results(i,3), results(i,4), results(i,5)] = count_holes(neur2d, false);
    [results(i,6), results(i,7), results(i,8)] = count_local_minimums(neur2d, 4, false);
    [results(i,9), results(i,10), results(i,11)] = count_local_minimums(neur2d, 8, false);
    results(i,12:14) = results(i,3:5)+results(i,6:8);
    results(i,15:17) = results(i,3:5)+results(i,9:11);
    [~,results(i,18)] = gaussian_test(neuron.A(:,i), 0.05, false);
end

%% save results
% header: 
% id,mse,similarity,nholes,herr,hnerr,nlm4,lm4err,lm4nerr,nlm8,lm8err,lm8nerr,nhlm4,hlm4err,hlm4nerr,nhlm8,hlm8err,hlm8nerr,pval
csvwrite(fullfile(out_dir, 'hresults.csv'), [(1:nNeurons)' results]);