function [ neuron, good_cells ] = remove_bad_cells( neuron, classifier_matname, pred_cutoff )
%REMOVE_BAD_CELLS Detects bad cells based on a trained logistic regression
%classifier.
%
% Input: 
%   neuron              :   a neuron object containing the spatial footprints
%   classifier_matname  :   the name of the classifier file (normally
% trained using script "train_bc_classifier"). It contains the indexes of
% the attributes used for classification (kept_vars) and the learned
% coefficients (coefficients) for the logistic regression.
%   pred_cutoff         :   spatial footprints receiving a prediction 
% higher than the specified prediction of being good cells, will be marked 
% as good cells and will not be removed.
%
% Output:
%   neuron              :   the neuron object with the bad neurons removed.
%   good_cells          :   the indexes of the good cells that were not
% deleted.


%% initialization
% use the default values if not specified
if ~ exist('classifier_matname','var')
    classifier_matname = 'classifier.mat';
end
if ~ exist('pred_cutoff','var')
    pred_cutoff = 0.5;
end

% load the logistic regression model
model = load(classifier_matname);

%% extract statistics for all neurons

nNeurons = size(neuron.C,1);    % number of neurons detected by CNMF-E

% initialize a matrix to save statistics per frame
results = zeros(nNeurons,18);   % mean squared error when compared to an ideal shape

for i=1:nNeurons % for each neuron
    % get the spatial footprints in 2-dimensional form
    neur2d = reshape(neuron.A(:,i), size(neuron.Cn));
    [x,y] = ind2sub(size(neur2d),find(neur2d>0)); % find the indexes to ...
    neur2d = neur2d(min(x):max(x),min(y):max(y)); % ... crop the boundaries
    
    % calculate statistics for this spatial footprint 
    if any(ismember(model.kept_vars, 1:2))
        [results(i,1), results(i,2)] = ideal_comparison(neur2d, false);
    end
    if any(ismember(model.kept_vars, [3:5 12:17])) 
        [results(i,3), results(i,4), results(i,5)] = count_holes(neur2d, false);
    end
    if any(ismember(model.kept_vars, [6:8 12:14]))
        [results(i,6), results(i,7), results(i,8)] = count_local_minimums(neur2d, 4, false);
    end
    if any(ismember(model.kept_vars, [9:11 15:17]))
    [results(i,9), results(i,10), results(i,11)] = count_local_minimums(neur2d, 8, false);
    end 
    results(i,12:14) = results(i,3:5)+results(i,6:8);
    results(i,15:17) = results(i,3:5)+results(i,9:11);
    if any(ismember(model.kept_vars, 18))
        [~,results(i,18)] = gaussian_test(neuron.A(:,i), 0.05, false);
    end
end


%% classify the neurons
pred_prob = mnrval(model.coefficients,results(:,model.kept_vars));

% convert prediction probabilities to binary predictions
good_cells = single(pred_prob(:,1)>pred_cutoff); 

%% remove bad cells
neuron.delete(~good_cells);

end

