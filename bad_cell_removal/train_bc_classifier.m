%% options
results_file = '~/Data/badcells/20170724_102838-014/hresults.csv';
gtruth_file = '~/Data/badcells/20170724_102838-014/neurqual.csv';
do_test = true;         % test the model (or use all data to train)
test_percent = 0.25;    % ratio of data used for testing (rest are used for training)
pred_cutoff = 0.5;      % cells with good cell probability higher than that will be classified as good cells
classifier_matname = 'classifier.mat';

%% Initialization
% load data
results = csvread(results_file);
gtruth = csvread(gtruth_file);

nNeurons = size(results,1);

%% Note: Variables in results
% (0. ID)
% 1. mean squared error     2. similarity
% 3. number of holes        4. holes error      5. holes normalized error
% 6. number of local mins   7. locmin error     8. locmin normalized error (for 4 connections)
% 9. number of local mins   10. locmin error    11. locmin normalized error (for 8 connections)
% 12. (3) plus (6)          13. (4) plus (7)    14. (5) plus (8)
% 15. (3) plus (9)          16. (4) plus (10)   17. (5) plus (11)
% 18. p value of gaussian test

%% prepare data for logistic regression
% shuffle the data
sort_order = randperm(nNeurons);
results = results(sort_order,:);
gtruth = gtruth(sort_order,:);
gtruth(find(gtruth==0))=2; % 1=good cell, 2=bad cell

if sum(results(:,1)~=gtruth(:,1)) ~= 0, disp('Something went wrong!'); end
gtruth = categorical(gtruth);

if do_test
    % split to training and test sets and remove the neuron ID column
    split_idx = floor((1 - test_percent) * nNeurons);
    training_vars = results(1:split_idx,2:end);
    training_class = gtruth(1:split_idx,2);
    test_vars = results(split_idx+1:end,2:end);
    test_class = gtruth(split_idx+1:end,2);
else
    % only remove the neuron ID column
    training_vars = results(:,2:end);
    training_class = gtruth(:,2);
end

%% train model
kept_vars = [1 2 3 9 17];  % mean accuracy (5-fold validation): 0.8963

B = mnrfit(training_vars(:,kept_vars),training_class,'model','nominal');

%% test model
if do_test
    pred_prob = mnrval(B,test_vars(:,kept_vars));

    % convert to same format with test_class
    prediction = single(pred_prob(:,1)>pred_cutoff); % 1=good cell, 0=bad cell
    prediction(find(prediction==0))=2;
    prediction = categorical(prediction);

    % evaluate
    trues = (prediction==test_class);
    positives = (test_class==categorical((1))); % positive = good cell
    tp = and( trues, positives);
    tn = and( trues,~positives);
    fn = and(~trues, positives); % found as negatives when trully positives
    fp = and(~trues,~positives); % found as positives when trully negatives
    totaltp = sum(tp);
    totalfp = sum(fp);
    totaltn = sum(tn);
    totalfn = sum(fn);
    fprintf('Hit rate (TPR): %f\n',totaltp/(totaltp+totalfn));
    fprintf('Specificity (TNR): %f\n',totaltn/(totaltn+totalfp));
    fprintf('Precision (PPV): %f\n',totaltp/(totaltp+totalfp));
    fprintf('Accuracy (ACC): %f\n',(totaltp+totaltn)/(totaltp+totaltn+totalfp+totalfn));
end

%% output model
coefficients = B;
save(classifier_matname, 'coefficients', 'kept_vars', '-v7.3');
