cnmfe_setup;
addpath('statistical_analysis');
%%
movieFiles={
    '20170724_094625', ...
    '20170724_100129', ...
    '20170724_101503', ...
    '20170724_102838', ...
    '20170724_104348', ...
    '20170725_092500', ...
    '20170725_093801', ...
    '20170725_095032', ...
    '20170725_100338', ...
    '20170725_101701', ...
    '20170726_103537', ...
    '20170726_104856', ...
    '20170726_110209', ...
    '20170726_111616', ...
    '20170726_113000', ...
    '20170727_094127', ...
    '20170727_095434', ...
    '20170727_100818', ...
    '20170727_102228', ...
    '20170727_103727', ...
    '20170728_095108', ... %session 2 starts:
    '20170717_092917', ...
    '20170717_094631', ...
    '20170717_100129', ...
    '20170717_101642', ...
    '20170717_103145', ...
    '20170718_093725', ...
    '20170718_095120', ...
    '20170718_100408', ...
    '20170718_101708', ...
    '20170718_103020', ...
    '20170719_100218', ...
    '20170719_101523', ...
    '20170719_102754', ...
    '20170719_104033', ...
    '20170719_105244', ...
    '20170720_095256', ...
    '20170720_100558', ...
    '20170720_101938', ...
    '20170720_103237', ...
    '20170720_104516', ...
    '20170721_093828'
};

path = '~/cb/Data/';
time_cut = 30; % (cut first 30 secs)
fs = 5; % Hz
maxIter=10000;

numFiles = length(movieFiles);
corrPerTrial = zeros(1,numFiles);
%corrPerTrial2 = zeros(1,numFiles);

%% calculate thresholds
for m=1:numFiles
    filename = fullfile([path movieFiles{m} '-014'],'data.mat');
    load(filename);
    
    % take the inferred spiking activity of each trial
    fs=5; time_cut=30;
    spikes = neuron.S(:,1+time_cut*fs:end);
    
    nn = size(spikes,1); % number of neurons and timesteps
    
    % detect problematic frames and remove them
    % yuste
    isSpike=spikes>0;  
    spikes = neuron.C(:,1+time_cut*fs:end);
    spikes=spikes(:,find(sum(isSpike)<=0.1*nn)); % remove the bad frames    
    
    % second
    %spikes=single(spikes>0);   
    %spikes=spikes(:,find(sum(spikes)<=0.1*nn)); % remove the bad frames
    
    % 1st
    %isSpike=spikes>0;   
    %spikes=spikes(:,find(sum(isSpike)<=0.1*nn)); % remove the bad frames
    spikes = spikes(:,1:1140+time_cut); % make sure all recordings are 3.9min long
    random_spikes = spikes; % backup of the original
    
    [nn, ts] = size(spikes); % number of neurons and timesteps
    
    % shift rows randomly
    fprintf('Spike shifting... '); tic;
    for i=1:maxIter
        shifts=randi(ts,nn,1); % random values for the shifting for each row
        r = rem(shifts,ts);
        c = [random_spikes,random_spikes];
        random_spikes = c(bsxfun(@plus,bsxfun(@plus,ts - r,0:ts-1)*nn,(1:nn)'));
        if mod(i,maxIter/10) == 0 
            fprintf('%d%% ',i/maxIter*100);
        end
    end
    fprintf('\nDone\n');toc;
    
    random_sindex = pairwise_sindex(random_spikes);
    sindex_threshold = prctile(random_sindex,99);
    
    real_sindex = pairwise_sindex(spikes);
    corrPerTrial(m)=sum(real_sindex>sindex_threshold)/length(random_sindex);    
    
    fprintf('Movie %s significance ratio: %f\n', movieFiles{m},corrPerTrial(m));
end

disp(corrPerTrial');
save('results_clean_statistical_ca2.mat','corrPerTrial','-v7.3');
%%
corrPerDay1 = [sum(corrPerTrial(1:5))/5 sum(corrPerTrial(6:10))/5 ... 
    sum(corrPerTrial(11:15))/5 sum(corrPerTrial(16:20))/5 corrPerTrial(21)];

corrPerDay2 = [sum(corrPerTrial(22:26))/5 sum(corrPerTrial(27:31))/5 ... 
    sum(corrPerTrial(32:36))/5 sum(corrPerTrial(37:41))/5 corrPerTrial(42)];

plot(1:5,corrPerDay1);hold on;plot(1:5,corrPerDay2);hold off;
legend('Overlapping','Random');
xlabel('Day');ylabel('Significant S-index ratio');
%axis([1 5 0 0.035]);

anova2([corrPerDay1;corrPerDay2]',1)
% columns: 2 conditions (first row=session 3, second=session 2)
% rows: 5 days
% we expect to find Prob>F < 0.01 (comparison between conditions) ???


