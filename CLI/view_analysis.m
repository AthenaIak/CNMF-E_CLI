%% view all steps of the analysis
% first select the folder that contains all relevant data (define path)
path = 'D:\�����������\CNMF_E\demos\data_endoscope\';

%% plot correlation image and peak-to-noise-ratio of the raw data
% can give you some ideas of how your data look like
load(sprintf('%sf01-cn&pnr.mat',path));

figure('position', [10, 500, 1776, 400]);
subplot(131);
imagesc(Cn, [0, 1]); colorbar;
axis equal off tight;
title('correlation image');

subplot(132);
imagesc(pnr,[0,max(pnr(:))*0.98]); colorbar;
axis equal off tight;
title('peak-to-noise ratio');

subplot(133);
imagesc(Cn.*pnr, [0,max(pnr(:))*0.98]); colorbar;
axis equal off tight;
title('Cn*PNR');
clear Cn; clear pnr;

%% print centroids of seed pixels after initialization of A, C.
load(sprintf('%sf02-cn_after_init.mat',path));
figure;
imagesc(Cn);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;
clear center; clear Cn;

%% display contours of the neurons
load(sprintf('%sf03-contours_after_init.mat',path));
figure;
neuron.viewContours(Cn, 0.95, 0, [], 2);
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');
clear neuron; clear Cn;

%% play the video data after subtracting the background components
load(sprintf('%sf04-Ysignal_background_subtracted.mat',path));
neuron.playMovie(Ysignal); 
clear neuron; clear Ysignal;

%% play the video data at the end
load(sprintf('%sf05-Ysignal_end.mat',path));
neuron.playMovie(Ysignal); 
clear neuron; clear Ysignal;

%% display neurons
load(sprintf('%sf06-neurons.mat',path));
if exist('dir_neurons', 'dir')
    temp = cd();
    cd(dir_neurons);
    delete *;
    cd(temp);
else
    mkdir(dir_neurons);
end
neuron.viewNeurons([], neuron.C_raw, dir_neurons);
clear dir_neurons; clear neuron;

%% display contours of the neurons
load(sprintf('%sf07-contours.mat',path));
figure;
% neuron.viewContours(Cn, 0.9, 0);  % correlation image computed with
% spatially filtered data
% [Cn, pnr] = neuron.correlation_pnr(Ysignal); % very slow 
Cnn = correlation_image(Ysignal(:, 1:5:end), 4, d1, d2);
neuron.viewContours(Cnn, 0.9, 0); % correlation image computed with background-subtracted data
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');
% plot contours with IDs
figure;
plot_contours(neuron.A, Cn, 0.9, 0, [], neuron.Coor);
clear neuron; clear Ysignal; clear d1; clear d2; clear Cn; clear Cnn;


