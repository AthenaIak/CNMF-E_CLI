%% view all steps of the analysis
% first select the folder that contains all relevant data (define path)
path = 'D:\Διπλωματική\CNMF_E\demos\data_endoscope\';

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

%% print centroids of seed pixels after initialization of A, C.
load(sprintf('%sf02-cn_after_init.mat',path));
figure;
imagesc(Cn);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;

%% display contours of the neurons
load(sprintf('%sf03-contours_after_init.mat',path));
figure;
neuron.viewContours(Cn, 0.95, 0, [], 2);
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');

