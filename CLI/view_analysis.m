%% view all steps of the analysis
% first select the folder that contains all relevant data (define path)
path = '~/tp/data/iHPC5 raw/recording_20160125_114832/output/mcorr_mosaic_128-recording_20160125_114832/';

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

%% view pair-wise highly spatially correlated neurons
load(sprintf('%sf06-neurons.mat',path));

corr_thres = .5;
neurons_detected = size(neuron.A,2);
d1 = neuron.options.d1;
d2 = neuron.options.d2;

T = size(neuron.C, 2);
t = 1:T;
if ~isnan(neuron.Fs)
    t = t/neuron.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

% create dir to save images
folder_nm = sprintf('%s%s',path,'neuron_comparison');
cur_cd = cd();

if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
else
    fprintf('The folder has been created and old results will be overwritten. \n');
end
cd(folder_nm);

for i=1:neurons_detected-1
    for j=i+1:neurons_detected
        pw_corr = corr(neuron.A(:,i),neuron.A(:,j));
        if pw_corr > corr_thres
            figure(1); 
            subplot('321');
            imagesc(reshape(neuron.A(:,i),[d1 d2]));
            subplot('322');
            imagesc(reshape(neuron.A(:,j),[d1 d2]));
            subplot('312');
            plot(t, neuron.C_raw(i, :)*max(neuron.A(:, i)), 'b', 'linewidth', 2); hold on;
            plot(t, neuron.C(i, :)*max(neuron.A(:, i)), 'r'); hold off;
            subplot('313');
            plot(t, neuron.C_raw(j, :)*max(neuron.A(:, j)), 'b', 'linewidth', 2); hold on;
            plot(t, neuron.C(j, :)*max(neuron.A(:, j)), 'r'); hold off;
            ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            title(sprintf('Neuron %d vs. %d: Spatial corr=%0.2f',i,j,pw_corr));
            saveas(gcf, sprintf('neurons_%03d_%03d_corr_%0.1f.png', i,j,pw_corr));
            %pause();
        end
    end
end


