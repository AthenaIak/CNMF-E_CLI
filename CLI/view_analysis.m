%% view all steps of the analysis
% first select the folder that contains all relevant data (define path)
path = '~/tr/athina/Data/32366/cnmf_output006/';
cd ~/tr/athina/CNMF-E_CLI;
run_setup;
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
    disp('The folder has been created and old results will be overwritten. \n');
end
cd(folder_nm);
%%
c1= colormap(gray);
c2 = colormap(pink);
c_dual = [c1(1:2:64,:); c2(1:2:64,:)];
colormap(c_dual);
%%
figure(1); 

for i=1:neurons_detected-1
    for j=i+1:neurons_detected
        pw_corr = corr(neuron.A(:,i),neuron.A(:,j));
        if pw_corr > corr_thres
            subplot('421');
            spatialNeuronA = reshape(neuron.A(:,i),[d1 d2]);
            normFactorA = 31/max(max(spatialNeuronA));
            spatialNeuronA = floor(normFactorA*spatialNeuronA);
            spatialNeuronA(1,1)=64;
            imagesc(spatialNeuronA);
            title(sprintf('Neuron %d', i));
            subplot('422');
            spatialNeuronB = reshape(neuron.A(:,j),[d1 d2]);
            normFactorB = 32/max(max(spatialNeuronB));
            spatialNeuronB = floor(32+normFactorB*spatialNeuronB);
            spatialNeuronB(1,1)=0;
            imagesc(spatialNeuronB);
            title(sprintf('Neuron %d', j));
            subplot('412');
            %spatialNeuronA = reshape(neuron.A(:,i),[d1 d2]);
            %normFactorA = 31/max(max(spatialNeuronA));
            %spatialNeuronA = floor(normFactorA*spatialNeuronA);
            %spatialNeuronB = reshape(neuron.A(:,j),[d1 d2]);
            %normFactorB = 31/max(max(spatialNeuronB));
            %spatialNeuronB = floor(32+normFactorB*spatialNeuronB);
            addedSpatial = spatialNeuronA+spatialNeuronB;
            addedSpatial(1,1) = 0;
            addedSpatial(1,2) = 128;
            imagesc(addedSpatial);
            title('Neurons superimposed');
            %subplot('424');
            %spatialNeuronA = reshape(neuron.A(:,i),[d1 d2]);
            %normFactorA = 31/max(max(spatialNeuronA));
            %spatialNeuronA = floor(32+normFactorA*spatialNeuronA);
            %spatialNeuronB = reshape(neuron.A(:,j),[d1 d2]);
            %normFactorB = 31/max(max(spatialNeuronB));
            %spatialNeuronB = floor(normFactorB*spatialNeuronB);
            %imagesc(spatialNeuronB+spatialNeuronA);
            %title(sprintf('Neuron %d+%d', j,i));
            ax1 = subplot('413');
            plot(t, neuron.C_raw(i, :)*max(neuron.A(:, i)), 'b', 'linewidth', 2); hold on;
            plot(t, neuron.C(i, :)*max(neuron.A(:, i)), 'r'); hold off;
            ax2 = subplot('414');
            plot(t, neuron.C_raw(j, :)*max(neuron.A(:, j)), 'b', 'linewidth', 2); hold on;
            plot(t, neuron.C(j, :)*max(neuron.A(:, j)), 'r'); hold off;
            ax=axes('Units','Normal','Position',[.075 .095 .85 .85],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            allYLim = get([ax1,ax2], 'YLim');
            allYLim = cat(2, allYLim{:});
            set(ax1, 'YLim', [min(allYLim), max(allYLim)]);
            set(ax2, 'YLim', [min(allYLim), max(allYLim)]);
            title(sprintf('Neuron %d vs. %d: Spatial corr=%0.2f',i,j,pw_corr));
            saveas(gcf, sprintf('neurons_%03d_%03d_corr_%0.1f.png', i,j,pw_corr));
        end
    end
end

cd(cur_cd);
clear d1 d2 corr_thres neuron neurons_detected pw_corr;
clear cur_cd dir_neurons folder_nm;
clear ax ax1 ax2 str_xlabel allYLim;
clear i j t T;

