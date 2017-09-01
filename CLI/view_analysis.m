%% view all steps of the analysis
% first select the folder that contains all relevant data (define path)
path = '~/tu/athina/Data/analyzed/32363/10.07.2017/ready_20170710_132114-008';
cd ~/tu/athina/CNMF-E_CI;
run_setup;
%% plot correlation image and peak-to-noise-ratio of the raw data
% can give you some ideas of how your data look like
load(fullfile(path,'f01-cn&pnr.mat'));

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
load(fullfile(path,'f02-cn_after_init.mat'));
figure;
imagesc(Cn);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;
clear center; clear Cn;

%% display contours of the neurons
load(fullfile(path,'f03-contours_after_init.mat'));
figure;
neuron.viewContours(Cn, 0.95, 0, [], 2);
colormap winter;
axis equal; axis off;
title('contours of estimated neurons');
clear neuron; clear Cn;

%% play the video data after subtracting the background components
load(fullfile(path,'f04-Ysignal_background_subtracted.mat'));
neuron.playMovie(Ysignal); 
clear neuron; clear Ysignal;

%% play the video data at the end
load(fullfile(path,'f05-Ysignal_end.mat'));
neuron.playMovie(Ysignal); 
clear neuron; clear Ysignal;

%% display neurons
load(fullfile(path,'f06-neurons.mat'));
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
load(fullfile(path,'f07-contours.mat'));
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
load(fullfile(path,'f06-neurons.mat'));
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

c1= colormap(gray);
c2 = colormap(pink);
c_dual = [c1(1:2:64,:); c2(1:2:64,:)];
colormap(c_dual);

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
            addedSpatial = spatialNeuronA+spatialNeuronB;
            addedSpatial(1,1) = 0;
            addedSpatial(1,2) = 128;
            imagesc(addedSpatial);
            title('Neurons superimposed');
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


%% measure hollowness of each neuron
load(fullfile(path,'f06-neurons.mat'));
figure;
for i=1:size(neuron.C,1)
neur = reshape(neuron.A(:,i), size(neuron.Cn));

neurClean = zeros(size(neur));
neurClean(find(neur>mean(neur(neur>0)))) = neur(find(neur>mean(neur(neur>0))));
%neurClean(find(neur>0)) = neur(find(neur>0));
[x,y] = ind2sub(size(neurClean),find(neurClean>0));
minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);
neurClean = neurClean(minx:maxx,miny:maxy);

[realCenterVal,centerX] = max(max(neurClean));
[~,centerY] = max(neurClean);
centerY = centerY(centerX);

ideal = zeros(size(neurClean));
halfSiz = min(min(centerY,size(neurClean,1)-centerY+1),min(centerX,size(neurClean,2)-centerX+1));
gSiz = halfSiz*2-1;
rowStart = centerY - halfSiz + 1;
colStart = centerX - halfSiz + 1;

%pd = fitdist(reshape(neurClean,[],1),'Normal');
%gSig = sqrt(std(pd));

ideal(rowStart:rowStart+gSiz-1,colStart:colStart+gSiz-1) = ...
 fspecial('gaussian', gSiz, gSiz/3);
gaussianCenterVal = ideal(rowStart+round(gSiz/2),colStart+round(gSiz/2));
gaussianFactor = realCenterVal/gaussianCenterVal;
ideal = ideal * gaussianFactor;

neur1d = reshape(neurClean,1,[]);
ideal1d = reshape(ideal,1,[]);
error = sqrt(sum(bsxfun(@minus,neur1d,ideal1d).^2,2));
%error = sqrt(sum(bsxfun(@minus,neur1d(neur1d==0),ideal1d(neur1d==0)).^2,2));
%error = error + ...
%     sqrt(sum(bsxfun(@minus,neur1d(ideal1d==0),ideal1d(ideal1d==0)).^2,2));
disp(sprintf('Neuron %d has error %.1f', i, error));

subplot('211');imagesc(neur);title(sprintf('Neuron %d', i));
subplot('223');imagesc(neurClean);title('Cleaned');
subplot('224');imagesc(ideal);title(sprintf('Error %.1f',error));
pause;
end

%% method 2: don't fit gaussian to center
load(fullfile(path,'f06-neurons.mat'));

% create dir to save images
folder_nm = fullfile(path,'neuron_hollowness');
cur_cd = cd();

if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
else
    disp('The folder has been created and old results will be overwritten. \n');
end
cd(folder_nm);

figure;
for i=1:size(neuron.C,1)
neur = reshape(neuron.A(:,i), size(neuron.Cn));

neurClean = zeros(size(neur));
%neurClean(find(neur>mean(neur(neur>0)))) = neur(find(neur>mean(neur(neur>0))));
neurClean(find(neur>max(max(neur))/5)) = neur(find(neur>max(max(neur))/5));
%neurClean(find(neur>0)) = neur(find(neur>0));
[x,y] = ind2sub(size(neurClean),find(neurClean>0));
minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);
neurClean = neurClean(minx:maxx,miny:maxy);

ideal = zeros(size(neurClean));
gSiz = min(size(neurClean));
rowStart = floor((size(ideal,1) - gSiz)/2) + 1;
colStart = floor((size(ideal,2) - gSiz)/2) + 1;

ideal(rowStart:rowStart+gSiz-1,colStart:colStart+gSiz-1) = ...
 fspecial('gaussian', gSiz, gSiz/4);
[realCenterVal,~] = max(max(neurClean));
realNormalFactor = 1/realCenterVal;
neurClean = neurClean * realNormalFactor;
gaussianCenterVal = ideal(rowStart+round(gSiz/2),colStart+round(gSiz/2));
gaussNormalFactor = 1/gaussianCenterVal;
ideal = ideal * gaussNormalFactor;

neur1d = reshape(neurClean,1,[]);
ideal1d = reshape(ideal,1,[]);
error = sqrt(sum(bsxfun(@minus,neur1d,ideal1d).^2,2));
if error>19
disp(sprintf('Neuron %d has error %.1f', i, error));
end

subplot('211');imagesc(neur);title(sprintf('Neuron %d', i));
subplot('223');imagesc(neurClean);title('Cleaned');
subplot('224');imagesc(ideal);title(sprintf('Error %.1f px:%d',error,numel(ideal)));

saveas(gcf, sprintf('neuron_%03d_holl_%0.1f.png', i,error));
end

cd(cur_cd);
clear x y maxx minx maxy miny;
clear neur neurClean ideal neur1d ideal1d;
clear colStart rowStart;
clear gSiz gaussianCenterVal realCenterVal realNormalFactor gaussNormalFactor;
clear i error folder_nm cur_cd;
clear dir_neurons neuron;
