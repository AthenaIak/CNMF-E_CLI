clear;clc;close;
disp('Done 1');

%%

name = '20170726_111616.mat';
%name = '20170726_113000.mat';
inDir = '~/cb/Data/';

filename = fullfile(inDir,name);
load(filename);

figure; colormap gray;

disp('Done 2');

%Y = reshape(double(Y), 1080*1440, []);
%imagesc(reshape(Y(:,1), 1080, 1440, []));
%errs = sum(bsxfun(@minus,Y(:,1),Y).^2)./size(Y,1);

%plot(errs);
%%
for i=1:size(Y,3)
    imagesc(Y(:,:,i));title(sprintf('Frame %04d',i));
    pause(0.02);
end