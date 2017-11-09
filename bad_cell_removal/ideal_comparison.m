function [ error ] = ideal_comparison( neur, doPlot )
%IDEAL_COMPARISON creates a gaussian filter with similar dimensions to the 
%given spatial footprint (ideal shape) and calculates the mean squared 
%error between the two matrices.
% Input:
%   neur    :   2-Dimensional spatial footprint
%   doPlot  :   if set true, creates a plot of the ideal spatial footprint
% Output: 
%   error   :   the minimum squared error of all pixels compared with an
%   ideal shape of the neuron    

    % new method:
    [x,y] = ind2sub(size(neur),find(neur>0));
    cropped_neur = neur(min(x):max(x),min(y):max(y));
    neur1d = reshape(cropped_neur,1,[]);
    
    % normalize the spatial footprint (0-1)
    %neur1d = neur1d/max(neur1d);
    ccc = normcdf(neur1d,mean(neur1d),std(neur1d));
    ccc(find(ccc==min(ccc)))=0;
    neur1d = neur1d/max(neur1d);
    error2 = sqrt(sum(bsxfun(@minus,neur1d,ccc).^2,2));
    
    % old method:
    % clear some noise
    neurClean = zeros(size(neur));
    neurClean(find(neur>max(max(neur))/5)) = neur(find(neur>max(max(neur))/5));

    
    %noiseLevel = mean(neur1d(find(neur1d>0)))-3*std(neur1d(find(neur1d>0)))
%neurClean(find(neur2d>noiseLevel)) = neur2d(find(neur2d>noiseLevel));
    [x,y] = ind2sub(size(neurClean),find(neurClean>0));
    neurClean = neurClean(min(x):max(x),min(y):max(y)); % crop blank boundaries

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

    if doPlot
        imagesc(ideal);
        title(sprintf('Ideal err: %.2f, new: %.2f', error,error2));
    end
end

