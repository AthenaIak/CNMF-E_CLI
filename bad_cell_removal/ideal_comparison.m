function [ error ] = ideal_comparison( neur, display )
%IDEAL_COMPARISON creates a gaussian filter with similar dimensions to the 
%given spatial footprint (ideal shape) and calculates the mean squared 
%error between the two matrices.
% Input:
%   neur    :   2-Dimensional spatial footprint
%   display :   if set true, creates a plot of the ideal spatial footprint

neurClean = zeros(size(neur));
neurClean(find(neur>max(max(neur))/5)) = neur(find(neur>max(max(neur))/5));

%neurClean(find(neur2d>noiseLevel)) = neur2d(find(neur2d>noiseLevel));
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

    if display
        imagesc(ideal);
        title(sprintf('Ideal err: %.2f', error));
    end
end

