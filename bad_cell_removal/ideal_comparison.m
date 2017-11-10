function [ error ] = ideal_comparison( neur, doPlot )
%IDEAL_COMPARISON creates a gaussian filter with similar dimensions to the 
%given spatial footprint (ideal shape) and calculates the mean squared 
%error between the two matrices.
% Input:
%   neur    :   2-Dimensional spatial footprint (assumes that boundaries
%   are cropped)
%   doPlot  :   if set true, creates a plot of the ideal spatial footprint
% Output: 
%   error   :   the minimum squared error of all pixels compared with an
%   ideal shape of the neuron    

    % create a square gaussian filter (almost) the same dimensions 
    ideal = zeros(size(neur));
    gSiz = min(size(neur));
    rowStart = floor((size(ideal,1) - gSiz)/2) + 1;
    colStart = floor((size(ideal,2) - gSiz)/2) + 1;
    ideal(rowStart:rowStart+gSiz-1,colStart:colStart+gSiz-1) = ...
     fspecial('gaussian', gSiz, gSiz/3);

    % normalize the spatial footprint's values to be in the 0-1 range
    [realCenterVal,~] = max(max(neur));
    realNormalFactor = 1/realCenterVal;
    neur = neur * realNormalFactor;
    
    % I don't remember how I decided this. It works pretty well, but it 
    % doesn't seem to have any mathematical premise
    %
    %gaussianCenterVal = ideal(rowStart+round(gSiz/2),colStart+round(gSiz/2));
    %gaussNormalFactor = 1/gaussianCenterVal;
    %ideal = ideal * gaussNormalFactor;
    %
    % so, after thinking about this, I think this was the rationale:
    % We optimize the gaussian filter to be as similar as possible to the
    % hgihest value ('centroid') of the spatial footprint. In that way, we
    % minimize error that comes from the higher values and we depend more
    % on errors away from the center. I think at the end, what this metric
    % really measured was the distance of the centroid from the center of
    % the picture (by dividing the idea with the centroid value, I pushed 4
    % values above the zero. One of them was above the centroid and so had
    % error=0, the other 3 were symmetrical to the centroid. The farther 
    % away the centroid was from the center of the picture, the higher the
    % error in the other 3 points. Although it kinda worked, I don't like
    % how arbitrary this metric is. But something that utilizes the
    % information about the centroid may make sense.
    
    % without the previous lines, it doesn't work at all
    %neur1d = reshape(neur,1,[]);
    %ideal1d = reshape(ideal,1,[]);
    %error = sqrt(sum(bsxfun(@minus,neur1d,ideal1d).^2,2));
    % ^ oime, i don't take the mean for the error D:D:D:D:
    
    % try 3 more metrics
    peaksnr = psnr(neur,ideal);
    ssimval = ssim(neur,ideal) ;
    error = immse(neur,ideal);
    % none worked!
    
    if doPlot
        imagesc(ideal);
        title(sprintf('n:%.2f s:%.4f e:%.3f', peaksnr, ssimval, error));
    end
end

