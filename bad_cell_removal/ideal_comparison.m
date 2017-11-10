function [ error,similarity ] = ideal_comparison( neur, doPlot )
%IDEAL_COMPARISON creates a gaussian filter with similar dimensions to the 
%given spatial footprint (ideal shape) and calculates the mean squared 
%error between the two matrices.
% Input:
%   neur    :   2-Dimensional spatial footprint (assumes that the
% boundaries are cropped)
%   doPlot  :   if set true, creates a plot of the ideal spatial footprint
% Output: 
%   error   :   the minimum squared error of all pixels compared with an
%   ideal shape of the neuron    

    % normalize the spatial footprint's values to be in the 0-1 range
    [realCenterVal,~] = max(max(neur));
    realNormalFactor = 1/realCenterVal;
    neur = neur * realNormalFactor;
    
    % create a gaussian filter with the same size and approximately same sigma
    mean_sigma = mean(size(neur))*std(reshape(neur,1,[]));
    ideal = fspecial('gaussian', size(neur), mean_sigma);
    % normalize its values to be in the 0-1 range
    normalFactor = max(max(ideal));
    ideal = ideal / normalFactor;
    
    % calculate error metrics
    error = immse(neur,ideal);
    similarity = ssim(neur,ideal);
    
    if doPlot
        % plot the ideal neuron
        imagesc(ideal);
        title(sprintf('err:%.2f sim:%.3f', error*10,similarity));
    end
end

