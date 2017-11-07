function [] = plot_cropped_neuron( neur, idx )
%PLOT_CROPPED_NEURON creates a plot of the spatial footprint after cropping
%the image boundaries.
% Input:
%   neur    :   2-Dimensional spatial footprint
%   idx     :   the index of the neuron

% find the coordinates to crop the boundaries of the mask
    [x,y] = ind2sub(size(neur),find(neur>0));
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    
    imagesc(neur(minx:maxx,miny:maxy));
    title(sprintf('Neuron: %d', idx));
end
