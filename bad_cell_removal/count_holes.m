function [ num_holes ] = count_holes( neur, doPlot )
%COUNT_HOLES Calculates the number of holes in a spatial footprint.
%
% Input:
%   neur        :   2-Dimensional spatial footprint
%   doPlot      :   creates a plot with the binary mask of the spatial 
% footprint
%
% Output:
%   num_holes   :   the total number of pixels that were empty in the 
% binary mask created based on the input spatial footprint

% binarize the spatial footprint
L = imbinarize(neur,graythresh(neur));

% calculate the number of holes in the binary mask
filled = imfill(L,'holes'); % fill holes in the shape
num_holes = sum(sum(filled-L)); % subtract the non-filled shape

% plot the binary mask of the spatial footprint
if doPlot
    % get rid of the blank boundaries
    [x,y] = ind2sub(size(L),find(L>0));
    imagesc(L(min(x):max(x),min(y):max(y)));
    title(sprintf('#Holes: %d', num_holes));
end
end

