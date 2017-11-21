function [ num_holes, error, norm_error ] = count_holes( neur, doPlot )
%COUNT_HOLES Calculates the number of holes in a spatial footprint.
%
% Input:
%   neur        :   2-Dimensional spatial footprint (assumes that the
% boundaries are cropped)
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

% find the indexes of the holes
[x,y] = ind2sub(size(L), find(filled-L == 1));

% calculate accumulated error (values of neighbouring pixels)
error = 0;
for idx = 1:length(x)
    xi = x(idx);
    yi = y(idx);
    error = error + sum(sum(neur(xi-1:xi+1,yi-1:yi+1)))-neur(xi,yi);
end

% calculate error when spatial footprint values are normalized to be in the range 0-1
norm_error = error/max(max(neur));

% plot the binary mask of the spatial footprint
if doPlot
    imagesc(L);
    title(sprintf('#Holes: %d', num_holes));
end
end

