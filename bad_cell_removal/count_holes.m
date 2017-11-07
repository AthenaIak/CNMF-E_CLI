function [ num_holes ] = count_holes( neur, display )
%COUNT_HOLES Calculates the number of holes in a spatial footprint
% Input:
%   neur    :   2-Dimensional spatial footprint
%   display :   creates a plot with the binary mask of the spatial footprint

% turn footprint to a binary matrix
L = bwlabel(neur);

% calculate the number of holes in the mask
filled = imfill(L,'holes'); % fill holes in the shape
num_holes = sum(sum(filled-L)); % subtract the non-filled shape

% plot the binary mask of the spatial footprint
if display
    % find the coordinates to crop the boundaries of the mask
    [x,y] = ind2sub(size(L),find(L>0));
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    imagesc(L(minx:maxx,miny:maxy));
    title(sprintf('#Holes: %d', num_holes));
end
end

