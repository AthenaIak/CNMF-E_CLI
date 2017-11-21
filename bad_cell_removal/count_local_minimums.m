function [ num_locmins, error, norm_error ] = count_local_minimums( neur, mn, doPlot )
%COUNT_LOCAL_MINS Calculates the number of local minimums in a spatial 
%footprint.
%
% Input:
%   neur        :   2-Dimensional spatial footprint (assumes that the
% boundaries are cropped)
%   mn          :   number of Neighbours examined for the locam Minimum
%   calculation (can be either 4 or 8).
%   doPlot      :   creates a plot with the binary mask of the spatial 
% footprint
%
% Output:
%   num_locmins :   the total number of local minimums found inside the
%   footprint

% find the local minimums
locmin = imregionalmin(neur, mn);

% count the number of local minimums in the spatial footprint
% 1. do not take into account the outer rows and columns)
% 2. it's a local minimum only if the pixel is non-negative in the spatial
% footprint
idx = find(and(locmin(2:end-1,2:end-1)==1,neur(2:end-1,2:end-1)~=0));
num_locmins = length(idx);

[x,y]=ind2sub(size(neur)-2,idx);

% calculate accumulated error (values of neighbouring pixels)
error = 0;
for idx = 1:length(x)
    xi = x(idx)+1;
    yi = y(idx)+1;
    error = error + sum(sum(neur(xi-1:xi+1,yi-1:yi+1)))-neur(xi,yi);
end

% calculate error when spatial footprint values are normalized to be in the range 0-1
norm_error = error/max(max(neur));

% plot the binary mask of the local minimums
if doPlot
    imagesc(mod(locmin+1,2)); hold on;
    [x,y]=ind2sub(size(neur)-2,idx);
    plot(y+1,x+1,'rs','MarkerSize',12); hold off;
    title(sprintf('#Local mins: %d', num_locmins));
end

end
