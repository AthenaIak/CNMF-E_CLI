function [ num_locmins ] = count_local_minimums( neur, mn, doPlot )
%COUNT_LOCAL_MINS Calculates the number of local minimums in a spatial 
%footprint.
%
% Input:
%   neur        :   2-Dimensional spatial footprint
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

% plot the binary mask of the local minimums
if doPlot
    imagesc(mod(locmin+1,2)); hold on;
    [x,y]=ind2sub(size(neur)-2,idx);
    plot(y+1,x+1,'rs','MarkerSize',12); hold off;
    title(sprintf('#Local mins: %d', num_locmins));
end

end
