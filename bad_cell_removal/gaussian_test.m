function [ rejected, P, stats ] = gaussian_test( neur, display )
%GAUSSIAN_TEST Checks if the normality hypothesis is rejected for the given
%spatial footprint
% Input:
%   neur    :   1-Dimensional spatial footprint
%   display :   creates a plot with the spatial footprint's histogram
% Output:
%   rejected    :   1 if the gaussian hypothesis was rejected (0 otherwise)
%   P           :   the P value
%   stats       :   the chi2gof statistics for the non-zero havues of the
%   neural spatial footprint

[rejected,P,stats] = chi2gof(neur(find(neur>0)));

% plot the histogram of the spatial footprint
if display
    % create a histogram from the chi2gof statistics
    bar(stats.edges(1:end-1)+diff(stats.edges),stats.O);
    
    % define title and axes text
    if rejected == 0
        conclusion = 'Possible';
    else
        conclusion = 'Rejected';
    end
    title(sprintf('%s, P: %.2f',conclusion,P));
    xlabel('Grey value (bins)');
    ylabel('Number of pixels');
end
end

