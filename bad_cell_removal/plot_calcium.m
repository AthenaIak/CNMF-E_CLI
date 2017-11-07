function [ output_args ] = plot_calcium( C_raw, C, maxA, t )
%PLOT_CALCIUM Creates a plot of the neuron's calcium concentration over
%time.
%   Detailed explanation goes here


    plot(t, C_raw*maxA, 'b', 'linewidth', 2); hold on;
    plot(t, C*maxA, 'r'); hold off;
    title('Calcium concentration over time');
end

