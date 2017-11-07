function [ t, str_xlabel ] = calculate_timeline( nframes, fs )
%CALCULATE_TIMELINE Determines the time points in seconds when possible.
% Commonly used as a plot helper.
% Input:
%   nframes     :   number of frames
%   fs          :   frequency
% Output:
%   t           :   all time points
%   str_xlabel  :   'Time (Sec.)' if frequency is defined, otherwise 'Frame'

    t = 1:nframes;
    if ~isnan(fs)
        t = t/fs;
        str_xlabel = 'Time (Sec.)';
    else
        str_xlabel = 'Frame';
    end

end
