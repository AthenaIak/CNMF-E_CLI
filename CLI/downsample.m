function [Y, neuron, neuron_raw] = downsample(neuron_raw, data,mat_filename)
%DOWNSAMPLE downsamples raw data
%	neuron_raw      The Sources2D object containing information about the
%	analysis of the data.
%	data            The 3D data (x*y*num_frames)
%   mat_filename    The filename of the mat file to be loaded (only used if
%   options.ssub ~=1 and options.ssub ~= 1). 

%% downsample data for fast and better initialization
sframe=1;						% user input: first frame to read (optional, default:1)
numFrame = data.Ysiz(3,1);
num2read= numFrame;             % user input: how many frames to read   (optional, default: until the end)

tic;
if and(neuron_raw.options.ssub==1, neuron_raw.options.tsub==1)
    neuron = neuron_raw;
    Y = double(data.Y(:, :, sframe+(1:num2read)-1));
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else
    [Y, neuron] = neuron_raw.load_data(mat_filename, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
end
Y = neuron.reshape(Y, 1);
neuron_raw.P.p = 2;      %order of AR model

fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);