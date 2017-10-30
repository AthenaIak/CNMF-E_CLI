function [ output_args ] = evaluate_parameters( movieFiles, parIDs )
%EVALUATE_PARAMETERS Summary of this function goes here
% input:
%   path    :   path to the parameter files 
%   parIDs  :   array containing the IDs of the parameter files that follow 
% the format parameters_cnmf_an[ID].m

run ('CellReg/setup')

% evaluation parameters
spatial_correlation_threshold = 0.3;

% CellReg parameters
results_dir = '/home/athina/tu/athina/CNMF-E_CLI/CellReg/SampleData/SampleResults';
file_names = {'/home/athina/tu/athina/CNMF-E_CLI/CellReg/SampleData/spatial_footprints_01.mat', ...
    '/home/athina/tu/athina/CNMF-E_CLI/CellReg/SampleData/spatial_footprints_02.mat', ...
    '/home/athina/tu/athina/CNMF-E_CLI/CellReg/SampleData/spatial_footprints_03.mat', ...
    '/home/athina/tu/athina/CNMF-E_CLI/CellReg/SampleData/spatial_footprints_04.mat', ...
    '/home/athina/tu/athina/CNMF-E_CLI/CellReg/SampleData/spatial_footprints_05.mat'};

for id = parIDs
    [results_dir, footprint_files] = run_analysis( movieFiles, tag );
    [ips, fps] = run_cellreg(results_dir, file_names); 
    % ips = cell index per session
    % fpr = spatial footprint of each neuron per session
    
    single_cells = ips(sum(ips>0,2)==1,:); % cells found only in one session
    
    for session=1:size(fps,2)
        d = size(fps{session});
        A = reshape(fps{session},d(1),d(2)*d(3))';
        
        % calculate the maximum correlation of each neuron with any other
        % neuron for the current session
        sp_corr = corr(A,A);
        sp_corr = sp_corr - eye(d(1)); % all neurons have correlation = 1 with themselves
        maxCorr = max(sp_corr); % max spatial correlation of each neuron with any other neuron
        
        % get the single cells in the current session
        curr_single = single_cells(single_cells(:,session)>0);
        
        allCells=1:d(1); % all cells for this session
        fs_count = 0; % false splitting count
        fd_count = 0; % false detection count (cell is actually noise)
        ts_count = 0;
        td_count = 0;
        for cell=allCells
            if sum(curr_single==cell) 
                isSingle = true;
            else
                isSingle = false;
            end
            if maxCorr(cell)> spatial_correlation_threshold
                isOverlapping = true;
            else
                isOverlapping = false;
            end
            
            if isOverlapping && ~isSingle
                ts_count = ts_count + 1;
            elseif isOverlapping && isSingle
                fs_count = fs_count + 1;
            elseif ~isOverlapping && ~isSingle
                td_count = td_count + 1;
            elseif ~isOverlapping && isSingle
                fd_count = fd_count + 1;
            end
        end
        fprintf('Session 1 Precision --> Detection: %f, Splitting: %f\n', ...
            td_count/(td_count+fd_count), ts_count/(ts_count+fs_count));
    end
    
end


end

