%% add path
CNMF_dir = fileparts(which('run_setup.m'));
addpath(sprintf('%s%sca_source_extraction%s', CNMF_dir, filesep, filesep));
addpath(genpath(sprintf('%s%sca_source_extraction%sutilities%s', CNMF_dir, filesep, filesep, filesep)));
addpath(genpath(sprintf('%s%sca_source_extraction%sendoscope', CNMF_dir, filesep, filesep, filesep)));
addpath(sprintf('%s%sGUI', CNMF_dir, filesep));
addpath(sprintf('%s%sGUI%sgui_callbacks', CNMF_dir, filesep, filesep));
addpath(sprintf('%s%sGUI%smodules', CNMF_dir, filesep, filesep));
addpath(sprintf('%s%soasis', CNMF_dir, filesep));
addpath(sprintf('%s%scnmfe_scripts', CNMF_dir, filesep));
addpath(sprintf('%s%sCLI', CNMF_dir, filesep));

%% setup cvx
if isempty(which('cvx_begin.m'))
    if ~exist('cvx', 'dir')
        %install cvx
        if ismac
            cvx_url = 'http://web.cvxr.com/cvx/cvx-maci64.zip';
        elseif isunix
            cvx_url = 'http://web.cvxr.com/cvx/cvx-a64.zip';
        elseif ispc
            cvx_url = 'http://web.cvxr.com/cvx/cvx-w64.zip';
        else
            fprints('Your platform is not supported by CVX\n');
            return;
        end
        fprintf('Downloading CVX...\n');
        unzip(cvx_url, CNMF_dir);
    end
    run(sprintf('%s%scvx%scvx_setup', CNMF_dir, filesep, filesep));
end
%% save path
%savepath();
