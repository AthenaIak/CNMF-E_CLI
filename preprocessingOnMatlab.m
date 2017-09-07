function [] = preprocessing(set_parameters)

TIFF_MAX_FRAMES = 1349;

run_setup;
clear CNMF_dir;

% set_parameters='~/tu/athina/Data/analyzed/parameters/parameters_an014_preproc';
run (set_parameters);
numMovies = length(movieFiles);

for m=1:numMovies
    numFiles = length(movieFiles{m});
    
    curr_id = 1;
    for f = 1:numFiles
        in_nam = fullfile(inDir,movieFiles{m}{f});
        
        %load movie file
        fprintf('Loading file %d/%d from movie %d/%d...\n', ...
            f,numFiles,m,numMovies); tic;
        tmpY = bigread2(in_nam);
        fprintf('Time cost of loading images: %.2f seconds\n', toc);
        
        %replace bad frames with weighted averages of previous and next
        %good frame.
        for badSet=1:size(badFrames{m}{f},1)
            prevFrame = badFrames{m}{f}{badSet}(1) - 1;
            nextFrame = badFrames{m}{f}{badSet}(end) + 1;
            numBad = size(badFrames{m}{f}{badSet},2);
            
            weights = linspace(1,0,numBad+2);
            weights = weights(2:end-1);
            for bf = 1:numBad
                tmpY(:,:,prevFrame+bf) = tmpY(:,:,prevFrame)*weights(bf)+...
                    tmpY(:,:,nextFrame)*(1-weights(bf));
            end
        end
        disp('Replaced bad frames.');
        
        numFrames = size(tmpY,3);
        Y(:,:,curr_id:curr_id+numFrames-1) = tmpY;
        clear tmpY;
        curr_id = curr_id + numFrames;
    end
    clear curr_id in_nam;
    clear badSet prevFrame nextFrame numBad weights bf;
    
    
    %temporally downsample
    downY = Y(:,:,1:down_factor:end);
    clear Y;
    
    %calculate indexes of replaced frames
    %first calculate the change resulting from file concatenation
    curr_id=1;
    for f=1:numFiles
        for badSet=1:size(badFrames{m}{f},1)
            vals = badFrames{m}{f}{badSet};
            numFr = size(vals,2);
            replaced(curr_id:curr_id+numFr-1) = vals+TIFF_MAX_FRAMES*(f-1);
            curr_id=curr_id+numFr;
        end
    end
    clear vals badSet f numFr;
    
    %now calculate the change resulting from downsampling
    replacedInUse = replaced(mod(replaced,down_factor)==1);
    newIdxs = floor(replacedInUse/down_factor)+1;
    clear replaced replacedInUse;
    
    disp('Indexes of frames with artificial data:');
    disp(newIdxs);
    
    %save as tiff
    [~,rec_nam,~] = fileparts(movieFiles{m}{1});
    curr_id=1;
    for img=1:min(size(downY,3))
        if mod(img,TIFF_MAX_FRAMES)==1
            outputFileName = fullfile(outDir,sprintf('pp_%s-%d.tif',rec_nam,curr_id));
            curr_id = curr_id+1;
        end
        imwrite(downY(:, :, img), outputFileName, ...
            'WriteMode', 'append', 'Compression','none');
    end
    
    %clear unnecessary memory
    pack;
end

clear numMovies

