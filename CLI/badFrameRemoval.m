function[] = badFrameRemoval(movieFiles, badFrames)

% old code (may be used in future development)

% load all files to tmpY and concatenate to Y
    
%         %replace bad frames with weighted averages of previous and next
%         %good frame.
%         %make sure it is not empty 
%         if ~isempty(badFrames{m}{f})
%             for badSet=1:size(badFrames{m}{f},1)
%                 prevFrame = badFrames{m}{f}{badSet}(1) - 1;
%                 nextFrame = badFrames{m}{f}{badSet}(end) + 1;
%                 numBad = size(badFrames{m}{f}{badSet},2);
% 
%                 weights = linspace(1,0,numBad+2);
%                 weights = weights(2:end-1);
%                 for bf = 1:numBad
%                     tmpY(:,:,prevFrame+bf) = tmpY(:,:,prevFrame)*weights(bf)+...
%                         tmpY(:,:,nextFrame)*(1-weights(bf));
%                 end
%             end
%             disp('Replaced bad frames.');
%         else
%             disp('No bad frames defined.');
%         end
        


%clear badSet prevFrame nextFrame numBad weights bf;

%calculate indexes of replaced frames
%first calculate the change resulting from file concatenation
% curr_id=1;
% for f=1:numFiles
% 	for badSet=1:size(badFrames{m}{f},1)
%     	vals = badFrames{m}{f}{badSet};
%         numFr = size(vals,2);
%         replaced(curr_id:curr_id+numFr-1) = vals+TIFF_MAX_FRAMES*(f-1);
%         curr_id=curr_id+numFr;
%     end
% end
% clear vals badSet f numFr;
%     
% %now calculate the change resulting from downsampling
% if exist('replacedInUse','var') && (~isempty(replacedInUse))
% 	replacedInUse = replaced(mod(replaced,down_factor)==1);
%     newIdxs = floor(replacedInUse/down_factor)+1;
%     clear replaced replacedInUse;
%     
%     disp('Indexes of frames with artificial data:');
%     disp(newIdxs);
% else
% disp('No frames were replaced in any file.');    
% end