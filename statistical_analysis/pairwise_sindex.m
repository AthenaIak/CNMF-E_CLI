function [ sindex ] = pairwise_sindex( spikes )
%PAIRWISE_SINDEX Computes the s-index for all possible pairs of neurons.
    nNeurons=size(spikes,1);
    nComb = sum(1:nNeurons-1);
    sindex=zeros(1,nComb,'double'); k=1;
    
    for i=1:nNeurons-1
        for j=i+1:nNeurons
            ca = spikes(i,:); cb = spikes(j,:);
            sindex(k) =  dot(ca,cb) / (norm(ca)^2+norm(cb)^2) * 2;
            k=k+1;
        end
    end


end

