function [dev,propdev] = uniformdev(condseq)

    if iscell(condseq)
       condseq = factor2num(condseq); 
    end
    
    counts = histcounts(condseq,[0 unique(condseq)+0.5]);
    
    %counts = wordCloudCounts(condseq);
    
    %counts = counts.Count;
    
    expected = ceil(sum(counts))./length(counts);
    
    dev = max(ceil(abs(counts-expected)));
    
    propdev = mean(abs(counts-expected)./expected);

end