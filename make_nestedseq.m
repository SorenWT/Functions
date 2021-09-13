function seqnest = make_nestedseq(condnumbers,nreps)

nel = prod(condnumbers);

thislength = nel;
for i = 1:length(condnumbers)
    thislength = thislength/condnumbers(i);
    seqnest(i,:) = repmat(elrepmat(1:condnumbers(i),thislength),1,nel/(thislength*condnumbers(i)));
end

seqnest = repmat(seqnest,1,nreps);