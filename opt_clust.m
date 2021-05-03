function clstruct = opt_clust(clfunc,X,kvals,nreps,dist,varargin)

argsin = varargin;

if strcmpi(dist,'sqeuclidean')
    dists = pdist(X,'euclidean');
else
    dists = pdist(X,dist);
end

clstruct.dists = squareform(dists);

for k = kvals
    [cl,ctrs] = consensus_clust(clfunc,X,k,nreps,'Distance',dist,argsin{:});
    sil = silhouette(X,cl,dist);
    
    for i = 1:nreps
        randset = randperm(size(X,2),round(size(X,2)/2));
        
        train = X(:,randset); test = X(:,randset);
        cltrain = clfunc(train,k,argsin{:}); cltest = clfunc(test,k,argsin{:});
        trainadj = cltrain==cltrain'; testadj = cltest==cltest';
        stab(i) = dice(belowDiag(trainadj),belowDiag(testadj));
    end
    clstruct.cl(:,k) = cl; clstruct.ctrs{k} = ctrs;
    clstruct.sil(:,k) = sil; clstruct.silscore(k) = nanmean(sil);
    clstruct.stab(:,k) = stab; clstruct.stabscore(k) = nanmean(stab);
end