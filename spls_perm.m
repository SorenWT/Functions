function mdl = spls_perm(X,Y,ncomps,nperms,varargin)

argsin = varargin;

% remove NaNs subjectwise
xnan = any(isnan(X),2); ynan = any(isnan(Y),2);
allnan = (xnan+ynan)>0;
X(allnan,:) = []; Y(allnan,:) = [];

mdl = spls_fit(X,Y,ncomps,argsin{:});

for i = 1:nperms
    permX = X(randperm(size(X,1)),:);
    permmdl(i) = spls_fit(permX,Y,ncomps,argsin{:});
end

permmdl = mergestructs(permmdl);

mdl.pperm = 1-nanmean(mdl.compcov>permmdl.compcov,1);