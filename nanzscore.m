function Z = nanzscore(X,W,dim)

Z = (X-nanmean(X,dim))./nanstd(X,W,dim);

