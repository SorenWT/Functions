function X = nancenter(X,dim)

X = X-nanmean(X,dim);