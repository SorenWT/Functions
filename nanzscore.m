function Z = nanzscore(X,W,dim)

if nargin < 2
   W = [];
end

if nargin < 3
   dim = 1; 
end

Z = (X-nanmean(X,dim))./nanstd(X,W,dim);

