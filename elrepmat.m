function newX = elrepmat(X,nrep,dim)
% like repmat, but instead of tiling the whole matrix, it makes siz copies
% of the elements in each dimension
% only works for 1D vectors for now

if nargin < 3
   dim = 2;
end

newX = [];
X = horz(X);
for c = 1:size(X,dim)
    newX = cat(dim,newX,ones(1,nrep).*X(c));
end
