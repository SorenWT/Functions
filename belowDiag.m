function [outmat,maskedmat,mask] = belowDiag(inmat)

mask = tril(ones(size(inmat)))-eye(size(inmat));
if ~iscell(inmat)
    maskedmat = inmat.*mask;
end
outmat = reshape(inmat(find(mask)),[],1);