function [outmat,maskedmat,mask] = belowDiag(inmat)

mask = tril(ones(size(inmat)))-eye(size(inmat));
maskedmat = inmat.*mask;
outmat = reshape(inmat(find(mask)),[],1);