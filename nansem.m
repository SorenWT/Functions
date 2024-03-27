function [se] = nansem(X,w,dim)

se = nanstd(X,w,dim)./sum(~isnan(X),dim);