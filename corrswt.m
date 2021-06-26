function [r,p,ci,n] = corrswt(x,y,varargin)

argsin = varargin;

argsin = setdefault(argsin,'rows','pairwise');
argsin = setdefault(argsin,'type','spearman');

[r,p] = corr(x,y,argsin{:});

for i = 1:size(x,2)
    for ii = 1:size(y,2)
        n(i,ii) = sum(~any(isnan([x(:,i) y(:,ii)]),2));
        ci(i,ii,:) = corrci(r(i,ii),n(i,ii));
    end
end