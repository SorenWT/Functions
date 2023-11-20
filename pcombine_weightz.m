function [p,zcomb] = pcombine_weightz(pvals,nvals)

z = norminv(1-pvals);

zcomb = sum(sqrt(nvals).*z)/sqrt(sum(nvals));

p = 1-normcdf(zcomb);