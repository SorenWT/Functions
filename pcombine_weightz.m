function p = pcombine_weightz(pvals,nvals)

z = norminv(1-pvals);

p = 1-normcdf(sum(sqrt(nvals).*z)/sqrt(sum(nvals)));