function [p,zcomb] = pcombine_weightz(pvals,nvals)
% assumes two-tailed tests. Negative p values are treated as p values in
% the opposite direction to the direction of interest

newp = pvals/2;
newp(newp>0) = 1-newp(newp>0);
newp(newp<0) = -newp(newp<0);

z = norminv(newp);

zcomb = sum(sqrt(nvals).*z)/sqrt(sum(nvals));

if zcomb > 0
p = (1-normcdf(zcomb))*2;
else
   p = normcdf(zcomb)*2; 
end