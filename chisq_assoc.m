function [p,chisq,conttbl] = chisq_assoc(x,y,ndisc)

if nargin < 3
ndisc = 2;
end

xdisc = discretize(x,prctile(x,linspace(0,100,ndisc+1)));
ydisc = discretize(y,prctile(y,linspace(0,100,ndisc+1)));

[conttbl,chisq,p] = crosstab(xdisc,ydisc);