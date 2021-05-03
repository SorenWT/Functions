function [XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress_boot(X,Y,ncomp,nperm)
% param is the parameter for the permutation test - default = chisq

xnan = any(isnan(X),2); ynan = any(isnan(Y),2);
allnan = (xnan+ynan)>0;
X(allnan,:) = []; Y(allnan,:) = [];

[XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress_swt(X,Y,ncomp);
r = corr(XS,YS); r = r(find(eye(size(r)))); r = r';

for i = 1:nperm
    permX = X(randperm(size(X,1)),:);
    permY = Y(randperm(size(Y,1)),:);
    [~,~,XSperm,YSperm,~,~,~,permstat] = plsregress_swt(permX,permY,ncomp);
    sings_perm(i,:) = permstat.sings;
    tmp = corr(XSperm,YSperm); rperm(i,:) = tmp(find(eye(size(tmp))));
end

stats.sings_perm = sings_perm; stats.rperm = rperm; stats.r = r;
stats.pperm = 1-nanmean(stats.sings>sings_perm); % one-tailed test