function [XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress_perm(X,Y,ncomp,nperm)
% param is the parameter for the permutation test - default = chisq


[XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress_swt(X,Y,ncomp);
r = corr(XS,YS); r = r(find(eye(size(r)))); r = r';

for i = 1:nperm
    permX = X(randperm(size(X,1)),:);
    [~,~,XSperm,YSperm,~,~,~,permstat] = plsregress_swt(permX,Y,ncomp);
    sings_perm(i,:) = permstat.sings;
    %tmp = corr(XSperm,YSperm); rperm(i,:) = tmp(find(eye(size(tmp))));
end

%stats.pperm = 1-value_prctile(rperm,r);
stats.pperm = 1-value_prctile(sings_perm,stats.sings);

%stats.pperm(stats.pperm>0.5) = 1-stats.pperm(stats.pperm>0.5);

%stats.pperm = 2*stats.pperm; %two-sided test