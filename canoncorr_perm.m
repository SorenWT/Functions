function [A,B,R,U,V,stats] = canoncorr_perm(X,Y,nperm,param)
% param is the parameter for the permutation test - default = chisq

if nargin < 4
    param = 'chisq';
end

[A,B,R,U,V,stats] = canoncorr(X,Y);

for i = 1:nperm
    permX = X(randperm(size(X,1)),:);
    [~,~,~,~,~,statperm(i)] = canoncorr(permX,Y);
end

statperm = mergestructs(statperm);

stats.pperm = 1-value_prctile(statperm.(param),stats.(param));

%stats.pperm(stats.pperm>0.5) = 1-stats.pperm(stats.pperm>0.5);

%stats.pperm = 2*stats.pperm; %two-sided test