function [canonmdl] = canoncorr_perm(X,Y,nperm,param,varargin)
% param is the parameter for the permutation test - default = chisq

if nargin < 4
    param = 'chisq';
end

xnan = any(isnan(X),2); ynan = any(isnan(Y),2);
allnan = (xnan+ynan)>0;
X(allnan,:) = []; Y(allnan,:) = [];

[A,B,R,U,V,stats] = canoncorr(X,Y);

for i = 1:nperm
    permX = X(randperm(size(X,1)),:);
    [~,~,Rperm(:,i),~,~,statperm(i)] = canoncorr(permX,Y);
end

statperm = mergestructs(statperm);
stats.pperm = 1-nanmean(R'>Rperm,2);

canonmdl = struct; 
canonmdl.pperm = horz(stats.pperm); 
canonmdl.A = A; canonmdl.B = B; canonmdl.R = R; canonmdl.U = U; canonmdl.V = V; canonmdl.stats = stats;
canonmdl.Rperm = Rperm; canonmdl.statperm = statperm;
%stats.pperm = 1-value_prctile(Rperm,R);

%stats.pperm = 1-value_prctile(statperm.(param),stats.(param));

%stats.pperm(stats.pperm>0.5) = 1-stats.pperm(stats.pperm>0.5);

%stats.pperm = 2*stats.pperm; %two-sided test