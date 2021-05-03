function [lassomdl] = lasso_perm(X,Y,nperms,varargin)

argsin = varargin;
argsin = setdefault(argsin,'CV',5);
argsin = setdefault(argsin,'NumLambda',10);
argsin = setdefault(argsin,'RelTol',1e-3);
argsin = setdefault(argsin,'Alpha',1);

xnan = any(isnan(X),2); ynan = isnan(Y);
allnan = (xnan+ynan)>0; goodindx = find(~allnan);
X(allnan,:) = []; Y(allnan) = [];


[lassoweights,stats] = lasso(X,Y,'CV',EasyParse(argsin,'CV'),...
    'NumLambda',EasyParse(argsin,'NumLambda'),'RelTol',EasyParse(argsin,'RelTol'),'Alpha',EasyParse(argsin,'Alpha'));

[~,stats.IndexMinMSE] = min(stats.MSE(any(lassoweights,1)));
lassoweights = lassoweights(:,stats.IndexMinMSE);
pred = X*lassoweights;
r_obs = corr(pred,Y,'rows','pairwise');

%Xshuf = X;
for q = 1:nperms
    %         for qq = 1:size(Xcomps,2)
    %             Xcompsshuf(:,qq) = Xcomps(randperm(size(X,1)),qq);
    %         end
    randindx = randperm(size(X,1));
    [lassoweightsshuf,stats] = lasso(X,Y(randindx),'CV',EasyParse(argsin,'CV'),...
        'NumLambda',EasyParse(argsin,'NumLambda'),'RelTol',EasyParse(argsin,'RelTol'),'Alpha',EasyParse(argsin,'Alpha'));
    [~,stats.IndexMinMSE] = min(stats.MSE(any(lassoweightsshuf,1)));
    lassoweightsshuf = lassoweightsshuf(:,stats.IndexMinMSE);
    pred_shuf = X*lassoweightsshuf;
    r_shuf(q) = corr(pred_shuf,Y(randindx),'rows','pairwise');
end
coeffs = lassoweights;
stats.pperm = 1-nanmean(r_obs>r_shuf);

lassomdl.coeffs = lassoweights; lassomdl.stats = stats; 
lassomdl.pred = NaN(length(allnan),1); lassomdl.pred(goodindx) = pred;

lassomdl.r = r_obs;
lassomdl.rperm = r_shuf; lassomdl.pperm = stats.pperm;
