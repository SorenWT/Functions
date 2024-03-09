function [lassomdl] = lasso_perm(X,Y,nperms,varargin)

argsin = varargin;
argsin = setdefault(argsin,'CV',5);
argsin = setdefault(argsin,'NumLambda',100);
argsin = setdefault(argsin,'RelTol',1e-4);
argsin = setdefault(argsin,'Alpha',0.001);
argsin = setdefault(argsin,'pca','off');

xnan = any(isnan(X),2); ynan = isnan(Y);
allnan = (xnan+ynan)>0; goodindx = find(~allnan);
X(allnan,:) = []; Y(allnan) = [];

if EasyParse(argsin,'pca','on')
    pc_reg = pca_swt(X,'compsel','kaiser','nperms',10); 
    oldX = X; 
    X = pc_reg.comps(:,1:pc_reg.ncomps);
end

[lassoweights,stats] = lasso(X,Y,'CV',EasyParse(argsin,'CV'),...
    'NumLambda',EasyParse(argsin,'NumLambda'),'RelTol',EasyParse(argsin,'RelTol'),'Alpha',EasyParse(argsin,'Alpha'));

[~,stats.IndexMinMSE] = min(stats.MSE(any(lassoweights,1)));
lassoweights = lassoweights(:,stats.IndexMinMSE);
pred = nancenter(X,1)*lassoweights;
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
    pred_shuf = nancenter(X,1)*lassoweightsshuf;
    r_shuf(q) = corr(pred_shuf,Y(randindx),'rows','pairwise');
end
coeffs = lassoweights;
stats.pperm = 1-nanmean(r_obs>r_shuf);

lassomdl.coeffs = lassoweights; lassomdl.stats = stats; 
lassomdl.pred = NaN(length(allnan),1); lassomdl.pred(goodindx) = pred;

lassomdl.r = r_obs;
lassomdl.rperm = r_shuf; lassomdl.pperm = stats.pperm;

if EasyParse(argsin,'pca','on')
    lassomdl.pccoeffs = lassomdl.coeffs;
    lassomdl.coeffs = pc_reg.weights(:,1:pc_reg.ncomps)*lassomdl.pccoeffs;
end
