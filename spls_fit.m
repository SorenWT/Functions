function mdl = spls_fit(X,Y,ncomps,varargin)
% fits a sparse PLS model based on Monteiro et al. (2016)
% uses the spls code from Monteiro et al. but implements their recommended
% cross-validation procedure to optimize cu and cv

argsin = varargin;
argsin = setdefault(argsin,'NumLambda',5);
argsin = setdefault(argsin,'Folds',5);

nl = EasyParse(argsin,'NumLambda');
kfolds = EasyParse(argsin,'Folds');

% remove NaNs subjectwise
xnan = any(isnan(X),2); ynan = any(isnan(Y),2);
allnan = (xnan+ynan)>0;
X(allnan,:) = []; Y(allnan,:) = [];

% center X and Y
X = X - mean(X,1); Y = Y - mean(Y,1);


thisX = X; thisY = Y;

cu = linspace(1,sqrt(size(X,2)),nl);
cv = linspace(1,sqrt(size(Y,2)),nl);

%cu = cu(2:end-1); cv = cv(2:end-1); % don't want the maximum sparsity that leaves only one variable, or the minimum which leaves everything

for comp = 1:ncomps    
    % find optimal cu and cv
    %cval = cvpartition(size(X,1),'KFold',kfolds);
    for q = 1:kfolds
        cval = cvpartition(size(X,1),'Holdout',0.2);
        testu = zeros(size(X,2),nl,nl); testv = zeros(size(Y,2),nl,nl);
        for i = 1:nl
            for ii = 1:nl
                trn = training(cval); tst = test(cval);
                [testu(:,i,ii),testv(:,i,ii)] = spls(thisX(trn,:),thisY(trn,:),cu(i),cv(ii),1e-4);
                testxs = thisX(tst,:)*testu(:,i,ii); testys = thisY(tst,:)*testv(:,i,ii);
                corrvals(i,ii,q) = corr(testxs,testys);
            end
        end
    end
    meancorrvals = mean(corrvals,3);
    [bestcu,bestcv] = ind2sub(size(meancorrvals),find(meancorrvals==max(meancorrvals),1));
    
    % make the final component
    [u(:,comp),v(:,comp)] = spls(thisX,thisY,cu(bestcu),cv(bestcv),1e-4);
    xs(:,comp) = thisX*u(:,comp); ys(:,comp) = thisY*v(:,comp);
    compcorr(comp) = corr(xs(:,comp),ys(:,comp));
    tmp = cov([xs(:,comp),ys(:,comp)]);
    compcov(comp) = tmp(1,2);
    optcu(comp) = cu(bestcu); optcv(comp) = cv(bestcv);
    
    [thisX,thisY] = proj_def(thisX,thisY,u(:,comp),v(:,comp));
end

mdl.u = u; mdl.v = v; mdl.xs = xs; mdl.ys = ys;
mdl.compcorr = compcorr; mdl.compcov = compcov;
mdl.optcu = optcu; mdl.optcv = optcv;


