function plsmdl = plsregress_perm(X,Y,ncomp,nperm,varargin)

argsin = varargin;

argsin = setdefault(argsin,'nboot',nperm);
argsin = setdefault(argsin,'stratify',ones(size(X,1),1));



xnan = any(isnan(X),2); ynan = any(isnan(Y),2);
allnan = (xnan+ynan)>0; goodindx = find(~allnan);
X(allnan,:) = []; Y(allnan,:) = [];

if size(X,1) >= 100
    argsin = setdefault(argsin,'permmethod','cv');
else
    argsin = setdefault(argsin,'permmethod','resub');
end
argsin = setdefault(argsin,'stratify',ones(size(X,1),1));

nboot = EasyParse(argsin,'nboot');
permmethod = EasyParse(argsin,'permmethod');
strat = EasyParse(argsin,'stratify');
strat(allnan) = [];

[XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress_swt(X,Y,ncomp);
r = corr(XS,YS); r = r(find(eye(size(r)))); r = r';
stats.r = r;

plsmdl = struct;

if strcmpi(permmethod,'cv')    
    for q = 1:nboot
        cvp = cvpartition(strat,'Holdout',0.25);
        Xtrain = X(training(cvp),:); Ytrain = Y(training(cvp),:);
        Xtest = X(except(1:size(X,1),find(training(cvp))),:); Ytest = Y(except(1:size(Y,1),find(training(cvp))),:); 
        [XLhold(:,:,q),YLhold(:,:,q)] = plsregress_swt(Xtrain,Ytrain,ncomp);
        
        [~,~,trans_x] = procrustes(XL,XLhold(:,:,q),'Scaling',false);
         [~,~,trans_y] = procrustes(YL,YLhold(:,:,q),'Scaling',false);
%         for i = 1:size(trans_x.T,2)
%            trans_x.T(find(abs(trans_x.T(:,i)) == max(abs(trans_x.T(:,i)))),i) = 1;
%            trans_x.T(find(trans_x.T(:,i) ~= 1),i) = 0;
%         end
%         for i = 1:size(trans_y.T,2)
%            trans_y.T(find(abs(trans_y.T(:,i)) == max(abs(trans_y.T(:,i)))),i) = 1;
%            trans_y.T(find(trans_y.T(:,i) ~= 1),i) = 0;
%         end
%         XLhold_rot(:,:,q) = XLhold(:,:,q)*trans_x.T; % reorder components, but don't change them otherwise
%         YLhold_rot(:,:,q) = YLhold(:,:,q)*trans_x.T; % reorder components, but don't change them otherwise
%         
        
        T = (abs(trans_x.T)+abs(trans_y.T))./2;
        for i = 1:size(T,2)
           T(find(abs(T(:,i)) == max(abs(T(:,i)))),i) = 1;
           T(find(T(:,i) ~= 1),i) = 0;
        end
        
        XLhold_rot(:,:,q) = XLhold(:,:,q)*T; % reorder components, but don't change them otherwise
        YLhold_rot(:,:,q) = YLhold(:,:,q)*T;
        allT(:,:,q) = T;

        
        %[~,XLhold_rot(:,:,q)] = procrustes(XLhold,XLhold,'Scaling',false);
        %[~,YLhold_rot(:,:,q)] = procrustes(YLhold,YLhold,'Scaling',false);
        
        %tmp = corr(Xtest*XLhold_rot(:,:,q),Ytest*YLhold_rot(:,:,q));
        tmp = corr(Xtest*XLhold(:,:,q),Ytest*YLhold(:,:,q));
        perf(:,q) = tmp(find(eye(size(tmp))));
        %cvres(:,q) = crossval(@(xtr,ytr,xts,yts)plspredict(xtr,ytr,xts,yts,ncomp),X,Y,'Holdout',0.25);
        
%         for i = 1:nperm
%             permX = X(randperm(size(X,1)),:);
%             permY = Y(randperm(size(Y,1)),:);
%             cvperm(:,q,i) = crossval(@(xtr,ytr,xts,yts)plspredict(xtr,ytr,xts,yts,ncomp),permX,permY,'Holdout',0.25);
%         end
    end
    
    stats.pperm = (1-nanmean(perf>0,2))';
    stats.allholdperf = perf;
    stats.meanholdperf = mean(perf,2)';
    
    plsmdl.holdperf = mean(perf,2)';
    
%     meancvres = mean(cvres,2);
%     meancvperm = squeeze(mean(cvperm,2));
%     
%     stats.cvres = cvres; stats.cvperm = cvperm;
%     
%     stats.pperm = 1-nanmean(meancvres>meancvperm');
    
elseif strcmpi(permmethod,'orig') | strcmpi(permmethod,'resub')
    
    for i = 1:nperm
        permX = X(randperm(size(X,1)),:);
        permY = Y(randperm(size(Y,1)),:);
        
        [~,~,XSperm,YSperm,~,~,mseperm(:,:,i),permstat] = plsregress_swt(permX,permY,ncomp);
        sings_perm(i,:) = permstat.sings;
        tmp = corr(XSperm,YSperm); rperm(i,:) = tmp(find(eye(size(tmp))));
    end
    stats.sings_perm = sings_perm; stats.rperm = rperm; 
    stats.pperm = 1-nanmean(stats.sings>sings_perm); % one-tailed test
end

plsmdl.XL = XL; plsmdl.YL = YL;
plsmdl.XS = NaN(length(allnan),size(XS,2)); plsmdl.XS(goodindx,:) = XS;
plsmdl.YS = NaN(length(allnan),size(YS,2)); plsmdl.YS(goodindx,:) = YS;
plsmdl.beta = beta; plsmdl.pctvar = pctvar;
plsmdl.mse = mse; plsmdl.stats = stats;
plsmdl.r = r; plsmdl.pperm = stats.pperm;

plsmdl.Xloads = corr(X,XS); plsmdl.Yloads = corr(Y,YS);

if nboot > 0
    for i = 1:nboot
        bootindX = ceil(rand(size(X,1),1)*size(X,1));
        bootX = X(bootindX,:);
        bootindY = ceil(rand(size(Y,1),1)*size(Y,1));
        bootY = Y(bootindY,:);
        [~,~,XSboot,YSboot] = plsregress_swt(bootX,bootY,ncomp);
        [~,XSboot] = procrustes(XS(bootindX,:),XSboot);
        [~,YSboot] = procrustes(YS(bootindY,:),YSboot);
        Xloads_boot(:,:,i) = corr(X(bootindX,:),XSboot);
        Yloads_boot(:,:,i) = corr(Y(bootindY,:),YSboot);
    end
    
    plsmdl.Xloads_boot = Xloads_boot; plsmdl.Yloads_boot = Yloads_boot;
    
    plsmdl.Xloads_bootz = plsmdl.Xloads./std(Xloads_boot,[],3);
    plsmdl.Xloads_bootp = ztop(plsmdl.Xloads_bootz);
    
    plsmdl.Yloads_bootz = plsmdl.Yloads./std(Yloads_boot,[],3);
    plsmdl.Yloads_bootp = ztop(plsmdl.Yloads_bootz);
end
end

function r = plspredict(Xtrain,Ytrain,Xtest,Ytest,ncomp)

[XL,YL] = plsregress(Xtrain,Ytrain,ncomp);

r = corr(Xtest*XL,Ytest*YL);
r = r(find(eye(size(r))));

end
