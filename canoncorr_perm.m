function ccmdl = canoncorr_perm(X,Y,nperm,varargin)

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

[XL,YL,r,XS,YS,stats] = canoncorr(X,Y);
%r = corr(XS,YS); r = r(find(eye(size(r)))); r = r';
stats.r = r;

ccmdl = struct;

switch permmethod
    case 'cv'
        for q = 1:nboot
            cvp = cvpartition(strat,'Holdout',0.33);
            Xtrain = X(training(cvp),:); Ytrain = Y(training(cvp),:);
            Xtest = X(except(1:size(X,1),find(training(cvp))),:); Ytest = Y(except(1:size(Y,1),find(training(cvp))),:);
            [XLhold(:,:,q),YLhold(:,:,q)] = canoncorr(Xtrain,Ytrain);
            
            %[~,~,trans_x] = procrustes(XL,XLhold(:,:,q),'Scaling',false);
            %[~,~,trans_y] = procrustes(YL,YLhold(:,:,q),'Scaling',false);
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
            
%             T = (abs(trans_x.T)+abs(trans_y.T))./2;
%             for i = 1:size(T,2)
%                 T(find(abs(T(:,i)) == max(abs(T(:,i)))),i) = 1;
%                 T(find(T(:,i) ~= 1),i) = 0;
%             end
%             
%             XLhold_rot(:,:,q) = XLhold(:,:,q)*T; % reorder components, but don't change them otherwise
%             YLhold_rot(:,:,q) = YLhold(:,:,q)*T;
%             allT(:,:,q) = T;
            
            
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
        
        ccmdl.holdperf = mean(perf,2)';
        
        %     meancvres = mean(cvres,2);
        %     meancvperm = squeeze(mean(cvperm,2));
        %
        %     stats.cvres = cvres; stats.cvperm = cvperm;
        %
        %     stats.pperm = 1-nanmean(meancvres>meancvperm');
        
    case {'orig','resub'}
        for i = 1:nperm
            permX = X(randperm(size(X,1)),:);
            permY = Y(randperm(size(Y,1)),:);
            
            [~,~,rperm(i,:),XSperm,YSperm,stats_perm] = canoncorr(permX,permY);
            %sings_perm(i,:) = permstat.sings;
            tmp = corr(XSperm,YSperm); rperm(i,:) = tmp(find(eye(size(tmp))));
            F_perm(i,:) = stats_perm.F;
        end
        stats.F_perm = F_perm; 
        stats.rperm = rperm;
        stats.pperm = 1-nanmean(stats.F>stats.F_perm); % one-tailed test
end

ccmdl.XL = XL; ccmdl.YL = YL;
ccmdl.XS = NaN(length(allnan),size(XS,2)); ccmdl.XS(goodindx,:) = XS;
ccmdl.YS = NaN(length(allnan),size(YS,2)); ccmdl.YS(goodindx,:) = YS;
%plsmdl.beta = beta; plsmdl.pctvar = pctvar;
%plsmdl.mse = mse; 
ccmdl.stats = stats;
ccmdl.r = r; ccmdl.pperm = stats.pperm;
ccmdl.fdr = fdr(ccmdl.pperm);

ccmdl.Xloads = corr(X,XS); ccmdl.Yloads = corr(Y,YS);

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
    
    ccmdl.Xloads_boot = Xloads_boot; ccmdl.Yloads_boot = Yloads_boot;
    
    ccmdl.Xloads_bootz = ccmdl.Xloads./std(Xloads_boot,[],3);
    ccmdl.Xloads_bootp = ztop(ccmdl.Xloads_bootz);
    
    ccmdl.Yloads_bootz = ccmdl.Yloads./std(Yloads_boot,[],3);
    ccmdl.Yloads_bootp = ztop(ccmdl.Yloads_bootz);
end
end

function r = plspredict(Xtrain,Ytrain,Xtest,Ytest,ncomp)

[XL,YL] = plsregress(Xtrain,Ytrain,ncomp);

r = corr(Xtest*XL,Ytest*YL);
r = r(find(eye(size(r))));

end
