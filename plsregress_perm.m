function plsmdl = plsregress_perm(X,Y,ncomp,nperm,varargin)

argsin = varargin;

argsin = setdefault(argsin,'nboot',nperm);
argsin = setdefault(argsin,'stratify',ones(size(X,1),1));
argsin = setdefault(argsin,'discrim',false); 

plsmdl = struct;

tmp = whos('X');
if tmp.bytes<1e8
   plsmdl.X = X;
end

tmp = whos('Y');
if tmp.bytes<1e8
   plsmdl.Y = Y;
end

xnan = any(isnan(X),2); ynan = any(isnan(Y),2);
allnan = (xnan+ynan)>0; goodindx = find(~allnan);
X(allnan,:) = []; Y(allnan,:) = [];

if size(X,1) >= 100
    argsin = setdefault(argsin,'permmethod','mholdout');
else
    argsin = setdefault(argsin,'permmethod','resub');
end
argsin = setdefault(argsin,'stratify',ones(size(X,1),1));

nboot = EasyParse(argsin,'nboot');
permmethod = EasyParse(argsin,'permmethod');
strat = EasyParse(argsin,'stratify');
strat(allnan) = [];

discrim = EasyParse(argsin,'discrim');

if discrim && (islogical(Y) | all(unique(Y(~isnan(Y)))==[0;1]))
   Y = Y*2-1;
end

[XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress_swt(X,Y,ncomp);
r = corr(XS,YS); r = r(find(eye(size(r)))); r = r';
stats.r = r; stats.mse = mse;
plsmdl.explained = stats.explained;

predY = [ones(size(X,1),1) X]*beta;

if discrim
    [ctab,chi2,p]=crosstab(predY>0,Y);
    plsmdl.ctab = ctab; 
    plsmdl.cmetrics.acc = sum(sum(ctab.*eye(size(ctab))))/sum(sum(ctab));
    plsmdl.cmetrics.sens = ctab(2,2)/sum(ctab(:,2));
    plsmdl.cmetrics.spec = ctab(1,1)/sum(ctab(:,1));
    plsmdl.cmetrics.balacc = mean([plsmdl.cmetrics.sens,plsmdl.cmetrics.spec]);
end


switch permmethod
    case 'mholdout'
        for q = 1:nboot
            cvp = cvpartition(strat,'Holdout',0.33);
            Xtrain = X(training(cvp),:); Ytrain = Y(training(cvp),:);
            Xtest = X(except(1:size(X,1),find(training(cvp))),:); Ytest = Y(except(1:size(Y,1),find(training(cvp))),:);
            
            % this is wrong! XL and YL aren't comparable to the
            % coefficients from matlab's PCA function
            %[XLhold(:,:,q),YLhold(:,:,q)] = plsregress_swt(Xtrain,Ytrain,ncomp);
            
            % fixed version
            [~,YLhold,XShold,~,betahold,~,~,stats] = plsregress_swt(Xtrain,Ytrain,ncomp);
            XWhold = stats.W;
                        
            XStest = nancenter(Xtest,1)*XWhold;
            
            % should this be centering across dimension 2?
            YStest = nancenter(Ytest,1)*YLhold;
            for qq = 2:size(XStest,2)
                for qqq = 1:qq-1
                    YStest(:,qq) = projout(YStest(:,qq),XStest(:,qqq));
                end
            end
            
            if discrim
                YpredTest = [ones(size(Xtest,1),1) Xtest]*betahold; 
                [ctab] = crosstab(YpredTest>0,Ytest);
                holdsens(q) = ctab(2,2)/sum(ctab(:,2));
                holdspec(q) = ctab(1,1)/sum(ctab(:,1));
                holdbalacc(q) = mean([holdsens(q),holdspec(q)]);
                
                holdacc(q) = sum(sum(ctab.*eye(size(ctab))))/sum(sum(ctab));
            end
            
            
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
            %tmp = corr(nancenter(Xtest,1)*XLhold,nancenter(Ytest,1)*YLhold);
            tmp = corr(XStest,YStest);
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
        if discrim
            %stats.allholdacc = holdacc;
            %stats.meanholdacc = mean(holdacc);

            plsmdl.holdmetrics.acc = mean(holdacc);
            plsmdl.holdmetrics.sens = mean(holdsens);
            plsmdl.holdmetrics.spec = mean(holdspec);
            plsmdl.holdmetrics.balacc = mean(holdbalacc);
        end
        %plsmdl.holdr = allholdperf;
        
        %     meancvres = mean(cvres,2);
        %     meancvperm = squeeze(mean(cvperm,2));
        %
        %     stats.cvres = cvres; stats.cvperm = cvperm;
        %
        %     stats.pperm = 1-nanmean(meancvres>meancvperm');
        
    case 'cvperm'
        % permutation test on cross-validation scores
        
        for i = 1:nperm
            permX = X(randperm(size(X,1)),:);
            permY = Y(randperm(size(Y,1)),:);
            
            [~,~,~,~,~,~,mseperm(:,:,i)] = plsregress_swt(permX,permY,ncomp,'cv',5);
        end
        
        
        stats.mseperm = mseperm;
        stats.pperm = 1-nanmean(mse(2,:)'<squeeze(mseperm(2,:,:)),2);
        stats.pperm(1) = [];
        
    case {'orig','resub'}
        for i = 1:nperm
            permX = X(randperm(size(X,1)),:);
            permY = Y(randperm(size(Y,1)),:);
            
            [~,~,XSperm,YSperm,betaperm,~,mseperm(:,:,i),permstat] = plsregress_swt(permX,permY,ncomp);
            sings_perm(i,:) = permstat.sings;
            tmp = corr(XSperm,YSperm); rperm(i,:) = tmp(find(eye(size(tmp))));
            
            if discrim
                Ypredperm = [ones(size(permX,1),1) permX]*betaperm; 
                [tmptbl] = crosstab(Ypredperm>0,permY);
                permacc(i) = sum(sum(tmptbl.*eye(size(tmptbl))))/sum(sum(tmptbl));
            end
        end
        stats.sings_perm = sings_perm; stats.rperm = rperm;
        stats.pperm = 1-nanmean(stats.sings>sings_perm); % one-tailed test
        %stats.pperm = 1-nanmean(stats.r>stats.rperm);
        
        if discrim
             stats.acc_perm = permacc; 
             stats.pperm_acc = 1-nanmean(plsmdl.cmetrics.acc>permacc);
        end
end

plsmdl.XL = XL; plsmdl.YL = YL;
plsmdl.XS = NaN(length(allnan),size(XS,2)); plsmdl.XS(goodindx,:) = XS;
plsmdl.YS = NaN(length(allnan),size(YS,2)); plsmdl.YS(goodindx,:) = YS;
plsmdl.predY = NaN(length(allnan),size(Y,2)); plsmdl.predY(goodindx,:) = predY;
plsmdl.beta = beta; plsmdl.pctvar = pctvar;
plsmdl.mse = mse; plsmdl.stats = stats;
plsmdl.r = r; plsmdl.pperm = stats.pperm;
%plsmdl.fdr = fdr(plsmdl.pperm);

plsmdl.pperm = horz(plsmdl.pperm);
plsmdl.pperm = plsmdl.pperm+1/nperm;
if isfield(stats,'sings_perm')
plsmdl.sings_perm = sings_perm; plsmdl.rperm = rperm;
end
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
