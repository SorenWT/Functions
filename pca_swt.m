function pcamdl = pca_swt(X,varargin)

argsin = varargin;
argsin = setdefault(argsin,'nperms',1000);
argsin = setdefault(argsin,'do_boot',true);
argsin = setdefault(argsin,'rotate','none');
argsin = setdefault(argsin,'ncomps',NaN);
argsin = setdefault(argsin,'compsel','permutation');
argsin = setdefault(argsin,'kfold',5);

do_boot = EasyParse(argsin,'do_boot');
nperms = EasyParse(argsin,'nperms');
rotate = EasyParse(argsin,'rotate');
ncompsprior = EasyParse(argsin,'ncomps');
compsel = EasyParse(argsin,'compsel');
k = EasyParse(argsin,'kfold');

if any(any(isnan(X),2))
    nanindx = any(isnan(X),2);
    Xorig = X;
    X(nanindx,:) = [];
else
    Xorig = X;
    nanindx = zeros(size(X,1),1);
end

if all(all(isnan(X)))
    pcamdl.weights = NaN(size(X,2),min(size(X,2),size(X,1)));
    pcamdl.comps = NaN(size(X,1),size(X,2));
    pcamdl.explained = NaN(size(X,2),1);
    pcamdl.expl_perm = NaN(nperms,size(X,2));
    pcamdl.pperm = NaN(size(X,2),1);
    pcamdl.ncomps = NaN;
    pcamdl.kaiser = NaN;
    pcamdl.loads = NaN(size(pcamdl.weights));
    if do_boot
        pcamdl.loads_boot = NaN(size(X,2),size(X,2),nperms);
        pcamdl.loads_bootz = pcamdl.loads;
        pcamdl.loads_bootp = pcamdl.loads;
    end
else
    argsin = removeargs(argsin,{'do_boot','rotate','ncomps','compsel','kfold','nperms'});
    %     argsin(find(strcmpi(argsin,'do_boot')):find(strcmpi(argsin,'do_boot'))+1) = [];
    %     argsin(find(strcmpi(argsin,'rotate')):find(strcmpi(argsin,'rotate'))+1) = [];
    %     argsin(find(strcmpi(argsin,'ncomps')):find(strcmpi(argsin,'ncomps'))+1) = [];
    %     argsin(find(strcmpi(argsin,'compsel')):find(strcmpi(argsin,'ncomps'))+1) = [];
    %     argsin(find(strcmpi(argsin,'kfold')):find(strcmpi(argsin,'kfold'))+1) = [];
    %             argsin(find(strcmpi(argsin,'compsel')):find(strcmpi(argsin,'compsel'))+1) = [];
    
    
    X = nancenter(X,1);
    [weights,comps,~,~,explained] = pca(X,argsin{:});
    
    pcamdl.weights = weights;
    pcamdl.comps = NaN(size(Xorig,1),size(comps,2)); pcamdl.comps(~nanindx,:) = comps;
    
    pcamdl.explained = explained;
    
    switch compsel
        case 'permutation'
            for i = 1:nperms
                Xperm = X;
                for q = 1:size(X,2)
                    Xperm(:,q) = X(randperm(size(X,1)),q);
                end
                
                [~,~,~,~,expl_perm(:,i)] = pca(Xperm,argsin{:});
            end
            
            pcamdl.expl_perm = expl_perm;
            pcamdl.pperm = 1-nanmean(explained>expl_perm,2);
            pcamdl.ncomps = find(pcamdl.pperm>0.05,1)-1;
        case 'perm_flip'
            for i = 1:nperms
                R = rand(size(X)); R(R>0.5)=1; R(R<0.5)=-1;
                Xperm = R.*X;
%                 for q = 1:size(X,2)
%                     Xperm(:,q) = X(randperm(size(X,1)),q);
%                 end
                
                [~,~,~,~,expl_perm(:,i)] = pca(Xperm,argsin{:});
            end
            
            pcamdl.expl_perm = expl_perm;
            pcamdl.pperm = 1-nanmean(explained>expl_perm,2);
            pcamdl.ncomps = find(pcamdl.pperm>0.05,1)-1;
        case 'crossval'
            % if there's a lot of data, ALS is going to be super slow. 
            % Restrict the number of components (since we know 
            % cross-validation is going to give us a small number)
            if size(X,1)*size(X,2) > 1e6
                ncompsrestrict = find(cumsum(explained)>95,1);
            else
                ncompsrestrict = min(size(X))-1;
            end
            
            cv = cvpartition(repmat([1:size(X,1)]',size(X,2),1),'KFold',k);
            for i = 1:k
                trainmask = reshape(training(cv,i),size(X,1),size(X,2));
                testmask = reshape(test(cv,i),size(X,1),size(X,2));
                %[comps_train,weights_train] = nipals(X.*nanmask(trainmask),ncompsrestrict);
                [weights_train,comps_train] = pca(X.*nanmask(trainmask),'algorithm','als','NumComponents',ncompsrestrict);
                %[weights_test,scr_test] = pca(X(test(cv,i),:));
                for ii = 1:size(weights_train,2)
                    %scr_test = nancenter(X(test(cv,i),:),1)*weights_train(:,:);
                    X_pred = comps_train(:,1:ii)*weights_train(:,1:ii)';
                    resid{i,ii} = X.*testmask-X_pred.*testmask;
                    %resid{i,ii} = nancenter(X(test(cv,i),:),1)-X_pred;
                end
            end
            
            
            err = cellfun(@(d)sum(sum(d.^2)),resid);
            err = mean(err,1);
            
            [~,pcamdl.ncomps] = findpeaks(-err);
            
            if length(pcamdl.ncomps)>1
                warning('Multiple peaks found - check by manual inspection. Lowest number of components selected')
                pcamdl.ncomps = pcamdl.ncomps(1);
            end
            pcamdl.cverr = err; pcamdl.resid = resid;
        case 'kaiser'
            pcamdl.kaiser = find(explained<100/size(X,2),1)-1;
            if isempty(pcamdl.kaiser)
                pcamdl.kaiser = length(pcamdl.explained);
            end
            pcamdl.ncomps = pcamdl.kaiser;
    end
    
    if isnan(ncompsprior)
        ncompsprior = pcamdl.ncomps;
    end
    
    pcamdl.loads = corr(X,comps);
    
    if strcmpi(rotate,'none')
        if do_boot
            for i = 1:nperms
                resamp = ceil(rand(size(X,1),1)*size(X,1));
                [~,comps_resam] = pca(X(resamp,:),argsin{:});
                [~,comps_resam] = procrustes(comps(resamp,:),comps_resam);
                loads_resam(:,:,i) = corr(X(resamp,:),comps_resam);
            end
            
            pcamdl.loads_boot = loads_resam;
            pcamdl.loads_bootz = pcamdl.loads./std(loads_resam,[],3);
            pcamdl.loads_bootp = ztop(pcamdl.loads_bootz);
        end
    else
        ncompsrot = ncompsprior;
        pcamdl.rotated.weights = rotatefactors(pcamdl.weights(:,1:ncompsrot),'method',rotate);
        tmpcomps = (X-nanmean(X,1))*pcamdl.rotated.weights;
        pcamdl.rotated.comps = NaN(size(Xorig,1),size(tmpcomps,2)); pcamdl.rotated.comps(~nanindx,:) = tmpcomps;
        pcamdl.rotated.explained = 100*nanvar(pcamdl.rotated.comps,[],1)./trace(cov(X));
        pcamdl.rotated.loads = corr(Xorig,pcamdl.rotated.comps,'rows','pairwise');
        if do_boot
            for i = 1:nperms
                resamp = ceil(rand(size(X,1),1)*size(X,1));
                [~,comps_resam] = pca(X(resamp,:),argsin{:});
                comps_resam = comps_resam(:,1:ncompsrot);
                [~,comps_resam] = procrustes(tmpcomps(resamp,:),comps_resam);
                loads_resam(:,:,i) = corr(X(resamp,:),comps_resam);
            end
            
            pcamdl.rotated.loads_boot = loads_resam;
            pcamdl.rotated.loads_bootz = pcamdl.rotated.loads./std(loads_resam,[],3);
            pcamdl.rotated.loads_bootp = ztop(pcamdl.rotated.loads_bootz);
        end
    end
end