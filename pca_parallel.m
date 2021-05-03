function pcamdl = pca_parallel(X,nperms,varargin)

argsin = varargin;
argsin = setdefault(argsin,'do_boot',true);

do_boot = EasyParse(argsin,'do_boot');

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
    
    argsin = varargin;
    argsin(find(strcmpi(argsin,'do_boot')):find(strcmpi(argsin,'do_boot'))+1) = [];
    
    [weights,comps,~,~,explained] = pca(X,argsin{:});
    
    pcamdl.weights = weights;
    pcamdl.comps = NaN(size(Xorig,1),size(comps,2)); pcamdl.comps(~nanindx,:) = comps;
    
    pcamdl.explained = explained;
    
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
    pcamdl.kaiser = find(explained<100/size(X,2),1)-1;
    if isempty(pcamdl.kaiser)
        pcamdl.kaiser = length(pcamdl.explained);
    end
    
    pcamdl.loads = corr(X,comps);
    
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
    
end