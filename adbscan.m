function allclusts = adbscan(X,k,minpts,res,reduceres)

if nargin < 4
    res = 100;
end

if nargin < 5
    reduceres = 1;
end

% assumes distance is precomputed for now
Xorig = X;
eps = min(min(min(X)),0); % start eps with a very low value
maxeps = max(max(X));

epsres  = (maxeps-eps)/res;


clustindx = 1:size(X,1);

allclusts = zeros(size(X,1),1)-1;
numcl = 0;

while numcl < k
    if eps >= maxeps
        if reduceres
            while eps >= maxeps
                warning('Failed to generate the requested number of clusters with this resolution. Lowering epsilon step and trying again...')
                res = res*0.8;
                allclusts = adbscan(Xorig,k,minpts,res);
                if (length(unique(allclusts))-1) >= k
                    return
                end
            end
        else
            allclusts = NaN(size(allclusts));
            return
        end
    end
    
    clusts = dbscan(X,eps,minpts,'distance','precomputed');
    if sum(clusts>0) > 0.1*length(clusts)
        allclusts(clustindx(clusts>0)) = clusts(clusts>0)+max(max(allclusts),0);
        numcl = length(unique(allclusts))-1;
        
        X(clusts>0,:) = [];
        X(:,clusts>0) = [];
        if isempty(X)
            eps = maxeps+1; % this is also a resolution problem - handle the same way
        end
        clustindx(clusts>0) = [];
    else
        eps = eps+epsres;
    end
end