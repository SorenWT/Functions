function [X,grpmeanstbl] = nancenter(X,dim)

if nargin < 2
    dim = 1;
end

if istable(X)
    numerics = varfun(@isnumeric,X); numerics = numerics{:,:};
    
    
    if ischar(dim)
        grpmeans = grpstats(X{:,numerics},X.(dim));
        
        [~,~,idc] = unique(X.(dim)) ;
        counts = accumarray(idc,ones(size(idc))) ;
        
        grpmeans = repelem(grpmeans,counts,1);
        grpmeanstbl = array2table(grpmeans,'VariableNames',strcat(X.Properties.VariableNames(numerics),['_' dim 'mean']));
 
        X{:,numerics} = X{:,numerics}-grpmeans;

    else 
    X{:,numerics} = X{:,numerics}-nanmean(X{:,numerics},dim);
    end
else
X = X-nanmean(X,dim);
end