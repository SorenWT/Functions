function result = cstat(T,statfield,statfun,condvalue,subj)
% assumes conditions are found in a table variable called overcond

if nargin < 5
   subj = []; 
end

if iscell(condvalue)
    for i = 1:length(condvalue)
        result(i,:) = cstat(T,statfield,statfun,condvalue{i},subj);
    end
else
    
    siz_cond = length(condvalue);
    
    for i = 1:siz_cond
        logvect(:,i) = strcmpi(cellfun(@(d)indexme(d,i),T.overcond,'UniformOutput',false),condvalue(i)) | strcmpi(condvalue(i),'*');
    end
    
    logvect = all(logvect,2);
    
    T = T(logvect,:);

    
    if nargin > 4 && ~isempty(subj)
        % if a grouping variable is supplied, do the requested function on
        % the group means
        result = statfun(grpstats(T.(statfield),T.(subj))); 
    else
        result = statfun(T.(statfield));
    end
    
end