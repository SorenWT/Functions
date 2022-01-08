function [result,logvect,T] = cstat(T,statfield,statfun,condvalue,subj,catdim)
% cstat calculates summary statistics by condition for the table T
%   it assumes conditions are found in a table variable called overcond
% Inputs:
%   T is the table containing all the behavioural data 
%   statfield is the field that you want to do statistics on (e.g. 'accuracy')
%   statfun is a function handle for the function you want to calculate (e.g. @nanmean)
%   condvalue is the value of overcond that you want to calculate
%   statistics for. You can use the following operators in condvalue: 
%      | (or)
%      ( ) (for order of operations)
%      - (for a difference of conditions. In this case, a subject grouping
%      must be supplied - see below).
%   subj is the table variable containing the subject information (e.g.
%   'subject'). If you input something here, cstat uses the grpstats 
%   function first to get the subject means of statfield, and then applies 
%   the given statfun to the subject means
%   catdim is mostly used internally, don't bother unless you have a
%   specific reason
%
% Outputs: 
%   result: the result of applying statfun to the specified conditions
%   logvect: a vector containing the table indices where the specified
%   condition is found
%   T: the input table, indexed by logvect



if nargin < 5
    subj = [];
end

if nargin < 6
    catdim = 1;
end

if isempty(statfun)
    statfun = @doNothing;
end



if iscell(condvalue)
    for i = 1:length(condvalue)
        result{i} = cstat(T,statfield,statfun,condvalue{i},subj);
        result{i} = horz(result{i});
    end
    result = cat(catdim,result{:});
else
    
    if contains(condvalue,'-') || contains(condvalue,'+')
        if isempty(subj)
            error('Cannot take a sum/difference of conditions if a grouping variable is not supplied!')
        end
        
        [tmpcond,operations] = strsplit(condvalue,{'-','+'});
        for i = 1:length(tmpcond)
            [tmpres{i},~,Tcond{i}] = cstat(T,statfield,@doNothing,tmpcond{i},subj,catdim);
        end

        subs = getfield_list(Tcond,subj);

        common = intersect(subs{:});

        for i = 1:length(tmpcond)
            m = match_str(unique(Tcond{i}.(subj)),common);
            tmpres{i} = tmpres{i}(m); 
        end

        tmpresstring = cellcat('tmpres{',cellstr(num2str([1:length(tmpres)]')),'',0);
        tmpresstring = cellcat('}',tmpresstring,'',1);
        
        dif = eval(strjoin(tmpresstring,operations));
        result = statfun(dif);
    elseif contains(condvalue,'(')
        tmpcond = extractBetween(condvalue,'(',')');
        result = cstat(T,statfield,@doNothing,tmpcond,subj,catdim);
        
    elseif contains(condvalue,'|')
        tmpcond = strsplit(condvalue,'|');
        for i = 1:length(tmpcond)
            [~,logvect(:,i)] = cstat(T,statfield,@doNothing,tmpcond{i},subj,catdim);
        end
        
        logvect = any(logvect,2);

        T = T(logvect,:);
        
        if ~isempty(subj)
            result = statfun(grpstats(T.(statfield),T.(subj)));
        else
            result = statfun(T.(statfield));
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
end