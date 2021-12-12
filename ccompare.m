function result = ccompare(T,statfield,statfun,conds,subj,catdim)

if nargin < 6
    % concatenate along 2nd dimension by default in case different conditions have different numbers of trials
    catdim = 2;
end

for i = 1:length(conds)
    cdat{i} = cstat(T,statfield,@(d)d.*1,conds{i},subj,catdim);
end

result = statfun(cdat);