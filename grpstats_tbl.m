function Xnew = grpstats_tbl(X,grpvar,grpfun)

if nargin < 3
    grpfun = @nanmean;
end


for i = 1:width(X)
    numerics(i) = isnumeric(X{:,i});
end

Xnew = X;

for i = 1:width(Xnew)
    if isnumeric(X{:,i})
        tmp = grpstats(X{:,i},X.(grpvar),grpfun);
    else
        tmp = grpstats_swt(X{:,i},X.(grpvar));
    end
    
    Xnew = Xnew(1:length(tmp),:);
    if iscell(tmp)
    Xnew(:,i) = vert(tmp);
    else
        Xnew{:,i} = vert(tmp);
    end
    
    
end
end

function grpout = grpstats_swt(dat,grpvar,fun)

if isnumeric(dat)
    grpout = grpstats(dat,grpvar,fun);
else
    ungrp = unique(grpvar,'stable');
    
    for i = 1:length(ungrp)
        grpout{i} = unique(dat(strcmpi(grpvar,ungrp{i})));
    end
end
end

