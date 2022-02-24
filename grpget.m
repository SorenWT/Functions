function grpdat = grpget(varin,grp,indx)

if nargin < 3
   indx = 1; 
end

grp = cellstr(categorical(grp));

ungrp = unique(grp);

for i = 1:length(ungrp)
    grpdat(i) = indexme(varin(strcmpi(grp,ungrp{i})),indx);
end

if iscell(grpdat) && all(cellfun(@isnumeric,grpdat))
   grpdat = [grpdat{:}]; 
end

grpdat = vert(grpdat);