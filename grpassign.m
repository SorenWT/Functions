function datout = grpassign(datin,grp)

grp = cellstr(categorical(grp));

ungrp = unique(grp);

if iscell(datin)
   datout = cell(length(grp),1);
else
    datout = NaN(length(grp),1);
end

for i = 1:length(ungrp)
    datout(strcmpi(grp,ungrp{i})) = datin(i);
end