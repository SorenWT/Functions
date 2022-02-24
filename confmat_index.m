function output = confmat_index(confmat,grp,indxmode)

if nargin < 3
   indxmode = 'selfself'; 
end

if iscell(grp)
    grp = factor2num(grp);
end

ungrp = unique(grp);

switch indxmode
    case 'selfself'
        for i = 1:length(ungrp)
            output(:,:,i) = confmat(grp==i,grp==i);
        end        
end