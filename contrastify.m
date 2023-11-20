function contr = contrastify(contr)

contr(contr~=0) = contr(contr~=0)-nanmean(contr(contr~=0));