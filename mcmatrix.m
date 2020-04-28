function [mcmat] = mcmatrix(siz)

mcmat = cell(siz);

for i = 1:siz(1)
    for ii = 1:siz(2)
        mcmat{i,ii} = [i ii];
    end
end