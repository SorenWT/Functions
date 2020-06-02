function matout = bdiagtomat(vectin,matsize)
% vectin is the output of belowDiag
% matsize is a scalar, the size of the original matrix


matout = zeros(matsize,matsize);
vectindx = 0;
count = 1;
for i = matsize:-1:2
    matindx = 1:i-1;
    vectindx = matindx+max(vectindx);
    matout(matindx+count,i) = vectin(vectindx);
    count = count+1;
end
matout = fliplr(matout);