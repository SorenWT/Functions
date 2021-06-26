function argsout = removeargs(argsin,argsrm)

for i = 1:length(argsrm)
    argsin(find(strcmpi(argsin,argsrm{i})):find(strcmpi(argsin,argsrm{i}))+1) = [];
end

argsout = argsin;