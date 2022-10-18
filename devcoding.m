function tblout = devcoding(tblin,var)

unvars = unique(tblin.(var));

dvarmat = dummyvar(categorical(tblin.(var)));

tblout = tblin;
for i = 1:length(unvars)
    dvarmat(:,i) = dvarmat(:,i)-(1/size(dvarmat,2));
    tblout.([var '_' unvars{i}]) = dvarmat(:,i);
end