function dat = rmnan(dat,dim)

if dim == 1
nandat = sum(dat,2);
elseif dim == 2
   nandat = sum(dat,1); 
end

rmindx = find(isnan(nandat));
if dim == 1
    dat(rmindx,:) = [];
elseif dim == 2
    dat(:,rmindx) = [];
end