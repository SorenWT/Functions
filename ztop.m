function p = ztop(z)

posmask = z>0;
negmask = z<0;
p = ones(size(z));

p(posmask) = 2*(1-normcdf(z(posmask)));
p(negmask) = 2*normcdf(z(negmask));
% if z>0
%     p = 2*(1-normcdf(z));
% else
%     p = 2*normcdf(z);
% end