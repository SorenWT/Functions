function [med] = nzmean(data,dim)
% currently only works with up to 2D data 

if nargin < 2
    dim = 1;
end

siz = size(data);

if sum(siz > 1) == 2
   for c = 1:siz(mod(dim,2)+1)
      if dim == 1
        tmp = data(:,c);
      else
          tmp = data(c,:);
      end
      med(c) = mean(tmp(find(tmp ~= 0))); 
      if isnan(med(c))
          med(c) = 0;
      end
   end
else
    for c = 1:siz(end)
       med{c} = nzmean(GetDim(data,siz(end)),dim) 
    end
end
