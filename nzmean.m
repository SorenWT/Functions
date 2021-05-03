function [m] = nzmean(data,dim)

if nargin < 2
   dim = 1; 
end

if isvector(data)
   data = vert(data); 
end

dat = data;
dat(isnan(dat)) = Inf; % convert actual NaNs to Infs to preserve actual nans in the data
dat(dat==0) = NaN;
m = nanmean(dat,dim);
m(isinf(m)) = NaN;




% if nargin < 2
%     dim = 1;
% end
% 
% siz = size(data);
% 
% if sum(siz > 1) == 2
%    for c = 1:siz(mod(dim,2)+1)
%       if dim == 1
%         tmp = data(:,c);
%       else
%           tmp = data(c,:);
%       end
%       med(c) = mean(tmp(find(tmp ~= 0))); 
%       if isnan(med(c))
%           med(c) = 0;
%       end
%    end
% else
%     for c = 1:siz(end)
%        med{c} = nzmean(GetDim(data,siz(end)),dim);
%     end
% end
