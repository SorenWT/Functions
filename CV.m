function cvout = CV(data,dim)

if nargin < 2
    dim = find(size(data) > 1,1); %acts along the first non-singleton dimension
end

cvout = nanstd(data,[],dim)./nanmean(data,dim);