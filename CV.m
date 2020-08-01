function cvout = CV(data,dim)

if nargin < 2
    dim = find(size(data) > 1,1); %acts along the first non-singleton dimension
end

if any(data < 0)
    %warning('Using geometric CV with base 10 exponentiation')
    cvout = geocv(10.^data,dim);
else
    cvout = nanstd(data,[],dim)./nanmean(data,dim);
end
