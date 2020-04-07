function cvout = CV(data,dim)

if nargin < 2
    dim = find(size(data) > 1,1); %acts along the first non-singleton dimension
end

cvout = std(data,[],dim)./mean(data,dim);