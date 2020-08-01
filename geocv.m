function cvout = geocv(dat,dim)

if nargin < 2
    dim = find(size(dat) > 1,1); %acts along the first non-singleton dimension
end

logdat = log(dat);
cvout = sqrt(exp(nanstd(logdat,[],dim).^2)-1);