function cvout = geocv(dat,dim)

if nargin < 2
    dim = find(size(dat) > 1,1); %acts along the first non-singleton dimension
end

if any(dat < 0)
    % warning('Negative values found: assuming log-transformed data, exponentiating with base 10');
    dat = 10.^dat;
end
logdat = log(dat);
cvout = sqrt(exp(nanstd(logdat,[],dim).^2)-1);