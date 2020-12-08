function z = rtoz(r,infopt)

z = 0.5*log((1+r)./(1-r));

if nargin > 1 && strcmpi(infopt,'zeroinf')
    z(isinf(z)) = 0;
elseif nargin > 1 && strcmpi(infopt,'naninf')
    z(isinf(z)) = NaN;
end