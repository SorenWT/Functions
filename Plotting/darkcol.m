function c = darkcol(c,palefact)

if nargin < 2
   palefact = 0.5; 
end

t = [0 0 0];
d = t - c;
c = c + (d * palefact);