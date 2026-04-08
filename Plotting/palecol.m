function c = palecol(c,palefact)

if nargin < 2
   palefact = 0.5; 
end

if length(palefact)>1
   for i = 1:length(palefact)
       cout(i,:) = palecol(c,palefact(i));
   end
   c = cout;
else

t = [1 1 1];
d = t - c;
c = c + (d * palefact);

end