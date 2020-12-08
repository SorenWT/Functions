function r = regresid(x,y)
% convenience function for getting regression residuals

if any(size(x)==1)
   x = vert(x); 
end

y = vert(y);

[~,~,r] = regress(y,x);