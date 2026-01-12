function [r,pred] = regresid(y,x,varargin)
% convenience function for getting regression residuals

if any(size(x)==1)
   x = vert(x); 
end

if size(y,2)==1
y = vert(y);
end

for i = 1:size(y,2)
    tmp = fitlm(x,y(:,i),varargin{:});
    r(:,i) = tmp.Residuals.Raw;
    pred(:,i) = tmp.predict;
    %[~,~,r(:,i)] = regress(y(:,i),[ones(size(x,1),1) x]);
end