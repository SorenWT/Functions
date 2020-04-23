function [r,p] = partcorr_pairwise(x,y,z,varargin)

if ~CheckInput(varargin,'type')
    type = 'pearson';
else
    type = EasyParse(varargin,'type');
end

for i = 1:size(x,2)
    [r(i,:),p(i,:)] = partialcorr(x(:,i),y,z(:,i),'type',type);
end