function [c,d,rmout] = RemoveOutliers(a,b,method)

if nargin < 3
   method = 'median';
end

if nargin < 2
    b = NaN(size(a));
end

rmout = [isoutlier(a,method) isoutlier(b,method)];
rmout = find(any(rmout,2));
a(rmout) = [];
b(rmout) = [];
c = a;
d = b;