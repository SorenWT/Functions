function [c,d,rmout] = RemoveOutliers(a,b)

if nargin < 2
    b = NaN(size(a));
end

rmout = [find(isoutlier(a))' find(isoutlier(b))'];
a(rmout) = [];
b(rmout) = [];
c = a;
d = b;