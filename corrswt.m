function [r,p,ci,n] = corrswt(x,y,varargin)

argsin = varargin;

argsin = setdefault(argsin,'rows','pairwise');
argsin = setdefault(argsin,'type','pearson');

x = x.*nanmask(~isinf(x)); y = y.*nanmask(~isinf(y));

%if any(any(isoutlier(x,'grubbs'))) || any(any(isoutlier(y,'grubbs')))
    argsin = setdefault(argsin,'removeoutliers',false);
%end

if EasyParse(argsin,'removeoutliers',true)
    x = x.*nanmask(~isoutlier(x)); y = y.*nanmask(~isoutlier(y));
    %[x,y,outls] = RemoveOutliers(x,y,'grubbs');
    
    warning(['Outliers removed. Check sample sizes for each correlation']);
else
    %if any(any(isoutlier(x,'grubbs'))) | any(any(isoutlier(y,'grubbs')))
        warning(['Outliers detected but not removed. Plot correlations to ensure outliers not driving']);
    %end
end

argsin = removeargs(argsin,{'removeoutliers'});

[r,p] = corr(x,y,argsin{:});

for i = 1:size(x,2)
    for ii = 1:size(y,2)
        n(i,ii) = sum(~any(isnan([x(:,i) y(:,ii)]),2));
        ci(i,ii,:) = corrci(r(i,ii),n(i,ii));
    end
end