function outliers = isoutlier_2d(X)

if any(all(isnan(X),2))
   Xorig = X; 
   nanindx = all(isnan(X),2);
    X(nanindx,:) = [];
else
    nanindx = zeros(size(X,1),1);
end


allouts = isoutlier(X,1);

allout_sum = sum(allouts,2);

outliers = isoutlier(allout_sum);

alloutliers = zeros(size(Xorig,1),1);
alloutliers(~nanindx) = outliers;
outliers = logical(alloutliers);


