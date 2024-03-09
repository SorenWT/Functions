function [h,p]=nanttest(data)

if all(isnan(data))
   h = NaN; p = NaN;
else
    [h,p] = ttest(data);
end