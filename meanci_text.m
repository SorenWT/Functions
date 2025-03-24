function [textout] = meanci_text(dat, cat, citype)

if nargin < 3
    citype = 'ci95';
end

if isnumeric(dat) && (length(unique(dat(~isnan(dat)))) > 2)
    m = nanmean(dat);
    switch citype
        case 'ci95'
            se = nanstd(dat)./sqrt(sum(~isnan(dat)));
            ci = [m-1.96*se m+1.96*se];
            textout = [num2str(niceround(m)) ', [' num2str(niceround(ci(1))) ',' num2str(niceround(ci(2))) ']'];
        case 'sd'
            ci = nanstd(dat);
            textout = [num2str(niceround(m)) ' (' num2str(niceround(ci)) ')'];
        case 'se'
            ci = nanstd(dat)./sqrt(sum(~isnan(dat)));
            textout = [num2str(niceround(m)) ' (' num2str(niceround(ci)) ')'];
    end
    
    
else
    m = nansum(dat);
    [~,ci] = binofit(m,sum(~isnan(dat)));
    ci = ci*sum(~isnan(dat));
    ci = round(ci);
    
    textout = [num2str(niceround(m)) ' ' cat ', [' num2str(niceround(ci(1))) ',' num2str(niceround(ci(2))) ']'];
    
end



