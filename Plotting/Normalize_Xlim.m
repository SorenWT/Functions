function [newlim] = Normalize_Xlim(fighandle,equal)

if nargin < 2
    equal = 0;
end

if isa(fighandle,'matlab.ui.Figure')
ax = findall(fighandle,'Type','Axes');
elseif isa(fighandle,'matlab.graphics.axis.Axes')
    ax = fighandle;
end

for c = 1:length(ax)
    try
        Xlims(c,:) = get(ax(c),'XLim');
    catch
    end
end

newlim = [min(Xlims(:,1)) max(Xlims(:,2))];

if equal
    maxabs = max(abs(newlim));
    newlim = [-maxabs maxabs];
end
for c = 1:length(ax)
    try
        set(ax,'XLim',newlim)
    catch
    end
end