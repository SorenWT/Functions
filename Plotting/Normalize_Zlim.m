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
        Zlims(c,:) = get(ax(c),'ZLim');
    catch
    end
end

newlim = [min(Zlims(:,1)) max(Zlims(:,2))];

if equal
    maxabs = max(abs(newlim));
    newlim = [-maxabs maxabs];
end
for c = 1:length(ax)
    try
        set(ax,'ZLim',newlim)
    catch
    end
end