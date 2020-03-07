function [b,r,p,n] = scatterWithLinearFit(x,y,xyLim,color,noPlot)
% scatterWithLinearFit
% Scatter plot with fitted line and correlation coefficient
%
% scatterWithLinearFit(x,y) creates a scatter plot with circles at the
% locations specified by the vectors x and y, and a least-squares fitted 
% line drawn in addition to the points. A box indicates the number of valid 
% data points, removing NaNs, the correlation coefficient and the
% corresponding p-value. 
%
% scatterWithLinearFit(x,y,xyLim) allows to set the x and y axis limits.
% Specify xyLim as a four-element vector of the form [xmin xmax ymin ymax].
%
% scatterWithLinearFit(x,y,xyLim,color) allows to set the color for the
% data points and fitted line.
%
% [b,r,p,n] = scatterWithLinearFit(__) returns the regression coefficients 
% (intercept and slope), correlation coefficient and corresponding p-value, 
% and the number of valid data pairs (excluding NaNs). Set noPlot to TRUE 
% to obtain the values without creating the plot.

% Anouk de Brouwer

cla reset % clear and reset figure axes
delete(findall(gcf,'type','annotation')) % delete any previous annotations

% make sure variables are in columns and remove NaNs
x = x(:);
y = y(:);
notNan = ~isnan(x) & ~isnan(y);
x = x(notNan);
y = y(notNan);
n = length(x);

% specify axis limits if not provided
if nargin==2 || isempty(xyLim)
    xyLim = [min(x)-0.05*range(x) max(x)+0.05*range(x),...
        min(y)-0.05*range(y) max(y)+0.05*range(y)];
end

% specify color if not provided
if nargin<4 || isempty(color)
    colors = get(gca,'colororder');
    color = colors(1,:);
end

% set noPlot to false if not set
if nargin<5
    noPlot = false;
end

% compute linear fit and correlation
b = polyfit(x,y,1);
[r,p] = corr(x,y);

% plot scatter
if ~noPlot
    plot(x,y,'o','color',color) 
    hold on
    
    % add linear fit
    xfit = xyLim(1:2);
    yfit = polyval(b,xfit);
    plot(xfit,yfit,'color',color)
    axis(xyLim); axis square
    plot(xyLim(1)-100,xyLim(2)-100,'ko') % add point outside of plot for legend hack
    
    % add correlation values
    str = {['n = ' num2str(n)],['r = ' num2str(r,'%1.2f')],['p = ' num2str(p,2)]};
    l = legend(str,'Location','best');
    a = annotation('textbox','String',str,'FitBoxToText','on');
    a.Position = [l.Position(1:2) a.Position(3:4)];
    delete(l);
    a.BackgroundColor = 'w';
    a.EdgeColor = color;
    a.Color = color;
end
