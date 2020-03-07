function p = plotMeansWithDataPoints(Y,mY,varY,color,labels)
% plotMeansWithDataPoints
% Boxplot-like plot with mean or median values, variance, and individual
% data points
%
% plotMeansWithDataPoints(Y) creates a boxplot-like plot with points
% displaying the individual values in Y, a horizontal line displaying the
% mean of the values in Y, and a shaded area displaying the standard error
% of the mean of the values in Y. If Y is a matrix, there is one 'box' per
% column; if Y is a vector, there is just one 'box'.
%
% plotMeansWithDataPoints(Y,mY) allows the use of an alternative central
% dispersion measure (e.g., the median) provided in mY.
%
% plotMeansWithDataPoints(Y,mY,varY) allows the use of an alternative
% variance measure (e.g., the interquartile range) provided in mY.
%
% plotMeansWithDataPoints(Y,mY,varY,color) allows to specify one or
% multiple RGB colors. Colors can be used to group variables.
%
% plotMeansWithDataPoints(Y,mY,varY,color,labels) allows set x axis labels,
% or, when multiple colors are provided, to specify the labels for a legend,
% with each cell corresponding to a specified color.
%
% p = plotMeansWithDataPoints(__) returns the figure handle.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

cla reset % clear and reset axes

% compute means if measure of central tendency is not provided
if ~isempty(Y) && (nargin==1 || isempty(mY))
    mY = mean(Y);
end

% compute standard error if measure of variance is not provided
if ~isempty(Y) && (nargin<3 || isempty(varY))
    varY = std(Y)./sqrt(size(Y,1));
end

% use grey if color is not provided
if nargin<4
    color = [0.2 0.2 0.2];
end
color = repmat(color,size(mY,2)/size(color,1),1);
fadedColors = color+(1-color)*0.75;

% create empty cell if labels are not provided
if nargin<5
    labels = {};
end

% plot
nRow = size(Y,1);
nCol = size(mY,2);
barWidth = 0.6;
x = [(1:nCol)-0.5*barWidth; (1:nCol)+0.5*barWidth]';  % x position of bars
xRandOffset = (rand(nRow,1)-0.5)*0.6*barWidth;        % x position of data points
for i = 1 : nCol
    % shading
    if ~isnan(varY(i))
        r = rectangle('Position',[x(i,1),mY(i)-varY(i),barWidth,2*varY(i)],...
            'EdgeColor','none','FaceColor',fadedColors(i,:));
    end
    
    % data points
    if ~isempty(Y)
        hold on
        plot(i+xRandOffset,Y(:,i),'.','markersize',10,'color',color(i,:))
    end
    
    % central tendency
    hold on
    p(i) = plot(x(i,:),[mY(i) mY(i)],'color',color(i,:),'linewidth',2);
    
    % set x axis labels or legend
    if i==length(labels)
        set(gca,'XTick',1:nCol)
        if length(labels)==nCol
            set(gca,'XTicklabels',labels)
        else
            legend(p,labels,'location','best');
        end
    end
end
xlim([1-barWidth nCol+barWidth])

end