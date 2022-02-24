function Plot_sigmask(ax,sigmask,type,varargin)
% Plot_sigmask plots either an overbar or a shaded area to indicate time
% points which are significant
% Inputs: 
%   ax: the MATLAB axis object to plot into (if the current axis, use gca)
%   sigmask: the significance mask to plot. This should be either a vector
%       with the same length as the time series you plotted, or a matrix of
%       size nchannels x ntimepoints if you want to use the cmapline option
%   type: the type of graphic to plot: 'bar','shade', or 'cmapline'. Bar
%       plots a line above your plot indicating which time points are
%       significant, shade colours the significant time points with a
%       transparent grey patch, and cmapline plots a graded overbar with the
%       number of channels which are significant at each time point
% Optional inputs (name-value pairs):
%   'color',mycolor: change the colour of the plotted object
%   'cmap',mycolormap: if using the cmapline option, change the colour map
%   of the plotted bar
%   'alpha',myalpha: change the alpha level of the plotted object. 
%   'linewidth', mylinewidth: change the width of the plotted line (applies
%   to bar and cmapline options only).

hold(gca,'on')

if nargin < 3
    type = 'bar';
end

if ~isempty(find(size(sigmask) == 1))
    sigmask = horz(sigmask);
end

%sigmask = sigmask > 0; % just to make sure it's zeros and ones

if strcmpi(type,'bar')
    varargin = setdefault(varargin,'color',[0.5 0.5 0.5]);
else
    varargin = setdefault(varargin,'color',[0.15 0.15 0.15]);
end
if strcmpi(type,'bar')
    varargin = setdefault(varargin,'alpha',0.5);
else
    varargin = setdefault(varargin,'alpha',0.10);
end
varargin = setdefault(varargin,'cmap',gray);
varargin = setdefault(varargin,'linewidth',2);
varargin = setdefault(varargin,'yoffset',0.04);

color = EasyParse(varargin,'color');
cmap = EasyParse(varargin,'cmap');
alpha = EasyParse(varargin,'alpha');
linewidth = EasyParse(varargin,'linewidth');
yoffset = EasyParse(varargin,'yoffset');

xlim = get(gca,'XLim');
xax = linspace(xlim(1),xlim(2),size(sigmask,2)); %time points must be the second dimension of sigmask

if strcmpi(type,'shade')
    yl = get(ax,'YLim');
    patchstep = sigmask*(yl(2)-yl(1));
    patchstep = patchstep + yl(1);
    area(xax,patchstep,yl(1),'FaceAlpha',alpha,'LineStyle','none','FaceColor',color,'HandleVisibility','off')
    set(ax,'YLim',yl)
elseif strcmpi(type,'bar')
    yl = get(ax,'YLim');
    ysize = yl(2)-yl(1);
    diffmask = diff([0 sigmask]);
    startindices = find(diffmask == 1);
    stopindices = find(diffmask == -1);
    if length(stopindices) < length(startindices)
        stopindices = [stopindices length(sigmask)];
    end
    indices = [vert(startindices) vert(stopindices)];
    
    if ~isempty(indices)
        for c = 1:size(indices,1)
            line([xax(indices(c,1)) xax(indices(c,2))],[yl(2)-ysize*0.04 yl(2)-ysize*0.04],...
                'color',[color alpha],'linewidth',linewidth);
        end
    end
elseif strcmpi(type,'cmapline') % for plotting the number of significant channels as the colour of the line
    yl = get(ax,'YLim');
    ymin = findMinY(get(ax,'XLim'));
    ysize = yl(2)-yl(1);
    y = repmat([ymin+ysize*yoffset],size(xax));
    %y = repmat([yl(2)-ysize*0.04],size(xax));
    z = zeros(size(xax));
    col = sum(sigmask,find(size(sigmask) ~= length(xax)));
    cindex = linspace(size(sigmask,1),0,size(cmap,1));
    col = cmap(FindClosest(cindex,col,1),:);
    col = permute(col,[3 1 2]);
    
    surface([xax;xax],[y;y],[z;z],[col;col],...
        'facecolor','no','edgecolor','interp','linewidth',linewidth)
    set(gca,'YLim',[yl(1) y(1)+ysize*yoffset])
    %colormap(cmap)
end

end


function Y=findMinY(x)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');

    axis(gca,'tight')
    
    %increase range of x values by 0.1 to ensure correct y max is used
    x(1)=x(1)-0.1;
    x(2)=x(2)+0.1;
    
    set(gca,'xlim',x) %Matlab automatically re-tightens y-axis

    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = max(yLim);

    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',oldYLim)

end %close findMinY