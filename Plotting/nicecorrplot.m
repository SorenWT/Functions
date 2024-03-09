function [corrrho,corrp,h,f,tb] = nicecorrplot(a,b,labels,varargin)
% nicecorrplot plots a scatter plot and a best-fit line for two variables
%
% nicecorrplot(a,b) plots the correlation between the variables a and b
%
% nicecorrplot(a,b,labels) takes a cell array "labels" and adds its
% elements as axis labels
%
% nicecorrplot(...,'Type',string) takes the input argument to corr
%    ('Pearson','Spearman', or 'Kendall')
% nicecorrplot(...,'Cov',x) takes another vector x and runs a partial
%    correlation between a and b, controlling for x

if EasyParse(varargin,'RemoveOutliers','on')
    rmoutliers = [find(isoutlier(a))' find(isoutlier(b))'];
    
    a(rmoutliers) = [];
    b(rmoutliers) = [];
end

l = lines;
argsin = varargin;
argsin = setdefault(argsin,'scatterclr',palecol(l(1,:),0.33));
scatterclr = EasyParse(argsin,'scatterclr');

% dark red line to show up better vs. dots
argsin = setdefault(argsin,'lineclr',[0.5 0 0]);
lineclr = EasyParse(argsin,'lineclr');

if size(a,1) == 1
    a = a';
end

if size(b,1) == 1
    b = b';
end

if nargin < 3
    labels = {'A','B'};
end

if mean(abs(a)) < 1e-6
    sclfact = mean(abs(a));
    a = a./sclfact;
else
    sclfact = 1;
end

anan = isnan(a); bnan = isnan(b);
allnan = (anan+bnan)>0;
a(allnan) = []; b(allnan) = [];

ainf = isinf(a); binf = isinf(b);
allinf = (ainf+binf)>0;
a(allinf) = []; b(allinf) = [];


B = regress(b,horzcat(ones(length(a),1),a));
B(2) = B(2)./sclfact;
a = a.*sclfact;

if CheckInput(varargin,'Type')
    type = EasyParse(varargin,'Type');
else
    type = 'Spearman';
end

if CheckInput(varargin,'Cov')
    [corrrho,corrp] = partialcorr(a,b,indexme(EasyParse(varargin,'Cov'),~allnan,1:size(EasyParse(varargin,'Cov'),2)),'Type',type,'rows','pairwise');
else
    [corrrho,corrp] = corr(a,b,'rows','pairwise','Type',type);
end

if CheckInput(varargin,'externalp')
    corrp = EasyParse(varargin,'externalp');
end

l = lines;

h = scatter(a,b,384./round(log(length(a))),scatterclr,'filled');
xl = xlabel(labels{1},'FontSize',16);
yl = ylabel(labels{2},'FontSize',16);
hold on;
f = plot(linspace(min(a),max(a),1000),B(1)+B(2)*linspace(min(a),max(a),1000),'color',lineclr);
ax = gca;
pos = ax.Position;

%ylim([-0.1 0.1])

set(f,'LineWidth',1.5)

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'fontsize',14,...
    'LineWidth'   , 1,...
    'XLim',[min(a)-0.05*range(a) max(a)+0.05*range(a)],...
    'YLim',[min(b)-0.05*range(b) max(b)+0.05*range(b)]);

set(xl,'FontSize',20)
set(yl,'FontSize',20)

if strcmpi(type,'spearman')
    rstring = '\rho';
else
    rstring = 'r';
end

if CheckInput(varargin,'Plot') && EasyParse(varargin,'Plot','r')
    tb = annotation('textbox','String',{[rstring ' = ' num2str(round(corrrho,2))];['p = ' num2str(round(corrp,2,'significant'))]},...
        'FitBoxToText','on','LineStyle','none','FontSize',14);
    tbsize = get(tb,'Position');
    delete(tb)
    edges = [pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)]; %left bottom right top
    tb = annotation('textbox','Position',[edges(3)-tbsize(3)-pos(3)*0.05 edges(4)-tbsize(4)-0.05*pos(4) tbsize(3) tbsize(4)],...
        'String',{[rstring ' = ' num2str(round(corrrho,2))]},'FitBoxToText','on','LineStyle','none','FontSize',14);
elseif CheckInput(varargin,'Plot') && EasyParse(varargin,'Plot','Beta')
    tb = annotation('textbox','String',{['\beta = ' num2str(round(B(2),3))]},'FitBoxToText','on','LineStyle','none','FontSize',14);
    tbsize = get(tb,'Position');
    delete(tb)
    edges = [pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)]; %left bottom right top
    tb = annotation('textbox','Position',[edges(3)-tbsize(3)-pos(3)*0.05 edges(4)-tbsize(4)-pos(4)*0.05 tbsize(3) tbsize(4)],...
        'String',{['\beta = ' num2str(round(B(2),2))]},'FitBoxToText','on','LineStyle','none','FontSize',14);
elseif CheckInput(varargin,'Plot') && EasyParse(varargin,'Plot','off')
elseif ~CheckInput(varargin,'Plot')
    tb = annotation('textbox','String',{[rstring ' = ' num2str(round(corrrho,2))];['p = ' num2str(round(corrp,2,'significant'))]},...
        'FitBoxToText','on','LineStyle','none','FontSize',14);
    tbsize = get(tb,'Position');
    delete(tb)
    edges = [pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)]; %left bottom right top
    tb = annotation('textbox','Position',[edges(3)-tbsize(3)-pos(3)*0.05 edges(4)-tbsize(4)-0.05*pos(4) tbsize(3) tbsize(4)],...
        'String',{[rstring ' = ' num2str(round(corrrho,2))];['p = ' num2str(round(corrp,2,'significant'))]},'FitBoxToText','on','LineStyle','none','FontSize',14);
end

if isnan(corrrho) && isnan(corrp)
    delete(tb)
end

FixAxes(gca,14)
