function tb = easyannotate(text,pos,varargin)

argsin = varargin;

argsin = setdefault(argsin,'LineStyle','none'); argsin = setdefault(argsin,'FontSize',14);

ax = gca;
edges = [ax.Position(1) ax.Position(2) ax.Position(1)+ax.Position(3) ax.Position(2)+ax.Position(4)];

tb = annotation('textbox','String',text,'FitBoxToText','on',argsin{:});
tbsize = get(tb,'Position');
delete(tb)
if nargin < 2 || isempty(pos)
    pos = ax.Position;
    pos = [edges(3)-tbsize(3)-pos(3)*0.05 edges(4)-tbsize(4)-pos(4)*0.05 tbsize(3) tbsize(4)];
end
%edges = [pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)]; %left bottom right top
tb = annotation('textbox','Position',pos,...
    'String',text,'FitBoxToText','on',argsin{:});