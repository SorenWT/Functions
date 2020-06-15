function patchhighlight(xlim,alpha)

if nargin < 2
   alpha = 0.15; 
end

yl = get(gca,'YLim');
patch([xlim(1) xlim(2) xlim(2) xlim(1)],[yl(1) yl(1) yl(2) yl(2)],[0 0 0],'EdgeColor','none','FaceAlpha',alpha);