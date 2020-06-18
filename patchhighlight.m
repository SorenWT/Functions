function patchhighlight(xlim,alpha,clr)

if nargin < 2
   alpha = 0.15; 
end

if nargin < 3
   clr = [0 0 0]; 
end

yl = get(gca,'YLim');
patch([xlim(1) xlim(2) xlim(2) xlim(1)],[yl(1) yl(1) yl(2) yl(2)],clr,'EdgeColor','none','FaceAlpha',alpha);