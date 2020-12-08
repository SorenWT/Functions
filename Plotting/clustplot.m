function cbar = clustplot(mat,clusts)

if length(clusts) ~= size(mat,1) 
    clusts = kmeans(mat,clusts,'distance','correlation');
end

[sortclusts,sortindx] = sort(clusts);

imagesc(mat(sortindx,sortindx));

hold on
l = diff(sortclusts);
lindx = find(l > 0)+0.5;

xl = [0.5 length(clusts)+0.5]; yl = xl;

for i = 1:length(lindx)
    line([lindx(i) lindx(i)],yl,'LineWidth',2,'color','k')
    line(xl,[lindx(i) lindx(i)],'LineWidth',2,'color','k')
end
set(gca,'XLim',xl,'YLim',yl,'YDir','reverse');
set(gca,'CLim',[min(belowDiag(mat)) max(belowDiag(mat))])
cbar = colorbar;