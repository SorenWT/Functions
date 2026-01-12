function heatmap_swt(mat,xlab,ylab)

imagesc(mat)
set(gca,'XTick',1:size(mat,2),'YTick',1:size(mat,1),...
    'XTickLabel',xlab,'YTickLabel',ylab)
xtickangle(45)