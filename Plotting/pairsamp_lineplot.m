function pairsamp_lineplot(datin,cols,linecols)

if nargin < 2
   cols = palecol(lines);
end

if nargin < 3
   linecols = palecol([1 0 0;0 0 1;0 0 0]); 
end

hold on
for i = 1:size(datin,2)
    scatter(ones(size(datin,1),1)*i,datin(:,i),96,cols(i,:),'d','filled')
end

for i = 1:(size(datin,2)-1)
    slope = datin(:,i+1)-datin(:,i);
    slopepos = slope>0;
    plot([ones(size(datin(slopepos,i),1),1)*i ones(size(datin(slopepos,i),1),1)*(i+1)]',[datin(slopepos,i) datin(slopepos,i+1)]',...
        'LineWidth',2,'Color',linecols(1,:))
    plot([ones(size(datin(~slopepos,i),1),1)*i ones(size(datin(~slopepos,i),1),1)*(i+1)]',[datin(~slopepos,i) datin(~slopepos,i+1)]',...
        'LineWidth',2,'Color',linecols(2,:))
    slopezero = slope==0;
    plot([ones(size(datin(slopezero,i),1),1)*i ones(size(datin(slopezero,i),1),1)*(i+1)]',[datin(slopezero,i) datin(slopezero,i+1)]',...
        'LineWidth',2,'Color',linecols(3,:))
end

FixAxes(gca,14)
set(gca,'XLim',[0.75 size(datin,2)+0.25]);