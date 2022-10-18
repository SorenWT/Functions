function plt = gscatter3(x,y,z,g,varargin)

ung = unique(g); ung(isnan(ung)) = [];

l = lines;

clrs = NaN(length(x),3);

for i = 1:length(ung)
    clrs(g==ung(i),:) = repmat(l(i,:),sum(g==ung(i)),1);
end


scatter3(x,y,z,36,clrs,'filled')