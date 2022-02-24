neighbs = ecc_getneighbs(struct('atlas',mmpgroup,'atlasname','mmp'))

% no neighbours for the ??? regions
neighbs(1).neighblabel = {};
neighbs(24).neighblabel = {};

tmplabs = erase('L_',mmpgroup.parcellationlabel(1:23));

midlineindx = 1+[1 2 3 4 6 7 8 13 16 18 19 20 22];
midlineregions = tmplabs(midlineindx);

for i = 1:length(midlineindx)
    neighbs(midlineindx(i)).neighblabel = [neighbs(midlineindx(i)).neighblabel ...
        {neighbs(midlineindx(i)+23).label}];
end

for i = 1:length(midlineindx)
    neighbs(midlineindx(i)+23).neighblabel = [neighbs(midlineindx(i)+23).neighblabel ...
        {neighbs(midlineindx(i)).label}];
end
