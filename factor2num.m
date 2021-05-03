function numout = factor2num(factorin)

u = unique(factorin);

for i = 1:length(factorin)
    numout(i) = find(strcmpi(factorin{i},u));
end