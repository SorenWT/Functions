function numout = factor2num(factorin,uniqueord)

if nargin < 2
u = unique(factorin);
else
    u = uniqueord;
end

for i = 1:length(factorin)
    numout(i) = find(strcmpi(factorin{i},u));
end