function cistr = ci2string(ci)

ci = niceround(ci);
cistr = ['[' num2str(ci(1)) ',' num2str(ci(2)) ']'];