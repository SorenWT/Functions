function ctrs = centerdesign(des)

if iscell(des)
   des = factor2num(des); 
end

tr = [find(diff(des)) length(des)];

ctrs = tr-diff([0 tr])/2;