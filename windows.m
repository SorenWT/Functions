function winsout = windows(datin,winsize,stepsize)

if nargin<3
   stepsize = 1; 
end

steps = 1:stepsize:(size(datin,2)-winsize+1);

for i = 1:length(steps)
    winsout(:,:,i) = datin(:,steps(i):steps(i)+(winsize-1));
end