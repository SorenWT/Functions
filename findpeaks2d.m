function [pks,linlocs,sublocs] = findpeaks2d(matin,edges)

if nargin < 2
   edges = 0; 
end

ispk = zeros(size(matin));

if edges == 1
    ispk = ones(size(matin));
    ispk(2:end-1,2:end-1) = 0;
end

for i = 1:size(matin,1)
    [~,locs] = findpeaks(matin(i,:));
    ispk(i,locs) = ispk(i,locs)+1;
end

for i = 1:size(matin,2)
    [~,locs] = findpeaks(matin(:,i));
    ispk(locs,i) = ispk(locs,i)+1;
end

pks = matin(find(ispk==2));
linlocs = find(ispk==2);
[sublocs(1,:),sublocs(2,:)] = ind2sub(size(matin),linlocs);

