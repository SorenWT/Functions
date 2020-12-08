function neighbplot(datamat,distance)

if nargin < 2
    distance = 'correlation';
end

[~,dists] = knnsearch(datamat,datamat,'distance',distance);

distsort = sort(dists);

plot(distsort)