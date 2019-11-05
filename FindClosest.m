function [indicesout,itemsout] = FindClosest(list,items,first)
% finds the closest value in list for every value of items
% if there are multiple of these values for a given item, the value of 
% first determines which one to return


if nargin < 3
    first = 1; 
end

for c = 1:length(items)
    %if first==1
        [~,tmp] = min(abs(list-items(c)));
        %tmp = find(abs(list-items(c)) == min(abs(list-items(c))),1);
    %else
        %tmp = find(abs(list-items(c)) == min(abs(list-items(c))));
    %end
    indicesout(c) = tmp(first);
end

itemsout = list(indicesout);