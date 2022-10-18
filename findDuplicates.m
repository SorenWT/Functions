function [dups,arr] = findDuplicates(arr,firstlast)

unarr = unique(arr);
dups = [];

for i = 1:length(unarr)
    if iscell(arr)
        if sum(strcmpi(unarr{i},arr))>1
            switch firstlast
                case 'first'
                    dups(end+1) = find(strcmpi(unarr{i},arr),2:sum(strcmpi(unarr{i},arr)));
                case 'last'
                    dups(end+1) = find(strcmpi(unarr{i},arr),1:sum(strcmpi(unarr{i},arr))-1);
            end
        end
    else
        if sum(arr==unarr(i))>1
            switch firstlast
                case 'first'
                    dups(end+1) = find(unarr(i)==arr,2:sum(unarr(i)==arr));
                case 'last'
                    dups(end+1) = find(unarr(i)==arr,1:sum(unarr(i)==arr)-1);
            end
        end
    end
end

if any(size(arr)==1)
arr(dups) = [];
end

arr(dups,:) = [];