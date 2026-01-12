function factout = num2factor(nums,facts,key)

uniquenums = unique(nums(~isnan(nums)));


if nargin < 2 || isempty(facts)
    % only works up to 676 right now
    alphabet = 'abcdefghijklmnopqrstuvwxyz';
    numletters = ceil(log(length(nums))/log(26));
    for i = 1:length(uniquenums)
        tmp = [];
        if i > 26
            tmp = [tmp alphabet(floor(uniquenums(i)/26))];
        end
        tmp = [tmp alphabet(1+uniquenums(i)-26*floor(uniquenums(i)/26))];
        facts{i} = tmp;
    end
    
end

if nargin < 3
    key = uniquenums;
end

factout = cell(size(nums));

for i = 1:length(nums)
    if ~isnan(nums(i))
        factout{i} = facts{find(key==nums(i))};
    else
        factout{i} = '';
    end
end