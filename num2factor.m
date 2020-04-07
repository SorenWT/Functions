function factout = num2factor(facts,key,nums)

factout = cell(size(nums));

for i = 1:length(nums)
    factout{i} = facts{find(key==nums(i))};
end