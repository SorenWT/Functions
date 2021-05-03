function tf = regexp_contains(str,pattern)

if ~iscell(str)
    str = {str};
end
if ~iscell(pattern)
    pattern = {pattern};
end

tf = zeros(size(str));
for i = 1:length(str)
    for ii = 1:length(pattern)
        if ~isempty(regexp(str{i},pattern{ii}))
            tf(i) = 1;
            break;
        end
    end
end