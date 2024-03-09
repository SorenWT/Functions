function [tf,tf2] = regexp_contains(str,pattern)

if ~iscell(str)
    str = {str};
end
if ~iscell(pattern)
    pattern = {pattern};
end

tf = zeros(size(str));
tf2 = zeros(size(pattern));
for i = 1:length(str)
    for ii = 1:length(pattern)
        if ~isempty(regexp(str{i},pattern{ii}))
            tf(i) = 1;
            tf2(ii) = 1;
            %break;
        end
    end
end