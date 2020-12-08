function [obj] = jsonread(filename)

fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw');

if any(str==newline)
    allstr = tokenize(str,newline);
    allstr(cellfun(@isempty,allstr,'uniformoutput',true)) = [];
    for i = 1:length(allstr)
        obj{i} = jsondecode(allstr{i});
    end
else
    obj = jsondecode(str);
end