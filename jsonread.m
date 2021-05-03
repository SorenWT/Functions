function [obj] = jsonread(filename,multiline)

if nargin < 2
   multiline = 1; 
end

fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw');

str = erase(str,'myid'); 

% fix this later
if any(str==newline) && multiline
    allstr = tokenize(str,newline);
    allstr(cellfun(@isempty,allstr,'uniformoutput',true)) = [];
    for i = 1:length(allstr)
        if allstr{i}(1) ~= '[' && allstr{i}(1) ~= '{'
        obj{i} = jsondecode(['["' allstr{i} '"]']);
        else
            obj{i} = jsondecode(allstr{i});
        end
    end
else
    obj = jsondecode(str);
end