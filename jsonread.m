function [obj] = jsonread(filename,multiline)

if nargin < 2
   multiline = 1; 
end

fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw');

str = erase(str,'myid'); 

%try
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
%catch
%   obj = jsonread(filename,~multiline); % if it doesn't work, try it with the other value of multiline
%end

fclose(fid);