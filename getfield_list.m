function [fieldarr] = getfield_list(listin,field)

fieldarr = {};
for i = 1:length(listin)
   if isfield_nest(listin{i},field)
      tmp = getfield_nest(listin{i},field);
      if any(size(tmp)~=1)
            fieldarr = [fieldarr {tmp}];
      else
          fieldarr = [fieldarr tmp];
      end
   end
end

if ~any(cellfun(@isnumeric,fieldarr,'uniformoutput',true)==0) && ~any(cellfun(@length,fieldarr,'uniformoutput',true)~=1)
   fieldarr = [fieldarr{:}]; 
end