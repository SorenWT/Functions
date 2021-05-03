function [fieldarr] = getfield_list(listin,field)

fieldarr = {};
for i = 1:length(listin)
    if isfield_nest(listin{i},field)
        tmp = getfield_nest(listin{i},field);
        %if any(size(tmp)~=1)
        fieldarr = [fieldarr {tmp}];
        %else
        %    fieldarr = [fieldarr tmp];
        %end
    else
        if istable(fieldarr{i-1})
            warning('Some structures don''t have this field - filling with NaNs instead')
            tmp = fieldarr{i-1};
            tmp{:,:} = NaN(size(tmp));
            fieldarr = [fieldarr {tmp}];
        elseif isnumeric(fieldarr{i-1})
            warning('Some structures don''t have this field - filling with NaNs instead')
            fieldarr = [fieldarr {NaN(size(fieldarr{i-1}))}];
        else
            error('Some structures don''t have this field - can''t fill with NaNs')
        end
    end
end

if ~any(cellfun(@isnumeric,fieldarr,'uniformoutput',true)==0) && ~any(cellfun(@length,fieldarr,'uniformoutput',true)~=1)
    fieldarr = [fieldarr{:}];
end