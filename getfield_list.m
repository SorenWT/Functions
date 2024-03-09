function [fieldarr] = getfield_list(listin,field)

if ~iscell(listin)
    listin = vert(listin);
    listin = mat2cell(listin,ones(length(listin),1));
end

fieldarr = {};
for i = 1:length(listin)
    if isfield_nest(listin{i},field)
        tmp = getfield_nest(listin{i},field);
        %if any(size(tmp)~=1)
        fieldarr{i} = tmp;
        %else
        %    fieldarr = [fieldarr tmp];
        %end
    end
end

for i = 1:length(listin)
    if ~isfield_nest(listin{i},field)
        template = fieldarr{find(cellfun(@(d)~isempty(d),fieldarr),1)};
        if istable(template)
            warning('Some structures don''t have this field - filling with NaNs instead')
            tmp = template;
            tmp{:,:} = NaN(size(template));
            fieldarr{i} = tmp;
        elseif isnumeric(template)
            warning('Some structures don''t have this field - filling with NaNs instead')
            fieldarr{i} = NaN(size(template));
        elseif ischar(template)
            warning('Some structures don''t have this field - filling with empty strings instead')
            fieldarr{i} = '';
        else
            error('Some structures don''t have this field - can''t fill with NaNs')
        end
    end
end

if ~any(cellfun(@isnumeric,fieldarr,'uniformoutput',true)==0) && ~any(cellfun(@length,fieldarr,'uniformoutput',true)~=1)
    fieldarr = [fieldarr{:}];
end