function [fieldarr] = getfield_list(listin,field)

if ~iscell(listin)
    listin = vert(listin);
    listin = mat2cell(listin,ones(length(listin),1));
end

fieldarr = {}; marker = zeros(1,length(listin));
for i = 1:length(listin)
    if isfield_nest(listin{i},field)
        try
        tmp = getfield_nest(listin{i},field);
        %if any(size(tmp)~=1)
        
        fieldarr{i} = tmp;
        catch
            marker(i) = 1;
        end
    end
end

for i = 1:length(listin)
    if ~any(isfield_nest(listin{i},field)) || marker(i)==1
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
        elseif iscell(template)
            warning('Some structures don''t have this field - filling with empty cells instead')
            fieldarr{i} = {};
        else
            error('Some structures don''t have this field - can''t fill with NaNs')
        end
    end
    %marker = 0;
end

if ~any(cellfun(@isnumeric,fieldarr,'uniformoutput',true)==0) && ~any(cellfun(@length,fieldarr,'uniformoutput',true)~=1)
    fieldarr = [fieldarr{:}];
end