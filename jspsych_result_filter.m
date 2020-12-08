function trlsout = jspsych_result_filter(listin,varargin)

argsin = varargin;
if length(varargin)>2
    argsin = reshape(argsin,[],2);
end
fields = argsin(:,1); values = argsin(:,2);


trlsout = {};
for i = 1:length(listin)
    allcmp = [];
    for q = 1:length(fields)
        allcmp(q) = checktrl(listin{i},fields{q},values{q});
    end
    cmp = ~any(allcmp==0);
    if cmp
        trlsout = [trlsout listin{i}];
    end
end

end

function [cmp] = checktrl(trl,field,value)
if isfield(trl,field)
    if ischar(value)
        cmp = strcmpi(trl.(field),value);
    else
        cmp = trl.field==value;
    end
else
    cmp = 0;
end
end
