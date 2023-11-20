function tableout = tablefill(tbl,template)

rnames = tbl.Properties.RowNames;

if height(tbl)>1
    for i = 1:height(tbl)
        tableout(i,:) = tablefill(tbl(i,:),template);
    end
    if ~isempty(tbl.Properties.RowNames)
       tableout.Properties.RowNames = tbl.Properties.RowNames; 
    end
    return
end

tofill = setdiff(template.Properties.VariableNames,tbl.Properties.VariableNames);

for i = 1:length(tofill)
    tmp = template.(tofill{i});
    if isnumeric(tmp)
        val = NaN;
    elseif ischar(tmp) || iscell(tmp)
        val = {''};
    elseif isdatetime(tmp)
        val = NaT;
    end
    tbl.(tofill{i}) = val;
end

tableout = tbl;
if ~isempty(rnames)
    tableout.Properties.RowNames = rnames;
end

[~,reorder] = match_str(template.Properties.VariableNames,tableout.Properties.VariableNames);
tableout = tableout(:,reorder);