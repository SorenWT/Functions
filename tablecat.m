function outtable = tablecat(tables)
% tables should be a cell array of tables you want to concatenate

for i = 1:length(tables)
    varnames{i} = tables{i}.Properties.VariableNames;
end

allvars = cat(2,varnames{:});
unvars = unique(allvars);

template = table;

for i = 1:length(unvars)
    for q = 1:length(tables)
       if ismember(unvars{i},tables{q}.Properties.VariableNames)
           tmp = tables{q}.(unvars{i});
            if isnumeric(tmp)
                  val = NaN;
            elseif ischar(tmp) || iscell(tmp)
                val = {''};
            elseif isdatetime(tmp)
                val = NaT;
            end
       end
    end
    template.(unvars{i}) = val;
end

for i = 1:length(tables)
    tables{i} = tablefill(tables{i},template);
    [~,match] = match_str(unvars,tables{i}.Properties.VariableNames);
    tables{i} = tables{i}(:,match);
end

outtable = cat(1,tables{:});