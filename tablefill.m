function tableout = tablefill(tbl,template)

if height(tbl)>1
   for i = 1:height(tbl)
      tableout(i,:) = tablefill(tbl(i,:),template); 
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