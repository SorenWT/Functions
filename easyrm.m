function [rm,widetbl] = easyrm(tbl,modelspec,timevar,untime,withinmodel)

if nargin < 5
   withinmodel = 'separatemeans'; 
end

if ~exist('untime','var')
untime = unique(tbl.(timevar));
end
widetbl = table;

for i = 1:length(untime)
    if isempty(widetbl)
        widetbl = tbl(strcmpi(tbl.(timevar),untime{i}),:);
        widetbl.Properties.VariableNames = strcat(replace(untime{i},'-','_'),'_',widetbl.Properties.VariableNames);
    else
        tmptbl = tbl(strcmpi(tbl.(timevar),untime{i}),:);
        tmptbl.Properties.VariableNames = strcat(replace(untime{i},'-','_'),'_',tmptbl.Properties.VariableNames);
        widetbl = cat(2,widetbl,tmptbl);
    end
end

untime = replace(untime,'-','_');

withinvar = extractBefore(modelspec,'~');

modelspec_rm = strcat(untime,'_',withinvar,',');
modelspec_rm = cat(2,modelspec_rm{:}); modelspec_rm(end) = [];
modelspec_rm = cat(2,modelspec_rm,'~',strcat(untime{1},'_',extractAfter(modelspec,'~')));

rm = fitrm(widetbl,modelspec_rm,'WithinModel',withinmodel);




