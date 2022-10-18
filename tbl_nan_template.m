function tblout = tbl_nan_template(template)

tblout = table('Size',size(template),'VariableNames',template.Properties.VariableNames,...
    'VariableTypes',varfun(@class,template,'OutputFormat','cell'));
tblout = standardizeMissing(tblout,0);

vartypes = varfun(@class,template,'OutputFormat','cell');

tblout(:,strcmpi(vartypes,'cell')) = {''};


