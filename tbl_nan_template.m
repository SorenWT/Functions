function tblout = tbl_nan_template(template,siz)

if nargin <2
   siz = size(template); 
end

tblout = table('Size',siz,'VariableNames',template.Properties.VariableNames,...
    'VariableTypes',varfun(@class,template,'OutputFormat','cell'));
tblout = standardizeMissing(tblout,0);

vartypes = varfun(@class,template,'OutputFormat','cell');

tblout(:,strcmpi(vartypes,'cell')) = {''};


