function tbldiff = tablediff(tblin,diffvar,diffvarorder)

if nargin < 3
u = unique(tblin.(diffvar)); 
else
    u = diffvarorder;
end

tbldiff = tbl_nan_template(tblin,[0 width(tblin)]);
varnames = tblin.Properties.VariableNames;


for i = 1:length(diffvarorder)-1
    tmp2 = tblin(strcmpi(tblin.(diffvar),u{i+1}),:); 
    tmp1 = tblin(strcmpi(tblin.(diffvar),u{i}),:); 
    
    tbldiff_tmp = table;
    for ii = 1:length(varnames)
       if isnumeric(tblin.(varnames{ii}))
          tbldiff_tmp.(varnames{ii}) = tmp2.(varnames{ii}) - tmp1.(varnames{ii});
       elseif iscell(tblin.(varnames{ii}))
          tbldiff_tmp.(varnames{ii}) = strcat(tmp2.(varnames{ii}),tmp1.(varnames{ii}));
       end
    end
    tbldiff = cat(1,tbldiff,tbldiff_tmp);
end
