function filled = scales_fill_nan(tofill,template)

filled = tofill; 

template_vars = template.Properties.VariableNames;
tofill_vars = tofill.Properties.VariableNames;

tmp = 