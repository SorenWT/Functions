function tableout = set_refcat(tablein,tablevar,refcat)

tablein.(tablevar) = categorical(tablein.(tablevar));
cats = categories(tablein.(tablevar));
tablein.(tablevar) = reordercats(tablein.(tablevar),[{refcat} horz(setdiff(cats,refcat))]);

tableout = tablein;

