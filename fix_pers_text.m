function persout = fix_pers_text(pers)

pers = rmfield(pers,'prettyitems');

if iscell(pers)
    for i = 1:length(pers)
        persout{i} = fix_pers_text(pers{i});
    end
else
    for i = 1:length(pers.items)
        goodtext = pers.items{i};
        goodtext = regexprep(goodtext, '([A-Z])', ' ${lower($1)}');
        % capitalize I's
        goodtext = regexprep(goodtext, '\si\s', ' I ');
        goodtext(1) = [];
        goodtext(1) = upper(goodtext(1));
        
        goodtext = regexprep(goodtext,'_\s',', ');
        goodtext = regexprep(goodtext,'_[^\s]',[char(39) '$0']);
        goodtext = erase(goodtext,'_');
        [s] = input([goodtext newline],'s');
        if ~isempty(s)
            goodtext = s;
        end
        
        pers.prettyitems{i} = goodtext;
    end
    persout = pers;
end