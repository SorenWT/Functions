function [torf] = isfield_nest(structin,field)

newfield = strsplit(field,'.');

fieldtrue = [];
for c = 1:length(newfield)
    if c > 1
        newstruct = getfield_nest(structin,strjoin(newfield(1:(c-1)),'.'));
    else
        newstruct = structin;
    end
    fieldtrue(c) = isfield(newstruct,newfield{c});
    %newstruct = getfield_nest([streval '.(''' newfield{c} ''')'];
end

torf = all(fieldtrue);

%eval(['fieldout = structin' streval ';']);