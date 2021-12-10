function [torf] = isfield_nest(structin,field)

newfield = strsplit(field,'.');

fieldtrue = [];
for c = 1:length(newfield)
    if c > 1
        try
            newstruct = getfield_nest(structin,strjoin(newfield(1:(c-1)),'.'));
        catch err
            if strcmpi(err.identifier,'MATLAB:nonExistentField')
                torf = false;
                return
            end
        end
    else
        newstruct = structin;
    end
    if isstruct(newstruct)
    fieldtrue(c) = isfield(newstruct,newfield{c});
    elseif istable(newstruct)
        fieldtrue(c) = ismember(newfield{c},newstruct.Properties.VariableNames); 
    end
    %newstruct = getfield_nest([streval '.(''' newfield{c} ''')'];
end

torf = all(fieldtrue);

%eval(['fieldout = structin' streval ';']);