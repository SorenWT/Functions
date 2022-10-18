function tcat = tblhorzcat(tblsin)

for i = 1:length(tblsin)
    vnames{i} = tblsin{i}.Properties.VariableNames;
end

unvnames = [];
tcat = table;

for i = 1:length(tblsin)
    for ii = 1:length(vnames{i})
        if ~any(strcmpi(unvnames,vnames{i}{ii}))
            unvnames{end+1} = vnames{i}{ii};
            tcat.(vnames{i}{ii}) = tblsin{i}.(vnames{i}{ii});
        else
            tmp = tcat.(vnames{i}{ii});
            if iscell(tmp) && all(strcmpi(tmp,tblsin{i}.(vnames{i}{ii})));
            elseif isnumeric(tmp) && all(tmp==tblsin{i}.(vnames{i}{ii}));
            else
                unvnames{end+1} = vnames{i}{ii};
                newvname = [vnames{i}{ii} '_' num2str(sum(strcmpi(vnames{i}{ii},unvnames)))];
                tcat.(newvname) = tblsin{i}.(vnames{i}{ii});
            end
            
        end
    end
end