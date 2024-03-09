function merged = mergestructs(structsin,tomerge,usetemplate)

if nargin < 2
    tomerge = [];
end

if nargin < 3
    usetemplate = 1;
end

if ~iscell(structsin)
    fields = cell_unpack(fieldnames_recurse(structsin(1)));
    if isempty(tomerge)
        tomerge = 1:length(fields);
    elseif islogical(tomerge)
        tomerge = find(tomerge);
    elseif iscell(tomerge)
        tomerge = match_str(fields,tomerge);
    end
    if usetemplate
        merged = structsin(1);
    else
        merged = struct;
    end
    tomerge = horz(tomerge);

    for c = tomerge
        dimn = size(getfield_nest(structsin(1),fields{c}));
        if ismember(1,dimn)
            dimn = find(dimn == 1,1);
        else
            dimn = length(dimn)+1;
        end
        alldata = cell(1,length(structsin));
        for cc = 1:length(structsin)
            alldata{cc} = getfield_nest(structsin(cc),fields{c});
        end
        merged = assignfield_nest(merged,fields{c},cat(dimn,alldata{:}));
    end
else
    fields = cell_unpack(fieldnames_recurse(structsin{1}));
    if isempty(tomerge)
        tomerge = 1:length(fields);
    elseif islogical(tomerge)
        tomerge = find(tomerge);
    elseif iscell(tomerge)
        tomerge = match_str(fields,tomerge);
    end
    if usetemplate
        merged = structsin{1};
    else
        merged = struct;
    end
    tomerge = horz(tomerge);
    for c = tomerge
        dimn = size(getfield_nest(structsin{1},fields{c}));
        if ismember(1,dimn)
            dimn = find(dimn == 1,1);
        else
            dimn = length(dimn)+1;
        end
        try
            alldata = cell(1,length(structsin));
            for cc = 1:length(structsin)
                alldata{cc} = getfield_nest(structsin{cc},fields{c});
            end
            disp([fields{c} ' merged'])
            merged = assignfield_nest(merged,fields{c},cat(dimn,alldata{:}));
        catch
            warning(['Failed to merge field "' fields{c} '"'])
        end
    end
end
