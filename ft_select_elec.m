function elecout = ft_select_elec(chans,elec)

if iscell(chans)
    chanindx = match_str(elec.label,chans);
elseif isnumeric(chans)
    chanindx = chans;
end

elecout = elec;
if isfield(elecout,'chanpos')
    elecout.chanpos = elecout.chanpos(chanindx,:);
end
if isfield(elecout,'elecpos')
    elecout.elecpos = elecout.elecpos(chanindx,:);
end
if isfield(elecout,'label')
    elecout.label = elecout.label(chanindx);
end
if isfield(elecout,'chantype')
    elecout.chantype = elecout.chantype(chanindx);
end
if isfield(elecout,'chanunit')
    elecout.chanunit = elecout.chanunit(chanindx);
end
