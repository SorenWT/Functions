function ft_topoplot_vec(layout,vect,label,varargin)

vect = vert(vect);

argsin = varargin;

tlock = [];
tlock.avg = vect;
tlock.dimord = 'chan_time';
tlock.label = label;
tlock.time = 1;

if ischar(layout)
    cfg = []; cfg.layout = layout;
    layout = ft_prepare_layout(cfg);
end
cfg = [];
cfg.layout = layout;
cfg.interpolateovernan = 'yes';
cfg.lay = layout;
cfg.channel = label;
cfg.comment = 'no'; cfg.interactive = 'no';
cfg.marker = 'labels';
%cfg.markersymbol = '.';
%cfg.markersize = 3;

for i = 1:2:length(argsin)
    cfg.(argsin{i}) = argsin{i+1}; 
end

ft_topoplotER(cfg,tlock);