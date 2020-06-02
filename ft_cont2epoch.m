function [data] = ft_cont2epoch(cont_data,epochlen)

if nargin < 2
   epochlen = 2; 
end

cfg = []; cfg.event = 1:epochlen*cont_data.fsample:floor(length(cont_data.time{1}));
cfg.event(end) = []; cfg.epoch = [0 (epochlen*cont_data.fsample)-1];
data = ft_epoch(cfg,cont_data);


