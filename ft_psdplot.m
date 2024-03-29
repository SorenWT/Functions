function [pxx,f]=ft_psdplot(data,chansel,range)
% Function for quickly plotting power spectra for fieldtrip structures
% Concatenates trials and calls psdplot

cfg = []; cfg.channel = chansel;
data = ft_selectdata(cfg,data);

catdata = cat(2,data.trial{:});

if nargin > 2
   [pxx,f]=psdplot(catdata,data.fsample,range)
else
    [pxx,f]=psdplot(catdata,data.fsample);
end