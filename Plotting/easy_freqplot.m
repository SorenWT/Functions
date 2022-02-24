function cbar = easy_freqplot(freq)
% Inputs:
%   freq: a struct with the following fields
%       toplot: a nfreqs x ntime array with the actual data to plot. This
%       should be already log-transformed, baseline-corrected, etc.
%       freq: a vector of length nfreqs with the frequencies corresponding
%          to the toplot field
%       time: a vector of length ntimes with the time points corresponding
%          to the toplot field
%       mask (optional): a nfreqs x ntimes array with the significance mask

freq.toplot = permute(freq.toplot,[3 4 1 2]);
freq.mask = permute(freq.mask,[3 1 2]);


freq = setdefault(freq,'plotunits','dB');

load('neuromag306all_helmet.mat')

plotdata = struct;
plotdata.freq = freq.freq;
plotdata.time = freq.time;
%plotdata.fourierspctrm(find(plotdata.fourierspctrm < 0)) = min(min(min(min(plotdata.fourierspctrm(find(plotdata.fourierspctrm > 0))))));
plotdata.fourierspctrm = repmat(freq.toplot,10,306,1,1);
%plotdata.fourierspctrm = squeeze(mean(plotdata.fourierspctrm,1));
plotdata.mask = repmat(freq.mask,306,1,1);
plotdata.fourierspctrm_dimord = 'rpt_chan_freq_time';
plotdata.dimord = 'rpt_chan_freq_time';
plotdata.label = layout.label(1:306);

cfg = []; cfg.parameter = 'fourierspctrm'; cfg.baseline = 'no';
cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
cfg.layout = 'neuromag306mag.lay'; cfg.channel = 1; cfg.latency = 'all';
cfg.interactive = 'no';  cfg.masknans = 'no';
ft_singleplotTFR(cfg,plotdata);
hold on
delete(findobj('parent',gca,'type','image'))

%db_plot = freq.toplot;
%db_plot = 10*log10(plotdata.fourierspctrm(:,1,:,plotdata.time>-0.22)./nanmean(plotdata.fourierspctrm(:,1,:,plotdata.time<0),4));
%db_plot(imag(db_plot)~=0) = NaN;
%db_plot = squeeze(nanmean(db_plot,1));
s=pcolor(freq.time,freq.freq,squeeze(freq.toplot));
s.EdgeAlpha = 0;
colorbar
hold on

maskline = findobj('parent',gca,'type','contour');
uistack(maskline,'top')
%maskline.XData = maskline.XData+0.5*mean(diff(freq.time));
%maskline.YData = maskline.YData+0.5*mean(diff(freq.freq));



xlabel('Time (s)')
ylabel('Frequency (Hz)')
hold on
colormap(jet)
cbar = colorbar;
cbar.Label.String = 'ERSP (dB)';
cbar.Label.FontSize = 14;


