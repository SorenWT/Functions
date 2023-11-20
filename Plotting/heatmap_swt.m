function heatmap_swt(toplot,plotmask,xvals,yvals)

if ~exist('plotmask','var')
   plotmask = ones(size(toplot)); 
end

if ~exist('xvals','var')
    xvals = 1:size(toplot,2);
end

if ~exist('yvals','var')
    yvals = 1:size(toplot,1);
end

toplot = squeeze(toplot);
plotmask = squeeze(plotmask);


%plotdata = meandata.(fields{c});
plotdata = struct;
plotdata.freq = yvals; plotdata.time = xvals;
%plotdata.fourierspctrm(find(plotdata.fourierspctrm < 0)) = min(min(min(min(plotdata.fourierspctrm(find(plotdata.fourierspctrm > 0))))));
plotdata.fourierspctrm = permute(toplot,[3 1 2]);
plotdata.fourierspctrm = repmat(plotdata.fourierspctrm,306,1,1);
plotdata.mask = logical(permute(plotmask,[3 1 2]));
plotdata.fourierspctrm_dimord = 'rpt_chan_freq_time';
plotdata.dimord = 'chan_freq_time';
plotdata.label = cellstr(num2str([1:306]'));

cfg = []; cfg.parameter = 'fourierspctrm'; cfg.baseline = 'no';
cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
cfg.layout = 'neuromag306mag.lay'; cfg.channel = 1; cfg.latency = [min(xvals) max(xvals)];
cfg.interactive = 'no';  cfg.masknans = 'no';

ft_singleplotTFR(cfg,plotdata);

set(gca,'XLim',[min(xvals)-0.5*mean(diff(xvals)) max(xvals)+0.5*mean(diff(xvals))],...
    'YLim',[min(yvals)-0.5*mean(diff(yvals)) max(yvals)+0.5*mean(diff(yvals))])
FixAxes(gca,16)
%hold on
%delete(findobj('parent',gca,'type','image'))

%db_plot = 10*log10(plotdata.fourierspctrm(:,1,:,plotdata.time>-0.22)./nanmean(plotdata.fourierspctrm(:,1,:,plotdata.time<0),4));
%db_plot(imag(db_plot)~=0) = NaN;
%db_plot = squeeze(nanmean(db_plot,1));
% s=pcolor(xvals,yvals,toplot);
% s.EdgeAlpha = 0;
% colorbar
% hold on
% 
% uistack(findobj('parent',gca,'type','contour'),'top')
%colormap(jet)