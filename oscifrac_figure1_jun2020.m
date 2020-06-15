
addpath(genpath('/project/def-gnorthof/sorenwt/MATLAB'))
addpath('/project/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

load('settings_camcan_1Hz.mat')

cd /scratch/sorenwt/camcan/Preprocessed/Task/Epoched/

files = dir('*tf.mat');

m = matfile(files(1).name);

fields = {'mixd','osci','frac'};

for c = 1:length(fields)
    meandata.(fields{c}) = m.(fields{c});
end

%% Reading in data
donefiles = [];
for i = 1:length(files)
    try
        %if ~ismember({files(i).name},donefiles)
            m = matfile(files(i).name);
            for c = 1:3
                tmp = m.(fields{c});
                tmp.fourierspctrm(find(tmp.fourierspctrm < 0)) = NaN;
                meandata.(fields{c}).fourierspctrm(i,:,:,:) = nanmean(tmp.fourierspctrm,1);
            end
            donefiles = [donefiles {files(i).name}];
        %end
    catch
    end
end

for c = 1:length(fields)
   meandata.(fields{c}).fourierspctrm = permute(meandata.(fields{c}).fourierspctrm,[1 3 2 4]); 
   meandata.(fields{c}).fourierspctrm_dimord = meandata.(fields{c}).dimord;
   meandata.(fields{c}).grad = settings.datasetinfo.grad;
end


%% Analysis and statistics

for c = 1:3
    post = find(meandata.(fields{c}).time >= 0);
    meanpost.(fields{c}) = meandata.(fields{c});
    meanpost.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(:,:,:,post);
    meanpost.(fields{c}).time = meandata.(fields{c}).time(post);
    
    bl = find(meandata.(fields{c}).time < 0);
    bl = bl((end-length(post)+1):end);
    meanbl.(fields{c}) = meandata.(fields{c});
    meanbl.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(:,:,:,bl);
    meanbl.(fields{c}).time = meandata.(fields{c}).time(bl);
end

post = cell(1,size(meanpost.frac.fourierspctrm,1));
bl = post;
postint = post;
blint = bl;


parfor i = 1:size(meanpost.frac.fourierspctrm,1)
    for ii = 1:size(meanpost.frac.fourierspctrm,2)
        for iii = 1:size(meanpost.frac.fourierspctrm,4)
            tmp = log10(meanpost.frac.freq);
            pow = log10(squeeze(meanpost.frac.fourierspctrm(i,ii,:,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            post{i}(1,ii,iii) = -p(1);
            postint{i}(1,ii,iii) = p(2);
            
            tmp = log10(meanbl.frac.freq);
            pow = log10(squeeze(meanbl.frac.fourierspctrm(i,ii,:,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            bl{i}(1,ii,iii) = -p(1);
            blint{i}(1,ii,iii) = p(2);
        end
    end
end

pledata.post = cat(1,post{:});
pledata.bl = cat(1,bl{:});

intdata.post = cat(1,postint{:});
intdata.bl = cat(1,blint{:});

for c = 1:3
    cfg = []; cfg.channel = {'MEG'}; cfg.avgoverchan = 'yes'; cfg.latency = [0 0.75];
    cfg.frequency = 'all'; cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_actvsblT';
    cfg.correctm = 'cluster'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum';
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; cfg.numrandomization = 2000;
    cfg.parameter = 'fourierspctrm';
    
    ntrials = size(meanpost.(fields{c}).fourierspctrm,1);
    design  = zeros(2,2*ntrials);
    design(1,1:ntrials) = 1;
    design(1,ntrials+1:2*ntrials) = 2;
    design(2,1:ntrials) = [1:ntrials];
    design(2,ntrials+1:2*ntrials) = [1:ntrials];
    
    cfg.design = design;
    cfg.ivar = 1;
    cfg.uvar = 2;
    
    cfg.parpool = 48;
    
    meanbl.(fields{c}).time = meanpost.(fields{c}).time;
    
    stats{c} = ft_freqstatistics(cfg,meanpost.(fields{c}),meanbl.(fields{c}));
end

% Cluster stats for PLE and broadband

datasetinfo = settings.datasetinfo;

opts = []; opts.nrand = 2000; opts.parpool = 48; opts.minnbchan = 0;
%opts.external = cfg.design;

stats_ple = EasyClusterCorrect({permute(pledata.post,[2 1 3]) repmat(mean(pledata.bl,3)',1,1,38)},...
    datasetinfo,'ft_statfun_fast_signrank',opts);

%stats_bb = EasyClusterCorrect({permute(squeeze(mean(meanpost.frac.fourierspctrm,3)),[2 1 3]) repmat(mean(squeeze(mean(meanbl.frac.fourierspctrm,3)),3)',1,1,38)},...
%    datasetinfo,'ft_statfun_fast_signrank',opts);

stats_int = EasyClusterCorrect({permute(intdata.post,[2 1 3]) repmat(mean(intdata.bl,3)',1,1,38)},...
    datasetinfo,'ft_statfun_fast_signrank',opts);

% for divisive baseline/gain model, set gc1 = 0, gc2 = 1
% for subtractive baseline, set gc1 = 1; gc2 = 0

% gc1 = 1;
% gc2 = 0; 
% 
% for c = find(extractfield(stats{1}.posclusters,'prob') < 0.05)
%     statmask = squeeze(stats{1}.posclusterslabelmat==c);
%     % for ple and broadband fractal, include those time points where a
%     % number of frequencies greater than the median value in the cluster
%     % are significant
%     sum_statmask = sum(statmask,1);
%     %statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
%     statmask_time = sum_statmask > 0;
%     
%     newindx = c;
%     % for ple and broadband fractal, include those time points where a
%     % number of frequencies greater than the median value in the cluster
%     % are significant
%     sum_statmask = sum(statmask,1);
%     statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
%     
%     tmp = mean(permute(permute(squeeze((mean(meanpost.mixd.fourierspctrm,2)...
%         -mean(mean(meanbl.mixd.fourierspctrm,4),2))./(gc1+gc2.*mean(mean(meanbl.mixd.fourierspctrm,4),2))),[2 3 1]).*statmask,[3 1 2]),3);
%     meanmix(:,newindx) = mean(tmp,2); % baseline correct, sum over time and freq for significant cluster
%     
%     tmp = mean(permute(permute(squeeze((mean(meanpost.osci.fourierspctrm,2)...
%         -mean(mean(meanbl.osci.fourierspctrm,4),2))./(gc1+gc2.*mean(mean(meanbl.osci.fourierspctrm,4),2))),[2 3 1]).*statmask,[3 1 2]),3);
%     meanosci(:,newindx) = mean(tmp,2);
%     
%     meanfrac(:,newindx) = mean(squeeze((mean(mean(meanpost.frac.fourierspctrm,2),3)...
%         -mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2))./(gc1+gc2*mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2))).*statmask_time,2); %use only the time statmask here, not freq - broadband power
%     
%     meanfrac2(:,newindx) = mean(mean(permute(permute(squeeze((mean(meanpost.frac.fourierspctrm,2)...
%         -mean(mean(meanbl.frac.fourierspctrm,4),2))./(gc1+gc2.*mean(mean(meanbl.frac.fourierspctrm,4),2))),[2 3 1]).*statmask,[3 1 2]),2),3);
%     
%     meanple(:,newindx) = mean(squeeze((mean(pledata.post,2)-mean(mean(pledata.bl,2),3))./(gc1+gc2.*mean(mean(pledata.bl,2),3))).*statmask_time,2);
% end
% 
% for c = find(extractfield(stats{1}.negclusters,'prob') < 0.05)
%     statmask = squeeze(stats{1}.negclusterslabelmat==c);
%     % for ple and broadband fractal, include those time points where a
%     % number of frequencies greater than the median value in the cluster
%     % are significant
%     sum_statmask = sum(statmask,1);
%     %statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
%     statmask_time = sum_statmask > 0;
%     
%     newindx = c+length(find(extractfield(stats{1}.posclusters,'prob') < 0.05));
%     % for ple and broadband fractal, include those time points where a
%     % number of frequencies greater than the median value in the cluster
%     % are significant
%     sum_statmask = sum(statmask,1);
%     statmask_time = sum_statmask > median(sum_statmask(find(sum_statmask>0))); 
%     
% 
%     tmp = mean(permute(permute(squeeze((mean(meanpost.mixd.fourierspctrm,2)...
%         -mean(mean(meanbl.mixd.fourierspctrm,4),2))./(gc1+gc2.*mean(mean(meanbl.mixd.fourierspctrm,4),2))),[2 3 1]).*statmask,[3 1 2]),3);
%     meanmix(:,newindx) = mean(tmp,2); % baseline correct, sum over time and freq for significant cluster
%     
%     tmp = mean(permute(permute(squeeze((mean(meanpost.osci.fourierspctrm,2)...
%         -mean(mean(meanbl.osci.fourierspctrm,4),2))./(gc1+gc2.*mean(mean(meanbl.osci.fourierspctrm,4),2))),[2 3 1]).*statmask,[3 1 2]),3);
%     meanosci(:,newindx) = mean(tmp,2);
%     
%     meanfrac(:,newindx) = mean(squeeze((mean(mean(meanpost.frac.fourierspctrm,2),3)...
%         -mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2))./(gc1+gc2.*mean(mean(mean(meanbl.frac.fourierspctrm,4),3),2))).*statmask_time,2); %use only the time statmask here, not freq - broadband power
%     
%     meanfrac2(:,newindx) = mean(mean(permute(permute(squeeze((mean(meanpost.frac.fourierspctrm,2)...
%         -mean(mean(meanbl.frac.fourierspctrm,4),2))./(gc1+gc2.*mean(mean(meanbl.frac.fourierspctrm,4),2))),[2 3 1]).*statmask,[3 1 2]),2),3);
%     
%     meanple(:,newindx) = mean(squeeze((mean(pledata.post,2)-mean(mean(pledata.bl,2),3))./(gc1+gc2.*mean(mean(pledata.bl,2),3))).*statmask_time,2);
% end

settings = NA_alpha_pf(settings);

for i = 2:settings.nfreqs
    tmp = meandata.osci.fourierspctrm;
    %tmp = abs(tmp);
    tmp(find(tmp <=0)) = NaN;
    for q = 1:size(meandata.osci.fourierspctrm,1)
        fbands{i}.post(q,:,:) = squeeze(nanmean(tmp(q,:,...
            intersect(find(meandata.osci.freq >= settings.tfparams.fbands{q,i}(1)),...
            find(meandata.osci.freq <= settings.tfparams.fbands{q,i}(2))),post),3));
        fbands{i}.bl(q,:,:) = squeeze(nanmean(tmp(q,:,...
            intersect(find(meandata.osci.freq >= settings.tfparams.fbands{q,i}(1)),...
            find(meandata.osci.freq <= settings.tfparams.fbands{q,i}(2))),bl),3));
    end
    stats_fbands{i} = EasyClusterCorrect({permute(fbands{i}.post,[2 1 3]) repmat(mean(fbands{i}.bl,3)',1,1,38)},...
        datasetinfo,'ft_statfun_fast_signrank',opts);
end


%% Plotting the figure

p = panel('no-manage-font');

p.pack('h',{1/2 1/2})


for c = 2:3
    p(c-1).select()
    
    plotdata = meandata.(fields{c});
    %plotdata.fourierspctrm(find(plotdata.fourierspctrm < 0)) = min(min(min(min(plotdata.fourierspctrm(find(plotdata.fourierspctrm > 0))))));
    plotdata.fourierspctrm(find(plotdata.fourierspctrm < 0)) = NaN;
    plotdata.fourierspctrm(:,1,:,:) = nanmean(plotdata.fourierspctrm,2); % make the first channel the mean so you can use plotting mask
    plotdata.mask = logical(cat(3,zeros(1,size(stats{c}.mask,2),length(plotdata.time)-size(stats{c}.mask,3)),stats{c}.mask));
    plotdata.fourierspctrm_dimord = 'rpt_chan_freq_time';
    plotdata.dimord = 'rpt_chan_freq_time';
    
    cfg = []; cfg.parameter = 'fourierspctrm'; cfg.baseline = [-Inf 0]; cfg.baselinetype = 'db';
    cfg.maskparameter = 'mask'; cfg.maskstyle = 'outline';
    cfg.layout = 'neuromag306mag.lay'; cfg.channel = 1; cfg.latency = [0 0.75];
    cfg.interactive = 'no';  cfg.masknans = 'no';
    ft_singleplotTFR(cfg,plotdata);
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    set(gca,'XLim',[-0.2 0.75])
    hold on
    line([0 0],get(gca,'YLim'),'Color',[0 0 0],'LineWidth',2)
    title(fields{c},'FontSize',14)
    ax(c) = gca;
    colormap(jet)
end

Set_Clim(ax,[-4 4])

%part 2

figure

p = panel('no-manage-font')

p.pack('v',{35 35 30})

p(1).pack('h',{1/6 1/3 1/3 1/6})
p(2).pack('h',repmat({1/(settings.nfreqs-1)},1,settings.nfreqs-1))
set(gcf,'units','normalized','position',[0 0 1 1])

%p(1,4).pack('v',{50 50})

plotdata = 100*(pledata.post-nanmean(pledata.bl,3))./nanmean(pledata.bl,3);
pleindx = nansum(nansum(plotdata.*repmat(permute(stats_ple.mask,[3 1 2]),49,1,1),2),3);

plot_tc_topo(meanpost.mixd.time*1000,permute(plotdata,[2 3 1]),settings,p,{1 2})
%stdshade(meanpost.mixd.time,squeeze(mean(pledata.post,2)),'k',0.15,2,'sem')
%Plot_sigmask(gca,stats_ple.mask,'cmapline','LineWidth',4)
p(1,2,1).select()
hold on
xlabel('Time (s)')
ylabel('Percent change from prestim (%)')
title('Fractal PLE')
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_ple.negclusters.prob,0,18);
delete(H(1));

p(1,3).marginleft = 25;

plotdata = 100*(10.^intdata.post-nanmean(10.^intdata.bl,3))./nanmean(10.^intdata.bl,3);
intindx = nansum(nansum(plotdata.*repmat(permute(stats_int.mask,[3 1 2]),49,1,1),2),3);

plot_tc_topo(meanpost.mixd.time*1000,permute(plotdata,[2 3 1]),settings,p,{1 3})
%Plot_sigmask(gca,stats_int.mask,'cmapline','LineWidth',4)
p(1,3,1).select()
FixAxes(gca,12)
xlabel('Time (s)')
%ylabel('Percent change from prestim (%)')
title('Fractal intercept')
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_int.negclusters.prob,0,18);
delete(H(1));


for i = 2:settings.nfreqs
   plotdata = 100*(fbands{i}.post-nanmean(fbands{i}.bl(:,:,end-4:end),3))./nanmean(fbands{i}.bl(:,:,end-4:end),3);
   fbandindx(:,i) = nansum(nansum(plotdata.*repmat(permute(stats_fbands{i}.mask,[3 1 2]),49,1,1),2),3);
   
   plot_tc_topo(meanpost.mixd.time*1000,permute(plotdata,[2 3 1]),settings,p,{2 i-1},'avgfunc',@nanmedian);
   p(2,i-1,1).select()
   FixAxes(gca,12)
   xlabel('Time (s)')
   if i ==2
    ylabel('Percent change from prestim (%)')
   end
   title(settings.tfparams.fbandnames{i})
   set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
   ax = gca;
   H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_int.negclusters.prob,0,18);
   delete(H(1));
end

p(3).select()
indxmat = abs([pleindx intindx fbandindx(:,2:end)]);
anovap = friedman([pleindx intindx fbandindx(:,2:end)]);
for i = 1:size(indxmat,2)
    for ii = 1:i
        mcp(i,ii) = signrank(indxmat(:,i),indxmat(:,ii)); 
    end
end

bar(mean(indxmat,1))
hold on
errorbar(1:8,mean(indxmat,1),std(indxmat,[],1)./sqrt(size(indxmat,1)),'LineStyle','none','Color',[0 0 0])
%H = notBoxPlot(indxmat);
%delete([H(:).data])
set(gca,'XTickLabel',[{'PLE','Fractal intercept'} settings.tfparams.fbandnames(2:end)])

colormap(lkcmap2)

savefig('Oscifrac_1b.fig'); export_fig


% p(1,4,2).select()
% stdshade(meanpost.mixd.time,squeeze(mean(mean(meanpost.frac.fourierspctrm,2),3)),'k',0.15,2,'sem')
% Plot_sigmask(gca,stats_bb.mask,'cmapline','LineWidth',3)
% FixAxes(gca,12)
% xlabel('Time (s)')
% ylabel('Fractal broadband power')
% set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)])




% 
% 
% p(2).pack('h',repmat({1/length(mdl)},1,length(mdl)));
% 
% for i = 1:length(mdl)
%    p(2,i).select()
%    scatter(1:3,mdl{i}.Coefficients.Estimate(2:end),36,[0 0 1],'x','LineWidth',2)
%    hold on
%    %scatter(1:3,partr(i,:),36,[1 0 0],'o','LineWidth',2)
%    er = errorbar(1:3,mdl{i}.Coefficients.Estimate(2:end),mdl{i}.Coefficients.SE(2:end)*1.96,...
%        'LineStyle','none','LineWidth',2,'Color',[0 0 1],'HandleVisibility','off');
%    xl = get(gca,'XLim');
%    line(xl+[-0.1 0.1],[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
%    set(gca,'XLim',xl + [-0.1 0.1],'XTickLabel',{'Oscilatory Power','PLE','Fractal Broadband'})
%    %legend({'Regression Coefficient','Partial R^2'})
%    %ylabel('Regression Coefficient')
%    FixAxes(gca,14)
%    fix_xticklabels(gca,0.1,{'FontSize',14}) 
%    ylabel('Coefficient')
% end

% p(1).de.marginleft = 30;
% p.marginright = 20;
% p.marginleft = 18;
% p.margintop = 8;



