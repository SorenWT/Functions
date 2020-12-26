
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
   meandata.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(1:length(files),:,:,:);
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
settings.nfreqs = 6;
settings.tfparams.fbands = repmat({[] [2 4] [4 8] [8 13] [13 30] [30 85]},49,1);

freq = meandata.mixd.freq;
fbands = cell(1,7);
for i = 2:settings.nfreqs
    tmp = meandata.osci.fourierspctrm;
    %tmp = abs(tmp);
    tmp(tmp<=0) = NaN;
    for q = 1:size(meandata.osci.fourierspctrm,1)
        tmpfrange = intersect(find(meandata.osci.freq >= settings.tfparams.fbands{q,i}(1)),...
            find(meandata.osci.freq <= settings.tfparams.fbands{q,i}(2)));
        fbands{i}.post(q,:,:) = squeeze(simps(freq(tmpfrange),tmp(q,:,tmpfrange,postindx),3));
        fbands{i}.bl(q,:,:) = squeeze(simps(freq(tmpfrange),tmp(q,:,tmpfrange,blindx),3));
    end
    stats_fbands{i} = EasyClusterCorrect({permute(fbands{i}.post,[2 1 3]) repmat(mean(fbands{i}.bl,3)',1,1,38)},...
        datasetinfo,'ft_statfun_fast_signrank',opts);
end

load('/group/northoff/share/camcan/Rest/Pseudo_epochs/IRASAtf/camcan_irasatf_rest_tresoutputs.mat');
outputs.meas = outputs.meas([2 1 3:7]); outputs.data = outputs.data(:,:,[2 1 3:7],:);

%% Plotting the figure

figure

p = panel('no-manage-font');

set(gcf,'units','normalized','position',[0.1 0 0.8 1])

p.pack('v',{1/2 1/2})
p(1).pack('h',{0.25 0.5 0.25})
p(2).pack('h',{1/2 1/2})


for c = 1:3
    if c == 1
        p(1,2).select()
    else
        p(2,4-c).select()
    end
    
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

Set_Clim(ax,[-3 3])

set(p(1,2).axis,'fontsize',14)
set(p(2,1).axis,'fontsize',14)
set(p(2,2).axis,'fontsize',14)

p.marginleft = 18;
p(2,2).marginleft = 40;
p(1).marginbottom = 40;
p.margintop = 8;
cbars = findall(gcf,'type','colorbar');
cbars(2).Label.String = 'ERSP (dB)';
cbars(2).Label.FontSize = 14;
set(gcf,'color','w')

savefig('Fig1a_erspplot.fig'); export_fig('Fig1a_erspplot.png','-m4')

%part 2

figure

p = panel('no-manage-font');

p.pack('h',{1/2 1/2})
p(1).pack('h',{1/2 1/2})
p(2).pack(2,3)
p.marginbottom = 22; p.marginleft = 22;
p.marginright = 12;
p(2).marginleft = 25;
p.margintop = 10;
set(gcf,'units','normalized','position',[0 0 1 1])

%p(1,4).pack('v',{50 50})

plotdata = effsizetc.ple;
%plotdata = 100*(pledata.post-nanmean(pledata.bl,3))./nanmean(pledata.bl,3);
%pleindx = nansum(nansum(plotdata.*repmat(permute(stats_ple.mask,[3 1 2]),49,1,1),2),3);

plot_tc_topo(meanpost.mixd.time*1000,plotdata+rand(102,38,49)*0.001,settings,p,{1 2},'color',{'r'})
%stdshade(meanpost.mixd.time,squeeze(mean(pledata.post,2)),'k',0.15,2,'sem')
%Plot_sigmask(gca,stats_ple.mask,'cmapline','LineWidth',4)
p(1,2,1).select()
hold on
xlabel('Time (s)')
ylabel('Percent change from baseline (%)')
%ylabel('Pseudo effect size')
title('Fractal PLE')
FixAxes(gca,16)
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_ple.negclusters.prob,0,18);
delete(H(1));


%plotdata = 100*(10.^intdata.post-nanmean(10.^intdata.bl,3))./nanmean(10.^intdata.bl,3);
plotdata = effsizetc.int;
%intindx = nansum(nansum(plotdata.*repmat(permute(stats_int.mask,[3 1 2]),49,1,1),2),3);

plot_tc_topo(meanpost.mixd.time*1000,plotdata+rand(102,38,49)*0.001,settings,p,{1 1},'color',{'r'})
%Plot_sigmask(gca,stats_int.mask,'cmapline','LineWidth',4)
p(1,1,1).select()
FixAxes(gca,16)
xlabel('Time (s)')
%ylabel('Percent change from prestim (%)')
title('Fractal intercept')
ylabel('Percent change from baseline (%)')
%ylabel('Pseudo effect size')
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_int.negclusters.prob,0,18);
delete(H(1));

%yl=Normalize_Ylim([p(1,1,1).axis,p(1,2,1).axis],0);
%Set_Ylim([p(1,1,1).axis,p(1,2,1).axis],[0 yl(2)]);


% 
tmpindx = [0 1 2 3 1 2 3];
for i = 2:settings.nfreqs
   plotdata = effsizetc.(fields{i+1});
   %plotdata = 100*(fbands{i}.post-nanmean(fbands{i}.bl(:,:,end-4:end),3))./nanmean(fbands{i}.bl(:,:,end-4:end),3);
   %fbandindx(:,i) = nansum(nansum(plotdata.*repmat(permute(stats_fbands{i}.mask,[3 1 2]),49,1,1),2),3);
   plot_tc_topo(meanpost.mixd.time*1000,plotdata+rand(102,38,49)*0.001,settings,p,{2,ceil((i-1)/3),tmpindx(i)},'avgfunc',@nanmedian);
   p(2,ceil((i-1)/3),tmpindx(i),1).select()
   FixAxes(gca,12)
   xlabel('Time (s)')
   if i ==2 || i == 5
    ylabel('Percent change from baseline (%)')
    %ylabel('Pseudo effect size')
   end
   title(settings.tfparams.fbandnames{i})
   set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
   ax = gca;
   H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_int.negclusters.prob,0,18);
   delete(H(1));
end

colormap(lkcmap2)


savefig('Fig1b_param_tc_cohen.fig'); export_fig('Fig1b_param_tc_cohen.png','-m4')

fbandnames = {'delta','theta','alpha','beta','gamma'};

%effsizefun = @(post,bl,mask) nansum(nansum((nanmean(post-nanmean(bl,3),1)./nanstd(post-nanmean(bl,3),[],1)).*mean(mask,1)));
%effsizefun = @(post,bl,mask) nansum(nansum((post-nanmean(bl,3))./nanstd(bl,[],3)));
effsizefun = @(post,bl,mask) nansum(nansum(abs(nanmedian(100*(post-nanmean(bl,3))./nanmean(bl,3),1)).*mean(mask,1),3),2);
effsizefun = @(post,bl,mask) nansum(nansum(abs(nanmean(post-nanmean(bl,3),1)./((repmat(nanstd(nanmean(bl,3),[],1),1,1,size(post,3))+nanstd(post,[],1))/2)).*mean(mask,1),3),2);

effsizedifffun = @(effsizefun,post1,bl1,mask1,post2,bl2,mask2) effsizefun(post1,bl1,mask1)-effsizefun(post2,bl2,mask2);

effsizetc = struct;
effsizetc.int = nanmean(10.^intdata.post-nanmean(10.^intdata.bl,3),1)./...
    ((repmat(nanstd(nanmean(10.^intdata.bl,3),[],1),1,1,38)+nanstd(10.^intdata.post,[],1))/2);
effsizetc.ple = nanmean(pledata.post-nanmean(pledata.bl,3),1)./...
    ((repmat(nanstd(nanmean(pledata.bl,3),[],1),1,1,38)+nanstd(pledata.post,[],1))/2);
for i = 2:settings.nfreqs
   %tmppost = log10(fbands{i}.post); tmppost(find(imag(tmppost)~=0)) = NaN;
   %tmpbl = log10(fbands{i}.post); tmpbl(find(imag(tmpbl)~=0)) = NaN;
   
   %outdata = outputs.data(:,:,i+1,:); outdata(outdata<=0) = NaN; outdata = log10(outdata);
   %effsizetc.(fbandnames{i-1}) = (tmppost-nanmean(tmpbl,3))./...
   %    nanstd(tmp,[],4);
    effsizetc.(fbandnames{i-1}) = nanmean(fbands{i}.post - nanmean(fbands{i}.bl,3),1)./...
        ((repmat(nanstd(nanmean(fbands{i}.bl,3),[],1),1,1,38)+nanstd(fbands{i}.post,[],1))/2);
end
effsizetc = structfun(@squeeze,effsizetc,'UniformOutput',false);

allstats = [{stats_int stats_ple} stats_fbands(2:end)];

fields = fieldnames(effsizetc);

for i = 1:length(fields)
    effsizesum(i) = nansum(nansum(abs(effsizetc.(fields{i})).*allstats{i}.mask));
end

effsizediff = ones(length(fields));
for i = 1:length(fields)
   if i == 1
       post1 = 10.^intdata.post;
       bl1 = repmat(nanmean(10.^intdata.bl,3),1,1,size(intdata.bl,3));
   elseif i == 2
       post1 = pledata.post;
       bl1 = repmat(nanmean(pledata.bl,3),1,1,size(pledata.bl,3));
   else
       post1 = fbands{i-1}.post;
       bl1 = repmat(nanmean(fbands{i-1}.bl,3),1,1,size(fbands{i-1}.bl,3));
   end
   mask1 = repmat(permute(allstats{i}.mask,[3 1 2]),size(pledata.post,1),1,1);
   effsizesum(i) = effsizefun(post1,bl1,mask1);
   effsizesum_ci(i,:) = bootci(10000,effsizefun,post1,bl1,mask1);
   for ii = 1:(i-1)
       if ii == 1
           post2 = 10.^intdata.post;
           bl2 = repmat(nanmean(10.^intdata.bl,3),1,1,size(intdata.bl,3));
       elseif ii == 2
           post2 = pledata.post;
           bl2 = repmat(nanmean(pledata.bl,3),1,1,size(pledata.bl,3));
       else
           post2 = fbands{ii-1}.post;
           bl2 = repmat(nanmean(fbands{ii-1}.bl,3),1,1,size(fbands{ii-1}.bl,3));
       end
       mask2 = repmat(permute(allstats{ii}.mask,[3 1 2]),size(pledata.post,1),1,1);
       [~,bootstat,S] = bootci_swt(10000,@(p1,b1,m1,p2,b2,m2)effsizedifffun(effsizefun,p1,b1,m1,p2,b2,m2),...
            post1,bl1,mask1,post2,bl2,mask2);
        effsizediff(i,ii) = ibootp(0,bootstat,S);
   end
end



indxmat = abs([intindx pleindx fbandindx(:,2:end)]);
anovap = friedman(indxmat);
for i = 1:size(indxmat,2)
    for ii = 1:i
        mcp(i,ii) = signrank(indxmat(:,i),indxmat(:,ii)); 
    end
end


% panel c : effect size
figure

p = panel('no-manage-font');

set(gcf,'units','normalized','position',[0 0.3 1 0.7])


p.pack('h',{60 40})
p(1).select()
b1 = bar([1 2],effsizesum(1:2),'FaceColor',palecol([1 0 0],0.3));
hold on
b2 = bar([4:(length(effsizesum)+1)],effsizesum(3:end),'FaceColor',palecol([0 0 1],0.3));
set(gca,'XTick',[1 2 4:(length(effsizesum)+1)],'XTickLabel',measnames_short)
e = errorbar([1 2 4:(length(effsizesum)+1)],effsizesum,...
    effsizesum'-effsizesum_ci(:,1),effsizesum_ci(:,2)-effsizesum',...
    'LineStyle','none','LineWidth',2,'Color','k');
%ylabel('Summed effect size')
ylabel('Summed absolute percent change')
FixAxes(gca,14)
set(gca,'XLim',[0 8.75])
%sigstar_frommat(1:size(depeffsize_sum,2),pdif_boot(:,:,i));
p(2).select()
[~,~,cbar]=imagesc_stars(horz(effsizesum)-vert(effsizesum),effsizediff,measnames_short)
FixAxes(gca,14)
cbar.Label.String = 'Summed absolute percent change difference';
cbar.Label.FontSize = 14;

p.marginleft = 26;
Normalize_Clim(gcf,1)
p.marginbottom = 20;
p.marginright = 30;

savefig('Fig1c_effsize_sum_cohend.fig'); export_fig('Fig1c_effsize_sum_cohend.png','-m5')

figure

p = panel('no-manage-font');

set(gcf,'units','normalized','position',[0.15 0.3 0.7 0.7])


p.pack()
p.pack({[0.7 0.7 0.25 0.3]})
p(1).select()
b1 = bar([1 2],mean(indxmat(:,1:2),1),'FaceColor',palecol([1 0 0],0.3));
hold on
b2 = bar([4:(length(effsizesum)+1)],mean(indxmat(:,3:end),1),'FaceColor',palecol([0 0 1],0.3));
set(gca,'XTick',[1 2 4:(length(effsizesum)+1)],'XTickLabel',measnames_short)
%e = errorbar([1 2 4:(length(effsizesum)+1)],effsizesum,...
%    effsizesum'-effsizesum_ci(:,1),effsizesum_ci(:,2)-effsizesum',...
%    'LineStyle','none','LineWidth',2,'Color','k');
ylabel('Summed percent change')
FixAxes(gca,14)
set(gca,'XLim',[0 8.75])
%sigstar_frommat(1:size(depeffsize_sum,2),pdif_boot(:,:,i));
p(2).select()
[~,~,cbar]=imagesc_stars(horz(mean(indxmat,1))-vert(mean(indxmat,1)),mcp,measnames_short)
cbar.Label.String = 'Summed percent change difference difference';
cbar.Label.FontSize = 12;

p.marginleft = 26;
Normalize_Clim(gcf,1)

savefig('Fig1c_prcchange_sum.fig'); export_fig('Fig1c_prcchange_sum.png','-m5')




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



