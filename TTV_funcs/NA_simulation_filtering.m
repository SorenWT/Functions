% Simulation test for nonadditivity

mvals = 2.^[-6:2:6]; % m is mean noise level
nvals = 2.^[-6:2:6]; % n is the variability of noise across trials
m = mvals(3); n = nvals(3);

datacalc_neg = cell(1,48);
parfor c = 1:48
    % Simulate trial data
    cfg = []; cfg.method = 'mix_oscifrac'; cfg.fsample = 500; cfg.trllen = 4;
    cfg.osci = struct; cfg.frac.ple = rand+0.5; cfg.frac.ampl = 1; cfg.frac.bpfreq = [0.5 50];
    %cfg.osci.s1.ampl = [ones(1,1000) -0.3*sin((1:1000)*pi/(1000))+1];
    cfg.numtrl = 128;
    cfg.noise.ampl = m + n*rand(1,128);
    for cc = 1:128
        spont = createFN(1.75/2,2000);
        spont = ft_preproc_lowpassfilter(spont,500,10,4);
        spont = NormOntoRange(spont,[1 2]);
        cfg.osci.s1.ampl{cc} = horz(spont)-[zeros(1,1000) rand*sin((1:50)*pi/50) zeros(1,950)];
    end
    task{c} = ft_freqsimulation_swt(cfg);
    
    tmpcfg = []; tmpcfg.hilbert = 'complex';
    task{c} = ft_preprocessing(tmpcfg,task{c});
    
    tmpcfg = []; tmpcfg.hpfreq = 1; tmpcfg.hpfilter = 'yes';
    task{c} = ft_preprocessing(tmpcfg,task{c});
    
    settings = struct;
    settings.units = 'prcchange';
    datacalc_tmp = Calc_sub(settings,task{c});
    datacalc_neg{c} = datacalc_tmp{1};

end

task_neg = task;
datacalc_neg = mergestructs(datacalc_neg);

datacalc_pos = cell(1,48);
n = 0; m = 0.5;
parfor c = 1:48
    % Simulate trial data
    cfg = []; cfg.method = 'mix_oscifrac'; cfg.fsample = 500; cfg.trllen = 4;
    cfg.osci = struct; cfg.frac.ple = 1.8; cfg.frac.ampl = 2; cfg.frac.bpfreq = [0.5 50];
    %cfg.osci.s1.ampl = [ones(1,1000) -0.3*sin((1:1000)*pi/(1000))+1];
    cfg.numtrl = 128;
    cfg.noise.ampl = m + n*rand(1,128);
    for cc = 1:128
        spont = createFN(1.75/2,2000);
        spont = ft_preproc_lowpassfilter(spont,500,10,4);
        spont = NormOntoRange(spont,[0.5 0.5]);
        cfg.osci.s1.ampl{cc} = horz(spont);
    end
    task{c} = ft_freqsimulation_swt(cfg);  
    prestim = zeros(1,128);
    poststim = zeros(1,128);
    for cc = 1:128
       prestim(cc) = mean(task{c}.trial{cc}(1,951:1000));
    end
    
    prestim = NormOntoRange(-prestim,[0 3]);
    ttv_spont = std(cat(3,task{c}.trial{:}),[],3);
    ttv_spont = ttv_spont(1,:).^2;
    evo = zeros(128,2000);
    for cc = 1:128
        evo(cc,:) = [zeros(1,1000) prestim(cc)*sin((1:50)*pi/50) zeros(1,950)];
    end
    
    %ttv_evo = std(evo,[],1).^2;
    
    for cc = 1:128
       task{c}.trial{cc}(1,:) = task{c}.trial{cc}(1,:)+evo(cc,:);
    end
    
    tmpcfg = []; tmpcfg.hilbert = 'complex';
    task{c} = ft_preprocessing(tmpcfg,task{c});
    
    settings = struct;
    settings.units = 'prcchange';
    datacalc_tmp = Calc_sub(settings,task{c});
    datacalc_pos{c} = datacalc_tmp{1};

end
task_pos = task;

datacalc_pos = mergestructs(datacalc_pos);

settings = parload('settings_camcan_1Hz.mat','settings');

datasetinfo = settings.datasetinfo;
datasetinfo.label = datasetinfo.label(1);

prestim_pseudo = 351:400; poststim_pseudo = 401:800; prestim_real = 951:1000; poststim_real = 1001:1400;

opts.minnbchan = 0; opts.nrand = 1000; opts.distmethod = 'distance';
stats = cell(1,length(nvals));

stats_erp_pt_neg = EasyClusterCorrect({permute(squeeze(datacalc_neg.nadderp.diff(1,:,1,:)),[3 2 1]),...
    permute(squeeze(datacalc_neg.nadderp.diff(1,:,2,:)),[3 2 1])},datasetinfo,'ft_statfun_fast_signrank',opts);
stats_erp_ttv_neg = EasyClusterCorrect({permute(datacalc_neg.ttv.real(1,:,:),[1 3 2]) 0.*permute(datacalc_neg.ttversp.real(1,:,:),[1 3 2])},...
    datasetinfo,'ft_statfun_fast_signrank',opts);

stats_erp_pt_pos = EasyClusterCorrect({permute(squeeze(datacalc_pos.nadderp.diff(1,:,1,:)),[3 2 1]),...
    permute(squeeze(datacalc_pos.nadderp.diff(1,:,2,:)),[3 2 1])},datasetinfo,'ft_statfun_fast_signrank',opts);
stats_erp_ttv_pos = EasyClusterCorrect({permute(datacalc_pos.ttv.real(1,:,:),[1 3 2]) 0.*permute(datacalc_pos.ttversp.real(1,:,:),[1 3 2])},...
    datasetinfo,'ft_statfun_fast_signrank',opts);


save('simulation_fig2_outputs.mat','datacalc_neg','datacalc_pos','stats_erp_pt_neg','stats_erp_ttv_neg',...
    'stats_erp_pt_pos','stats_erp_ttv_pos','task_neg','task_pos','-v7.3')

%% Making the figure


p = panel('no-manage-font')

p.pack('h',{1/2 1/2})
p(1).pack('v',{1/3 1/3 1/3})
p(2).pack('v',{1/3 1/3 1/3})

p(1,1).select()
t = linspace(0,0.8,400);
plot(t(1:150),[2*ones(1,150)],'k','LineWidth',2)
hold on
plot(t(1:150),[1+0.5*sin(pi*(1:50)/50) ones(1,100)],'k','LineWidth',2)
plot(t(1:150),[sin(pi*(1:50)/50) 0*ones(1,100)],'k','LineWidth',2)
set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none')

p(1,2).select()
t = linspace(0,0.8,400);
plot(t(1:150),mean(datacalc_pos.ttv.real(1,1:150,:),3),'LineWidth',2)
ylabel('Percent change in TTV')
FixAxes(gca,14)

p(1,3).select()
plot(t(1:150),squeeze(mean(datacalc_pos.nadderp.diff(1,1:150,1,:),4)),'b','LineWidth',2)
hold on
plot(t(1:150),squeeze(mean(datacalc_pos.nadderp.diff(1,1:150,2,:),4)),'r','LineWidth',2)
ylabel('Voltage (\muV)')
xlabel('Time (s)')
legend({'Corrected prestim low','Corrected prestim high'})
FixAxes(gca,14)


p(2,1).select()
t = linspace(0,0.8,400);
plot(t(1:150),sin(2*pi*t(1:150)/0.1).*[1-sin(pi*(1:50)/50) ones(1,100)],'k','LineWidth',2)
set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none')

p(2,2).select()
t = linspace(0,0.8,400);
plot(t(1:150),mean(datacalc_neg.ttv.real(1,1:150,:),3),'LineWidth',2)
ylabel('Percent change in TTV')
FixAxes(gca,14)

p(2,3).select()
plot(t(1:150),squeeze(mean(datacalc_neg.nadderp.diff(1,1:150,1,:),4)),'b','LineWidth',2)
hold on
plot(t(1:150),squeeze(mean(datacalc_neg.nadderp.diff(1,1:150,2,:),4)),'r','LineWidth',2)
ylabel('Voltage (\muV)')
xlabel('Time (s)')
legend({'Corrected prestim low','Corrected prestim high'})
FixAxes(gca,14)

p.marginleft = 22;
p.de.marginbottom = 10; 
p.marginbottom = 18;
set(gcf,'color','w')




function [datacalc] = Calc_sub(settings,sim)

timefreq_data{1} = sim;

numbands = 1;

aucindex = 1:400;

datacalc = cell(1,1);
datacalc{1} = struct;

prestim_pseudo = 351:400; poststim_pseudo = 401:800; prestim_real = 951:1000; poststim_real = 1001:1400;

for q = 1:numbands
    nbchan = length(timefreq_data{q}.label);
    timefreq_data{q}.matrix = cat(3,timefreq_data{q}.trial{:});
    
    datacalc{q}.raw.ersp = mean(abs(timefreq_data{q}.matrix),3);
    
    
    %% ERSP NA
    datacat = timefreq_data{q}.matrix;
    
    split_real = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_real,:)),2));
    split_pseudo = squeeze(mean(abs(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
    
    for c = 1:nbchan
        splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
        
        datacalc{q}.naddersp.raw.pseudo(c,:,1) = mean(abs(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.naddersp.raw.pseudo(c,:,2) = mean(abs(datacat(c,:,find(splitindex))),3);
        
        datacalc{q}.naddersp.pseudo(c,:,1) = (mean(abs(datacat(c,poststim_pseudo,find(~splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),3),2));
        datacalc{q}.naddersp.pseudo(c,:,2) = (mean(abs(datacat(c,poststim_pseudo,find(splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_pseudo,find(splitindex))),3),2));
        
        tmp = abs(datacat(c,poststim_pseudo,find(~splitindex)))...
            -mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),2);
        tmppseudo{1} = squeeze(trapz(tmp,2));
        
        tmp = abs(datacat(c,poststim_pseudo,find(splitindex)))...
            -mean(abs(datacat(c,prestim_pseudo,find(splitindex))),2);
        tmppseudo{2} = squeeze(trapz(tmp,2));
        
        switch settings.units
            case 'prcchange'
                datacalc{q}.naddersp.pseudo(c,:,1) = 100*datacalc{q}.naddersp.pseudo(c,:,1)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
                datacalc{q}.naddersp.pseudo(c,:,2) = 100*datacalc{q}.naddersp.pseudo(c,:,2)./mean(mean(abs(datacat(c,prestim_pseudo,:)),3),2);
            case 'zscore'
                datacalc{q}.naddersp.pseudo(c,:,1) = zscore(datacalc{q}.naddersp.pseudo(c,:,1),0,2);
                datacalc{q}.naddersp.pseudo(c,:,2) = zscore(datacalc{q}.naddersp.pseudo(c,:,2),0,2);
            case 'log'
                datacalc{q}.naddersp.pseudo(c,:,1) = 10*log10(datacalc{q}.naddersp.pseudo(c,:,1));
                datacalc{q}.naddersp.pseudo(c,:,2) = 10*log10(datacalc{q}.naddersp.pseudo(c,:,2))
        end
        
        splitindex = split_real(c,:) > median(split_real(c,:));
        
        datacalc{q}.naddersp.raw.real(c,:,1) = mean(abs(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.naddersp.raw.real(c,:,2) = mean(abs(datacat(c,:,find(splitindex))),3);
        
        datacalc{q}.naddersp.real(c,:,1) = (mean(abs(datacat(c,poststim_real,find(~splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_real,find(~splitindex))),3),2));
        datacalc{q}.naddersp.real(c,:,2) = (mean(abs(datacat(c,poststim_real,find(splitindex))),3)...
            -mean(mean(abs(datacat(c,prestim_real,find(splitindex))),3),2));
        
        datacalc{q}.naddersp.diff = datacalc{q}.naddersp.real- datacalc{q}.naddersp.pseudo;
        
        switch settings.units
            case 'prcchange'
                datacalc{q}.naddersp.real(c,:,1) = 100*datacalc{q}.naddersp.real(c,:,1)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
                datacalc{q}.naddersp.real(c,:,2) = 100*datacalc{q}.naddersp.real(c,:,2)./mean(mean(abs(datacat(c,prestim_real,:)),3),2);
            case 'zscore'
                datacalc{q}.naddersp.real(c,:,1) = zscore(datacalc{q}.naddersp.real(c,:,1),0,2);
                datacalc{q}.naddersp.real(c,:,2) = zscore(datacalc{q}.naddersp.real(c,:,2),0,2);
            case 'log'
                datacalc{q}.naddersp.real(c,:,1) = 10*log10(datacalc{q}.naddersp.real(c,:,1));
                datacalc{q}.naddersp.real(c,:,2) = 10*log10(datacalc{q}.naddersp.real(c,:,2));
        end
    end
    %% TTV of ERSP
    
    datacalc{q}.raw.ttversp(:,:) = std(abs(datacat),[],3);
    datacalc{q}.ttversp.pseudo(:,:) = (std(abs(datacat(:,poststim_pseudo,:)),[],3)...
        -mean(std(abs(datacat(:,prestim_pseudo,:)),[],3),2));
    datacalc{q}.ttversp.real(:,:) = (std(abs(datacat(:,poststim_real,:)),[],3)...
        -mean(std(abs(datacat(:,prestim_real,:)),[],3),2));
    
    switch settings.units
        case 'prcchange'
            datacalc{q}.ttversp.real(:,:) = 100*datacalc{q}.ttversp.real(:,:)./mean(datacalc{q}.raw.ttversp(:,prestim_real),2);
            datacalc{q}.ttversp.pseudo(:,:) = 100*datacalc{q}.ttversp.pseudo(:,:)./mean(datacalc{q}.raw.ttversp(:,prestim_pseudo),2);
        case 'zscore'
            datacalc{q}.ttversp.real(:,:) = zscore(datacalc{q}.ttversp.real(:,:),0,2);
            datacalc{q}.ttversp.pseudo(:,:) = zscore(datacalc{q}.ttversp.pseudo(:,:),0,2);
        case 'log'
            datacalc{q}.ttversp.real(:,:) = 10*log10(datacalc{q}.ttversp.real(:,:));
            datacalc{q}.ttversp.pseudo(:,:) = 10*log10(datacalc{q}.ttversp.pseudo(:,:));
    end
    
    %% ERP nonadditivity
    datacalc{q}.erp(:,:) = mean(real(timefreq_data{q}.matrix),3);
    
    split_real = squeeze(mean(real(timefreq_data{q}.matrix(:,prestim_real,:)),2));
    split_pseudo = squeeze(mean(real(timefreq_data{q}.matrix(:,prestim_pseudo,:)),2));
        
    for c = 1:nbchan
        splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
        
        datacalc{q}.nadderp.raw.pseudo(c,:,1) = mean(real(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.nadderp.raw.pseudo(c,:,2) = mean(real(datacat(c,:,find(splitindex))),3);
        
        
        datacalc{q}.nadderp.pseudo(c,:,1) = (mean(real(datacat(c,poststim_pseudo,find(~splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_pseudo,find(~splitindex))),3),2));
        datacalc{q}.nadderp.pseudo(c,:,2) = (mean(real(datacat(c,poststim_pseudo,find(splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_pseudo,find(splitindex))),3),2));
        
        splitindex = split_real(c,:) > median(split_real(c,:));
        
        datacalc{q}.nadderp.raw.real(c,:,1) = mean(real(datacat(c,:,find(~splitindex))),3);
        datacalc{q}.nadderp.raw.real(c,:,2) = mean(real(datacat(c,:,find(splitindex))),3);
        
        
        datacalc{q}.nadderp.real(c,:,1) = (mean(real(datacat(c,poststim_real,find(~splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_real,find(~splitindex))),3),2));
        datacalc{q}.nadderp.real(c,:,2) = (mean(real(datacat(c,poststim_real,find(splitindex))),3)...
            -mean(mean(real(datacat(c,prestim_real,find(splitindex))),3),2));
        
        datacalc{q}.nadderp.diff = datacalc{q}.nadderp.real-datacalc{q}.nadderp.pseudo;
    end
        
    SD_all = zeros(length(timefreq_data{q}.label),size(timefreq_data{q}.matrix,2));
    for c = 1:length(timefreq_data{q}.label)
        SD_all(c,:) = std(real(timefreq_data{q}.matrix(c,:,:)),[],3);
    end
    datacalc{q}.raw.sd(:,:) = SD_all;
    
    % normalize by prestim for both real and pseudotrial
    datacalc{q}.ttv.pseudo(:,:) = (SD_all(:,poststim_pseudo)-...
        mean(SD_all(:,prestim_pseudo),2));
    datacalc{q}.ttv.real(:,:) = (SD_all(:,poststim_real)-...
        mean(SD_all(:,prestim_real),2));
    switch settings.units
        case 'prcchange'
            datacalc{q}.ttv.pseudo(:,:) = 100*datacalc{q}.ttv.pseudo(:,:)./...
                mean(SD_all(:,prestim_pseudo),2);
            datacalc{q}.ttv.real(:,:) = 100*datacalc{q}.ttv.real(:,:)./...
                mean(SD_all(:,prestim_real),2);
        case 'zscore'
            datacalc{q}.ttv.pseudo(:,:) = zscore(datacalc{q}.ttv.pseudo(:,:),0,2);
            datacalc{q}.ttv.real(:,:) = zscore(datacalc{q}.ttv.real(:,:),0,2);
        case 'log'
            datacalc{q}.ttv.pseudo(:,:) = 10*log10(datacalc{q}.ttv.pseudo(:,:));
            datacalc{q}.ttv.real(:,:) = 10*log10(datacalc{q}.ttv.real(:,:));
    end
end
end
