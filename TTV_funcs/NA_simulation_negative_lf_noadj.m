% Simulation test for nonadditivity

mvals = 2.^[-6:2:6]; % m is mean noise level
nvals = 2.^[-6:2:6]; % n is the variability of noise across trials

settings = struct;
settings.pseudo.prestim = 701:750; 
settings.pseudo.poststim= 901:1300; 
settings.real.prestim = 1301:1350; 
settings.real.poststim = 1501:1900;
settings.prestim_pseudo = settings.pseudo.prestim;
settings.prestim_real = settings.real.prestim;
settings.poststim_pseudo = settings.pseudo.poststim;
settings.poststim_real = settings.real.poststim;
settings.srate = 500;


freqs = {[] [2 4]};

for m = mvals
    for n = nvals
        datacalc = cell(1,48);
        varspont = cell(1,48);
        varevo = cell(1,48);
        spont = cell(1,48);
        evo = cell(1,48);
        corrfact = cell(1,48);
        
        
        parfor i = 1:48
            % Simulate resting-state
            cfg = []; cfg.method = 'mix_oscifrac'; cfg.fsample = 500; cfg.trllen = 4;
            cfg.osci = struct; cfg.frac.ple = rand+0.5; cfg.frac.ampl = 1; cfg.frac.bpfreq = [0.5 50];
            %cfg.osci.s1.ampl = [ones(1,1000) -0.3*sin((1:1000)*pi/(1000))+1];
            cfg.numtrl = 128;
            cfg.noise.ampl = m + n*rand(1,128);
            
            
            %subAmpl = rand+0.5;
            
%             for cc = 1:75 % 5-min resting state
%                 spont = createFN(1.75/2,2000)*subAmpl;
%                 spont = ft_preproc_lowpassfilter(spont,500,10,4);
%                 spont = NormOntoRange(spont,[1 2]);
%                 cfg.osci.s1.ampl = horz(spont); % no task-evoked activity
%                 %cfg.osci.s1.ampl{cc} = horz(spont)-[zeros(1,1000) rand*sin((1:500)*pi/500) zeros(1,500)];
%             end
%             rest{c} = ft_freqsimulation_swt(cfg);
            
            
            % Simulate trial data
            cfg = []; cfg.method = 'mix_oscifrac'; cfg.fsample = 500; cfg.trllen = 5;
            cfg.osci = struct; cfg.frac.ple = rand+0.5; cfg.frac.ampl = 1; cfg.frac.bpfreq = [0.5 50];
            %cfg.osci.s1.ampl = [ones(1,1000) -0.3*sin((1:1000)*pi/(1000))+1];
            cfg.numtrl = 128;
            cfg.noise.ampl = m + n*rand(1,128);
	    cfg.osci.s1.freq = 3;

            
            spont{i} = zeros(128,2500);
            evo{i} = zeros(128,2500);
            for cc = 1:128
                spont{i}(cc,:) = createFN(0.625,2500);%*subAmpl;
                spont{i}(cc,:) = ft_preproc_lowpassfilter(spont{i}(cc,:),500,1,4);
                spont{i}(cc,:) = NormOntoRange(spont{i}(cc,:),[1.25 1.75]);
                %corrfact{i}(cc) = 0.5*(mean(spont{i}(cc,951:1000))-0.5);
                corrfact{i}(cc) = 0.5*(rand+0.625);
                evo{i}(cc,:) = [zeros(1,1500) corrfact{i}(cc)*sin((1:150)*pi/150) zeros(1,850)]
            end
            
            varspont{i} = std(spont{i},[],1);
            varevo{i} = std(evo{i},[],1);
            
            for cc = 1:128
                cfg.osci.s1.ampl{cc} = spont{i}(cc,:)+evo{i}(cc,:);
            end
            
            task{i} = ft_freqsimulation_swt(cfg);
            
            data_allrange = 1:2000;
           % data_allrange = (settings.pseudo.prestim(1)-ceil(settings.srate/5)):(settings.real.poststim(end));
            cfg = []; cfg.method = 'wavelet'; cfg.output = 'fourier'; cfg.foi = exp(linspace(log(2),log(4),6));
            cfg.keeptrials = 'yes'; cfg.toi = task{i}.time{1}(data_allrange); cfg.width = 3;
            freqdata = ft_freqanalysis(cfg,task{i});
            
            timefreq_data = cell(1,length(freqdata.freq)+1);
            timefreq_data{1} = task{i};
            for c = 1:length(timefreq_data{1}.trial)
                timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
            end
            
            for c = 1:length(freqdata.freq)
                for cc = 1:length(task{i}.trial)
                    timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,1,:));
                end
                timefreq_data{c+1}.time = freqdata.time;
                timefreq_data{c+1}.label = task{i}.label;
                for cc = 1:length(freqs)
                    if ~isempty(freqs{cc}) && freqdata.freq(c) >= freqs{cc}(1) && freqdata.freq(c) <= freqs{cc}(2)
                        timefreq_data{c+1}.parent = cc;
                    end
                end
                if ~isfield(timefreq_data{c+1},'parent')
                    if freqs{2}(1) > freqdata.freq(c)
                        timefreq_data{c+1}.parent = 2;
                    elseif freqs{end}(2) < freqdata.freq(c)
                        timefreq_data{c+1}.parent = length(freqs);
                    end
                end
                freqdata.fourierspctrm(:,:,1,:) = []; %remove bits of the matrix each time to save memory
            end
            freqdata = [];
            
           % tmpcfg = []; tmpcfg.bpfilter = 'yes'; tmpcfg.bpfreq = [8 13]; tmpcfg.hilbert = 'complex';
          %  task_bp{i} = ft_preprocessing(tmpcfg,task{i});
            
            newsettings = struct;
            newsettings.units = 'prcchange';
            newsettings.tfparams.fbands = {'Broadband','Delta'};
            newsettings.pseudo = settings.pseudo; newsettings.real = settings.real;
            datacalc_tmp = Calc_sub(newsettings,timefreq_data);
            datacalc{i} = datacalc_tmp{2};
            
            %restmeas_tmp = Rest_calc([8 13],[0.5 50],rest{c});
            %restmeas{c} = restmeas_tmp;
        end
        all_datacalc{find(mvals == m,1),find(nvals == n,1),:} = datacalc;
%        all_restmeas{find(mvals == m,1),find(nvals == n,1),:} = restmeas;
    end
end

%diagnostics - comment out when running for real
% 
% varspont = cat(1,varspont{:});
% varevo = cat(1,varevo{:});
% spont = cat(3,spont{:});
% evo = cat(3,evo{:});
% corrfact = cat(1,corrfact{:});
% 
% for i = 1:48
%     for ii = 1:2000
%         r(i,ii) = corr(spont(:,ii,i),evo(:,ii,i));
%     end
% end
% 
% critvar = -varevo./(2*varspont);
% 
% plot(mean(critvar,1)); hold on; plot(mean(r,1));
% 
% tmp = mergestructs(datacalc);
% 
% figure
% plot(squeeze(mean(tmp.ttversp.real(1,:,:),3)))
% 
% figure
% plot(squeeze(mean(tmp.naddersp.diff(1,:,1,:),4)),'b');
% hold on
% plot(squeeze(mean(tmp.naddersp.diff(1,:,2,:),4)),'r')
% 
% figure 
% hold on
% for i = 1:48
% ft_psdplot(task{i},1,[0.5 50]);
% end
% set(gca,'XScale','log','YScale','log')

%end of diagnostics

for c = 1:length(mvals)
    for cc = 1:length(nvals)
        all_datacalc{c,cc} = mergestructs(all_datacalc{c,cc,:});
%        all_restmeas{c,cc} = mergestructs(all_restmeas{c,cc,:});
    end
end

settings = parload('settings_camcan_1Hz.mat','settings');

datasetinfo = settings.datasetinfo;
datasetinfo.label = datasetinfo.label(1);

allstats_pt = cell(length(mvals),length(nvals));
allstats_ttv = allstats_pt;

prestim_pseudo = 351:400; poststim_pseudo = 401:800; prestim_real = 951:1000; poststim_real = 1001:1400;

        opts2.display_mod = 0;

%opts = rmfield(opts,'parpool');

delete(gcp('nocreate'))
for c = 1:length(mvals)
    opts.minnbchan = 0; opts.nrand = 1000; opts.distmethod = 'distance';
    stats_ersp_pt = cell(1,length(nvals));
    stats_ersp_ttv = cell(1,length(nvals));
    parfor cc = 1:length(nvals)
        stats_ersp_pt{cc} = EasyClusterCorrect({permute(squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,1,:)),[3 2 1]),...
            permute(squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,2,:)),[3 2 1])},datasetinfo,'ft_statfun_fast_signrank',opts);
        stats_ersp_ttv{cc} = EasyClusterCorrect({permute(all_datacalc{c,cc}.ttversp.real(1,:,:),[1 3 2]) 0.*permute(all_datacalc{c,cc}.ttversp.real(1,:,:),[1 3 2])},...
            datasetinfo,'ft_statfun_fast_signrank',opts)
        
        %stats_erp_pt{cc} = EasyClusterCorrect({permute(squeeze(all_datacalc{c,cc}.nadderp.diff(1,:,1,:)),[3 2 1]),...
        %    permute(squeeze(all_datacalc{c,cc}.nadderp.diff(1,:,2,:)),[3 2 1])},datasetinfo,'ft_statfun_fast_signrank',opts);
        %stats_erp_ttv{cc} = EasyClusterCorrect({permute(all_datacalc{c,cc}.ttv.real(1,:,:),[1 3 2]) 0.*permute(all_datacalc{c,cc}.ttversp.real(1,:,:),[1 3 2])},...
        %    datasetinfo,'ft_statfun_fast_signrank',opts)
        % Testing the mediation model
%         med{c,cc} = mediation_covariates(vert(squeeze(mean(all_datacalc{c,cc}.raw.ersp(1,poststim_real,:),2))),...
%             all_restmeas{c,cc}.bp(:,1),vert(squeeze(mean(all_datacalc{c,cc}.raw.ersp(1,prestim_real,:),2))),...
%             vert(squeeze(mean(all_datacalc{c,cc}.raw.ersp(1,poststim_pseudo,:),2))),opts2);

    end
    allstats_ersp_pt(c,:) = stats_ersp_pt;
    allstats_ersp_ttv(c,:) = stats_ersp_ttv;
    %allstats_erp_pt(c,:) = stats_erp_pt;
    %allstats_erp_ttv(c,:) = stats_erp_ttv;
end


save('simulation_negative_lf_alloutputs.mat','allstats_ersp_pt','allstats_ersp_ttv','all_datacalc','-v7.3')

p = panel('no-manage-font');

pos = get(gcf,'position');

set(gcf,'position',[pos(1) pos(2) pos(3)*3 pos(4)*3]);

p.pack(7,7)

p.de.margin = [7 7 7 7];
p.marginleft = 18;

t = linspace(0,800,400);
warning('Hard-coded')
for c = 1:length(mvals)
    for cc = 1:length(nvals)
        p(c,cc).select()
        ylim1 = stdshade(t,squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,1,:)),[0 0 1],0.15,1,'std');
        hold on
        ylim2 = stdshade(t,squeeze(all_datacalc{c,cc}.naddersp.diff(1,:,2,:)),[1 0 0],0.15,1,'std');
        set(gca,'XTickLabel',{},'XLim',[0 800],'YLim',[min(ylim1(1),ylim2(1)) max(ylim1(2),ylim2(2))])
                ytick = get(gca,'YTickLabel');
        if length(ytick) > 4
           set(gca,'YTickLabel',ytick([1:2:length(ytick)])) 
        end
        Plot_sigmask(gca,allstats_ersp_pt{c,cc}.prob < 0.05,'bar')
        set(gca,'XLim',[0 800])
        if c == 1
            title(['Noise TTV ' num2str(round(0.3*nvals(cc),3,'significant'))],'FontSize',10)
        end
        if cc == 1
            AddFigureLabel(gca,['SNR ' newline num2str(round(1/mvals(c),3,'significant'))],'middle_left','FontSize',10)
        
        end
        %         if c == length(mvals)
        %            xlabel('Time (ms)')
        %         end
        %         if cc == 1
        %            ylabel('% change from prestim')
        %         end
    end
end

set(gcf,'Color','w','units','normalized','position',[0 0 1 1])

savefig('Simulation_negative_fig_pseudotrial.fig')
export_fig('Simulation_negative_fig_pseudotrial.png','-m4')



p = panel('no-manage-font');

pos = get(gcf,'position');

set(gcf,'position',[pos(1) pos(2) pos(3)*3 pos(4)*3]);

p.pack(7,7)

p.de.margin = [7 7 7 7];

p.marginleft = 18;

t = linspace(0,800,400);
warning('Hard-coded')
for c = 1:length(mvals)
    for cc = 1:length(nvals)
        p(c,cc).select()
        ylim1 = stdshade(t,squeeze(all_datacalc{c,cc}.ttversp.real(1,:,:)),[0 0 1],0.15,1,'std');
        set(gca,'XTickLabel',{},'XLim',[0 800],'YLim',ylim1)
        ytick = get(gca,'YTickLabel');
        if length(ytick) > 4
           set(gca,'YTickLabel',ytick([1:2:length(ytick)])) 
        end
        Plot_sigmask(gca,allstats_ersp_ttv{c,cc}.prob < 0.05,'bar')
        set(gca,'XLim',[0 800])
        if c == 1
            title(['Noise TTV ' num2str(round(0.3*nvals(cc),3,'significant'))],'FontSize',10)
        end
        if cc == 1
            AddFigureLabel(gca,['SNR ' newline num2str(round(1/mvals(c),3,'significant'))],'middle_left','FontSize',10)
        end
        %         if c == length(mvals)
        %            xlabel('Time (ms)')
        %         end
        %         if cc == 1
        %            ylabel('% change from prestim')
        %         end
    end
end

set(gcf,'Color','w','units','normalized','position',[0 0 1 1])

savefig('Simulation_negative_fig_ttv.fig')
export_fig('Simulation_negative_fig_ttv.png','-m4')

function [restmeas] = Rest_calc(foi,bandpass,sim)
    sim.fsample = 500;
    sim = ft_concat(sim);
    for c = 1:size(sim.trial{1},1)
        restmeas.bp(c) = bandpower(sim.trial{1}(c,:),sim.fsample,foi);
        restmeas.rel_bp(c) = restmeas.bp(c)./bandpower(sim.trial{1}(c,:),sim.fsample,bandpass);
    end
end


function [datacalc] = Calc_sub(settings,sim)

timefreq_data = sim;

numbands = length(timefreq_data);

aucindex = 1:400;

datacalc = cell(1,1);
datacalc{1} = struct;
prestim_pseudo = settings.pseudo.prestim; 
prestim_real = settings.real.prestim;
poststim_pseudo = settings.pseudo.poststim;
poststim_real = settings.real.poststim;
%prestim_pseudo = 351:400; poststim_pseudo = 401:800; prestim_real = 951:1000; poststim_real = 1001:1400;

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
    
    % normalize by 50ms prestim ifor both real and pseudotrial
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


%if ~strcmpi(settings.tfparams.method,'hilbert')
    %timefreq_data = parload(files(1).name,'timefreq_data');
    for c = 2:length(timefreq_data)
        parents(c) = timefreq_data{c}.parent;
    end
    newmeas = cell(1,length(settings.tfparams.fbands));
    for c = 2:length(settings.tfparams.fbands)
        children{c} = find(parents == c);
        fields = fieldnames_recurse(datacalc{c});
        fields = cell_unpack(fields);
        for cc = 1:length(fields)
            tmp = [];
            for ccc = 1:length(children{c})
                dimn = length(size(getfield_nest(datacalc{children{c}(ccc)},fields{cc})));
                tmp = cat(dimn+1,tmp,getfield_nest(datacalc{children{c}(ccc)},(fields{cc})));
            end
            newmeas{c} = assignfield_nest(newmeas{c},fields{cc},nanmean(tmp,dimn+1));
        end
    end
    
    dimn = [];
    children{1} = 2:length(timefreq_data);
    fields = fieldnames_recurse(datacalc{1});
    fields = cell_unpack(fields);
    
    for cc = 1:length(fields)
        tmp = [];
        for ccc = 1:length(children{1})
            dimn = length(size(getfield_nest(datacalc{children{1}(ccc)},fields{cc})));
            tmp = cat(dimn+1,tmp,getfield_nest(datacalc{children{1}(ccc)},(fields{cc})));
        end
        newmeas{1} = assignfield_nest(newmeas{1},fields{cc},nanmean(tmp,dimn+1));
    end
    
    % for broadband, put the ERP and TTV stuff back to the original values
    newmeas{1}.erp = datacalc{1}.erp;
    newmeas{1}.nadderp = datacalc{1}.nadderp;
    newmeas{1}.ttv = datacalc{1}.ttv;
    
    settings.nfreqs = length(settings.tfparams.fbands);
    numbands = settings.nfreqs;
    
    datacalc = newmeas;
    newmeas = [];
%end

end
