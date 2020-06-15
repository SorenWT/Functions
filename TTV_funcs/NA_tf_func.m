function [settings] = NA_tf_func(settings)

files = dir(settings.files);

pool = gcp('nocreate');

if settings.pool > 1
    if ~isempty(pool) && pool.NumWorkers ~= settings.pool
        delete(pool)
        parpool(settings.pool)
    elseif isempty(pool)
        parpool(settings.pool)
    end
end

foi = cell(1,length(files));

if strcmpi(settings.tfparams.pf_adjust,'yes')
    allfreqs = cell(1,length(files));
    pf = zeros(1,length(files));
end

if settings.pool > 1
    parfor i = 1:length(files)
        if strcmpi(settings.datatype,'EEG')
            EEG = pop_loadset('filename',files(i).name,'filepath',settings.inputdir);
            data = eeglab2fieldtrip(EEG,'preprocessing','none');
            data = ft_struct2single(data);
        else
            data = parload(files(i).name,'data');
        end
        %data = ft_struct2single(data);
        
        if strcmpi(settings.tfparams.pf_adjust,'yes')
            [freqs pf(i)] = NA_convert_alpha_pf(settings,ft_concat(data));
            allfreqs{i} = horz(freqs);
        else
            freqs = settings.tfparams.fbands;
        end
        
        timefreq_data = cell(1,length(freqs));
        if strcmpi(settings.tfparams.continue,'no') ||  ~exist(fullfile(settings.outputdir,[settings.datasetname '_' files(i).name '_calc.mat']),'file')
            switch settings.tfparams.method
                case 'hilbert'
                    if data.fsample ~= settings.srate
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        data = ft_resampledata(cfg,data);
                    end
                    
                    if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                        cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                        data = ft_selectdata(cfg,data);
                    end
                    
                    for q = 1:settings.nfreqs
                        if ~isempty(freqs{q})
                            cfg = [];
                            if isnan(freqs{q}(1))
                                cfg.lpfilter = 'yes'; cfg.lpfreq = freqs{q}(2); cfg.lpinstabilityfix = 'split';
                            else
                                cfg.bpfilter = 'yes'; cfg.bpfreq = freqs{q}; cfg.bpinstabilityfix = 'split';
                            end
                            timefreq_data{q} = ft_preprocessing(cfg,data);
                            %                    timefreq_data{q} = ft_struct2single(timefreq_data{q});
                        else
                            timefreq_data{q} = data;
                            %                    timefreq_data{q} = ft_struct2single(timefreq_data{q});
                        end
                        cfg = []; cfg.hilbert = 'complex';
                        timefreq_data{q} = ft_preprocessing(cfg,timefreq_data{q});
                    end
                    
                case 'wavelet'
                    if settings.srate ~= data.fsample
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        data = ft_resampledata(cfg,data);
                    end
                    
                    %                if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                    %                    cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                    %                    data = ft_selectdata(cfg,data);
                    %                end
                    
                    data_allrange = (settings.pseudo.prestim(1)-ceil(settings.srate/5)):(settings.real.poststim(end));
                    cfg = []; cfg.method = 'wavelet'; cfg.output = 'fourier'; cfg.foi = exp(linspace(log(freqs{2}(1)),log(freqs{end}(2)),50));
                    cfg.keeptrials = 'yes'; cfg.toi = data.time{1}(data_allrange); cfg.width = 3;
                    freqdata = ft_freqanalysis(cfg,data);
                    foi{i} = freqdata.freq;
                    
                    timefreq_data{1} = data;
                    for c = 1:length(timefreq_data{1}.trial)
                        timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
                    end
                    
                    for c = 1:length(foi{i})
                        for cc = 1:length(data.trial)
                            timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,1,:));
                        end
                        timefreq_data{c+1}.time = freqdata.time;
                        timefreq_data{c+1}.label = data.label;
                        for cc = 1:length(freqs)
                            if ~isempty(freqs{cc}) && foi{i}(c) >= freqs{cc}(1) && foi{i}(c) <= freqs{cc}(2)
                                timefreq_data{c+1}.parent = cc;
                            end
                        end
                        if ~isfield(timefreq_data{c+1},'parent')
                            if freqs{2}(1) > foi{i}(c)
                                timefreq_data{c+1}.parent = 2;
                            elseif freqs{end}(2) < foi{i}(c)
                                timefreq_data{c+1}.parent = length(freqs);
                            end
                        end
                        freqdata.fourierspctrm(:,:,1,:) = []; %remove bits of the matrix each time to save memory
                    end
                    freqdata = [];
                case 'fft'
                    cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                    data = ft_resampledata(cfg,data);
                    
                    if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                        cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                        data = ft_selectdata(cfg,data);
                    end
                    
                    data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
                    cfg = []; cfg.method = 'mtmconvol'; cfg.output = 'fourier'; cfg.foi = exp(linspace(log(freqs{2}(1)),log(freqs{end}(2)),50));
                    cfg.keeptrials = 'yes'; cfg.taper = 'hanning'; cfg.t_ftimwin = ones(length(cfg.foi))*1.5;
                    cfg.toi = [data.time{1}(data_allrange)];
                    freqdata = ft_freqanalysis(cfg,data);
                    foi{i} = freqdata.freq;
                    
                    timefreq_data{1} = data;
                    for c = 1:length(timefreq_data{1}.trial)
                        timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
                    end
                    
                    for c = 1:length(foi{i})
                        for cc = 1:length(data.trial)
                            timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,1,:));
                        end
                        timefreq_data{c+1}.time = freqdata.time;
                        timefreq_data{c+1}.label = data.label;
                        for cc = 1:length(freqs)
                            if ~isempty(freqs{cc}) && foi{i}(c) >= freqs{cc}(1) && foi{i}(c) <= freqs{cc}(2)
                                timefreq_data{c+1}.parent = cc;
                            end
                        end
                        if ~isfield(timefreq_data{c+1},'parent')
                            if freqs{2}(1) > foi{i}(c)
                                timefreq_data{c+1}.parent = 2;
                            elseif freqs{end}(2) < foi{i}(c)
                                timefreq_data{c+1}.parent = length(freqs);
                            end
                        end
                        freqdata.fourierspctrm(:,:,1,:) = [];
                    end
                    freqdata = [];
                    
                case 'irasa'
                    if ~exist([files(i).name '_IRASAtf.mat'],'file')
                        if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                            cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                            data = ft_selectdata(cfg,data);
                        end
                        
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        newdata = ft_resampledata(cfg,data);
                        
                        data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
                        
                        cfg = []; cfg.oscifrac = settings.tfparams.oscifrac; cfg.winsize = 2;
                        cfg.toi = newdata.time{1}(data_allrange); cfg.foi = freqs(2:end);
                        freqdata = IRASA_tf(cfg,data);
                    else
                        freqdata = parload([files(i).name '_IRASAtf.mat']',settings.tfparams.oscifrac);
                        
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        newdata = ft_resampledata(cfg,data);

                        cfg = []; cfg.latency = [min(freqdata.time) max(freqdata.time)];
                        newdata = ft_selectdata(cfg,newdata);
                        
                        %data_allrange = round((settings.pseudo.prestim(1)-settings.srate/5)):(settings.real.poststim(end));
                    end
                    
                    
                    
                    foi{i} = freqdata.freq;
                    cfg.foi = freqdata.freq;
                    
                    freqdata.fourierspctrm = permute(freqdata.fourierspctrm,[1 3 2 4]);
                    
                    if strcmpi(settings.tfparams.naddple,'yes')
                        Calc_naddPLE(settings,freqdata,files(i).name)
                    end
                    
                   % freqdata.fourierspctrm = freqdata.fourierspctrm + 0.001*min(min(min(min(freqdata.fourierspctrm))))*j; % add a tiny imaginary component for compatibility
                    
                    timefreq_data{1} = newdata;
                    for c = 1:length(timefreq_data{1}.trial)
                        timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c};
                    end
                    
                    for c = 1:length(foi{i})
                        for cc = 1:length(data.trial)
                            timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,c,:));
                        end
                        timefreq_data{c+1}.time = freqdata.time;
                        timefreq_data{c+1}.label = data.label;
                        
                        for cc = 1:length(freqs)
                            if ~isempty(freqs{cc}) && foi{i}(c) >= freqs{cc}(1) && foi{i}(c) <= freqs{cc}(2)
                                timefreq_data{c+1}.parent = cc;
                            end
                        end
                        
                        if ~isfield(timefreq_data{c+1},'parent')
                            timefreq_data{c+1}.parent = length(settings.tfparams.fbandnames)+1; % don't assign frequencies not in a specified band
                        end
                    end
                    freqdata = [];
                    
            end
            
            
            for q = 1:length(freqs)
                if ~isempty(freqs{q}) && isnan(freqs{q}(1)) && isfield(settings,'rest')
                    freqs{q}(1) = settings.rest.bandpass(1);
                end
            end
            
            Calc_sub(settings,timefreq_data,files(i).name)
        end
        %    parsave([settings.outputdir '/' files(i).name '_timefreq_filtered.mat'],'timefreq_data',timefreq_data);
    end
else
    for i = 1:length(files)
        if strcmpi(settings.datatype,'EEG')
            EEG = pop_loadset('filename',files(i).name,'filepath',settings.inputdir);
            data = eeglab2fieldtrip(EEG,'preprocessing','none');
            data = ft_struct2single(data);
        else
            data = parload(files(i).name,'data');
        end
        %data = ft_struct2single(data);
        
        if strcmpi(settings.tfparams.pf_adjust,'yes')
            [freqs pf(i)] = NA_convert_alpha_pf(settings,ft_concat(data));
            allfreqs{i} = horz(freqs);
        else
            freqs = settings.tfparams.fbands;
        end
        
        timefreq_data = cell(1,length(freqs));
        if strcmpi(settings.tfparams.continue,'no') ||  ~exist(fullfile(settings.outputdir,[settings.datasetname '_' files(i).name '_calc.mat']),'file')
            switch settings.tfparams.method
                case 'hilbert'
                    if data.fsample ~= settings.srate
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        data = ft_resampledata(cfg,data);
                    end
                    
                    if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                        cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                        data = ft_selectdata(cfg,data);
                    end
                    
                    for q = 1:settings.nfreqs
                        if ~isempty(freqs{q})
                            cfg = [];
                            if isnan(freqs{q}(1))
                                cfg.lpfilter = 'yes'; cfg.lpfreq = freqs{q}(2); cfg.lpinstabilityfix = 'split';
                            else
                                cfg.bpfilter = 'yes'; cfg.bpfreq = freqs{q}; cfg.bpinstabilityfix = 'split';
                            end
                            timefreq_data{q} = ft_preprocessing(cfg,data);
                            %                    timefreq_data{q} = ft_struct2single(timefreq_data{q});
                        else
                            timefreq_data{q} = data;
                            %                    timefreq_data{q} = ft_struct2single(timefreq_data{q});
                        end
                        cfg = []; cfg.hilbert = 'complex';
                        timefreq_data{q} = ft_preprocessing(cfg,timefreq_data{q});
                    end
                    
                case 'wavelet'
                    if settings.srate ~= data.fsample
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        data = ft_resampledata(cfg,data);
                    end
                    
                    %                if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                    %                    cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                    %                    data = ft_selectdata(cfg,data);
                    %                end
                    
                    data_allrange = (settings.pseudo.prestim(1)-ceil(settings.srate/5)):(settings.real.poststim(end));
                    cfg = []; cfg.method = 'wavelet'; cfg.output = 'fourier'; cfg.foi = exp(linspace(log(freqs{2}(1)),log(freqs{end}(2)),50));
                    cfg.keeptrials = 'yes'; cfg.toi = data.time{1}(data_allrange); cfg.width = 3;
                    freqdata = ft_freqanalysis(cfg,data);
                    foi{i} = freqdata.freq;
                    
                    timefreq_data{1} = data;
                    for c = 1:length(timefreq_data{1}.trial)
                        timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
                    end
                    
                    for c = 1:length(foi{i})
                        for cc = 1:length(data.trial)
                            timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,1,:));
                        end
                        timefreq_data{c+1}.time = freqdata.time;
                        timefreq_data{c+1}.label = data.label;
                        for cc = 1:length(freqs)
                            if ~isempty(freqs{cc}) && foi{i}(c) >= freqs{cc}(1) && foi{i}(c) <= freqs{cc}(2)
                                timefreq_data{c+1}.parent = cc;
                            end
                        end
                        if ~isfield(timefreq_data{c+1},'parent')
                            if freqs{2}(1) > foi{i}(c)
                                timefreq_data{c+1}.parent = 2;
                            elseif freqs{end}(2) < foi{i}(c)
                                timefreq_data{c+1}.parent = length(freqs);
                            end
                        end
                        freqdata.fourierspctrm(:,:,1,:) = []; %remove bits of the matrix each time to save memory
                    end
                    freqdata = [];
                case 'fft'
                    cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                    data = ft_resampledata(cfg,data);
                    
                    if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                        cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                        data = ft_selectdata(cfg,data);
                    end
                    
                    data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
                    cfg = []; cfg.method = 'mtmconvol'; cfg.output = 'fourier'; cfg.foi = exp(linspace(log(freqs{2}(1)),log(freqs{end}(2)),50));
                    cfg.keeptrials = 'yes'; cfg.taper = 'hanning'; cfg.t_ftimwin = ones(length(cfg.foi))*1.5;
                    cfg.toi = [data.time{1}(data_allrange)];
                    freqdata = ft_freqanalysis(cfg,data);
                    foi{i} = freqdata.freq;
                    
                    timefreq_data{1} = data;
                    for c = 1:length(timefreq_data{1}.trial)
                        timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c}(:,data_allrange);
                    end
                    
                    for c = 1:length(foi{i})
                        for cc = 1:length(data.trial)
                            timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,1,:));
                        end
                        timefreq_data{c+1}.time = freqdata.time;
                        timefreq_data{c+1}.label = data.label;
                        for cc = 1:length(freqs)
                            if ~isempty(freqs{cc}) && foi{i}(c) >= freqs{cc}(1) && foi{i}(c) <= freqs{cc}(2)
                                timefreq_data{c+1}.parent = cc;
                            end
                        end
                        if ~isfield(timefreq_data{c+1},'parent')
                            if freqs{2}(1) > foi{i}(c)
                                timefreq_data{c+1}.parent = 2;
                            elseif freqs{end}(2) < foi{i}(c)
                                timefreq_data{c+1}.parent = length(freqs);
                            end
                        end
                        freqdata.fourierspctrm(:,:,1,:) = [];
                    end
                    freqdata = [];
                    
                case 'irasa'
                    if ~exist([files(i).name '_IRASAtf.mat'],'file')
                        if isfield(settings.tfparams,'condition') && ~strcmpi(settings.tfparams.condition,'all')
                            cfg = []; cfg.trials = find(ismember(data.trialinfo(:,1),settings.tfparams.condition));
                            data = ft_selectdata(cfg,data);
                        end
                        
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        newdata = ft_resampledata(cfg,data);
                        
                        data_allrange = (settings.pseudo.prestim(1)-settings.srate/5):(settings.real.poststim(end));
                        
                        cfg = []; cfg.oscifrac = settings.tfparams.oscifrac; cfg.winsize = 2;
                        cfg.toi = newdata.time{1}(data_allrange); cfg.foi = freqs(2:end);
                        freqdata = IRASA_tf(cfg,data);
                    else
                        freqdata = parload([files(i).name '_IRASAtf.mat']',settings.tfparams.oscifrac);
                        
                        cfg = []; cfg.resamplefs = settings.srate; cfg.detrend = 'no';
                        newdata = ft_resampledata(cfg,data);
                        
                        cfg = []; cfg.latency = [min(freqdata.time) max(freqdata.time)];
                        newdata = ft_selectdata(cfg,newdata);
                        
                        %data_allrange = round((settings.pseudo.prestim(1)-settings.srate/5)):(settings.real.poststim(end));
                    end
                    
                    
                    
                    foi{i} = freqdata.freq;
                    cfg.foi = freqdata.freq;
                    
                    freqdata.fourierspctrm = permute(freqdata.fourierspctrm,[1 3 2 4]);
                    
                    if strcmpi(settings.tfparams.naddple,'yes')
                        Calc_naddPLE(settings,freqdata,files(i).name)
                    end
                    
                   % freqdata.fourierspctrm = freqdata.fourierspctrm + 0.001*min(min(min(min(freqdata.fourierspctrm))))*j; % add a tiny imaginary component for compatibility
                    
                    timefreq_data{1} = data;
                    for c = 1:length(timefreq_data{1}.trial)
                        timefreq_data{1}.trial{c} = timefreq_data{1}.trial{c};
                    end
                    
                    for c = 1:length(foi{i})
                        for cc = 1:length(data.trial)
                            timefreq_data{c+1}.trial{cc} = squeeze(freqdata.fourierspctrm(cc,:,c,:));
                        end
                        timefreq_data{c+1}.time = freqdata.time;
                        timefreq_data{c+1}.label = data.label;
                        
                        for cc = 1:length(freqs)
                            if ~isempty(freqs{cc}) && foi{i}(c) >= freqs{cc}(1) && foi{i}(c) <= freqs{cc}(2)
                                timefreq_data{c+1}.parent = cc;
                            end
                        end
                        
                        if ~isfield(timefreq_data{c+1},'parent')
                            timefreq_data{c+1}.parent = length(settings.tfparams.fbandnames)+1; % don't assign frequencies not in a specified band
                        end
                    end
                    freqdata = [];
                    
            end
            
            
            for q = 1:length(freqs)
                if ~isempty(freqs{q}) && isnan(freqs{q}(1)) && isfield(settings,'rest')
                    freqs{q}(1) = settings.rest.bandpass(1);
                end
            end
            
            Calc_sub(settings,timefreq_data,files(i).name)
        end
        %    parsave([settings.outputdir '/' files(i).name '_timefreq_filtered.mat'],'timefreq_data',timefreq_data);
    end
end

if strcmpi(settings.tfparams.pf_adjust,'yes')
    allfreqs = cat(1,allfreqs{:});
    settings.tfparams.fbands = allfreqs;
    settings.alpha_pf = pf;
end

if ~strcmpi(settings.tfparams.method,'hilbert')
    settings.nfreqs = length(foi{1})+1;
    
    prestim_pseudo = settings.pseudo.prestim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    prestim_real = settings.real.prestim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    poststim_pseudo = settings.pseudo.poststim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    poststim_real = settings.real.poststim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    
    settings.pseudo.prestim = prestim_pseudo;
    settings.real.prestim = prestim_real;
    settings.pseudo.poststim = poststim_pseudo;
    settings.real.poststim = poststim_real;
end


end

function [settings] = Calc_sub(settings,timefreq_data,filename)

numbands = length(timefreq_data);

aucindex = settings.aucindex;

%if isempty(gcp('nocreate'))
%    parpool(numbands)
%end

datacalc = cell(1,length(settings.tfparams.fbands));
for c = 1:length(datacalc)
    datacalc{c} = struct;
end

if strcmpi(settings.tfparams.parameter,'power')
    if strcmpi(settings.tfparams.method,'irasa')
        pfunc = @doNothing;
    else
        pfunc = @sq;
    end
elseif strcmpi(settings.tfparams.parameter,'amplitude')
    if strcmpi(settings.tfparams.method,'irasa')
        pfunc = @sqrt;
    else
        pfunc = @doNothing;
    end
end

if strcmpi(settings.tfparams.method,'hilbert') || strcmpi(settings.tfparams.method,'irasa')
    prestim_pseudo = settings.pseudo.prestim;
    prestim_real = settings.real.prestim;
    poststim_pseudo = settings.pseudo.poststim;
    poststim_real = settings.real.poststim;
else
    prestim_pseudo = settings.pseudo.prestim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    prestim_real = settings.real.prestim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    poststim_pseudo = settings.pseudo.poststim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    poststim_real = settings.real.poststim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
end

datacat = cell(size(timefreq_data));
for c = 1:length(datacat)
    datacat{c} = struct;
end

nbchan = length(timefreq_data{1}.label);

for q = 1:numbands
    trlmat = cat(3,timefreq_data{q}.trial{:});
       
    if strcmpi(settings.tfparams.method,'irasa')
       trlmat(find(trlmat <= 0)) = NaN; 
    end
    
    datacat{q}.pow = pfunc(abs(trlmat));
    datacat{q}.real = real(trlmat);
    
    datacat{q}.raw.itc = abs(nanmean(trlmat./abs(trlmat),3));

    datacat{q}.itc.real(:,:) = (datacat{q}.raw.itc(:,poststim_real)-nanmean(datacat{q}.raw.itc(:,prestim_real),2));
    datacat{q}.itc.pseudo(:,:) = (datacat{q}.raw.itc(:,poststim_pseudo)-nanmean(datacat{q}.raw.itc(:,prestim_pseudo),2))./...
        nanmean(datacat{q}.raw.itc(:,prestim_pseudo),2);
    
    switch settings.units
        case 'prcchange'
            datacat{q}.itc.real(:,:) = 100*datacat{q}.itc.real(:,:)./nanmean(datacat{q}.raw.itc(:,prestim_real),2);
            datacat{q}.itc.pseudo(:,:) = 100*datacat{q}.itc.pseudo(:,:)./nanmean(datacat{q}.raw.itc(:,prestim_pseudo),2);
        case 'zscore'
            datacat{q}.itc.real(:,:) = zscore(datacat{q}.itc.real(:,:),0,2);
            datacat{q}.itc.pseudo(:,:) = zscore(datacat{q}.itc.pseudo(:,:),0,2);
        case 'log'
            datacat{q}.itc.real(:,:) = 10*log10(datacat{q}.itc.real(:,:));
            datacat{q}.itc.pseudo(:,:) = 10*log10(datacat{q}.itc.pseudo(:,:));
    end
end

if ~strcmpi(settings.tfparams.method,'hilbert')
    %timefreq_data = parload(files(1).name,'timefreq_data');
    for c = 2:length(timefreq_data)
        parents(c) = timefreq_data{c}.parent;
    end
    newmeas = cell(1,length(settings.tfparams.fbands));
    
    for c = 2:length(settings.tfparams.fbands)
        children{c} = find(parents == c);
        fields = fieldnames_recurse(datacat{c});
        fields = cell_unpack(fields);
        for cc = 1:length(fields)
            tmp = [];
            for ccc = 1:length(children{c})
                dimn = length(size(getfield_nest(datacat{children{c}(ccc)},fields{cc})));
                tmp = cat(dimn+1,tmp,getfield_nest(datacat{children{c}(ccc)},(fields{cc})));
            end
            newmeas{c} = assignfield_nest(newmeas{c},fields{cc},nanmean(tmp,dimn+1));
        end
    end
    
    dimn = [];
    children{1} = 2:length(timefreq_data);
    fields = fieldnames_recurse(datacat{1});
    fields = cell_unpack(fields);
    
    for cc = 1:length(fields)
        tmp = [];
        for ccc = 1:length(children{1})
            dimn = length(size(getfield_nest(datacat{children{1}(ccc)},fields{cc})));
            tmp = cat(dimn+1,tmp,getfield_nest(datacat{children{1}(ccc)},(fields{cc})));
        end
        newmeas{1} = assignfield_nest(newmeas{1},fields{cc},nanmean(tmp,dimn+1));
    end
    
    % for broadband, put the ERP and TTV stuff back to the original values
    newmeas{1}.real = datacat{1}.real;
    
    settings.nfreqs = length(settings.tfparams.fbands);
    numbands = settings.nfreqs;
    
    datacat = newmeas;
    newmeas = [];
end

for q = 1:settings.nfreqs
    SD_all = zeros(nbchan,size(datacat{q}.real,2));
    for c = 1:length(timefreq_data{q}.label)
        SD_all(c,:) = nanstd(datacat{q}.real(c,:,:),[],3);
    end
    
    datacalc{q}.raw.sd = SD_all;
    
    % normalize by prestim for both real and pseudotrial
    datacalc{q}.ttv.pseudo(:,:) = (SD_all(:,poststim_pseudo)-...
        mean(SD_all(:,prestim_pseudo),2));
    datacalc{q}.ttv.real(:,:) = (SD_all(:,poststim_real)-...
        mean(SD_all(:,prestim_real),2));
    switch settings.units
        case 'prcchange'
            datacalc{q}.ttv.pseudo(:,:) = 100*datacalc{q}.ttv.pseudo(:,:)./...
                nanmean(SD_all(:,prestim_pseudo),2);
            datacalc{q}.ttv.real(:,:) = 100*datacalc{q}.ttv.real(:,:)./...
                nanmean(SD_all(:,prestim_real),2);
        case 'zscore'
            datacalc{q}.ttv.pseudo(:,:) = zscore(datacalc{q}.ttv.pseudo(:,:),0,2);
            datacalc{q}.ttv.real(:,:) = zscore(datacalc{q}.ttv.real(:,:),0,2);
        case 'log'
            datacalc{q}.ttv.pseudo(:,:) = 10*log10(datacalc{q}.ttv.pseudo(:,:));
            datacalc{q}.ttv.real(:,:) = 10*log10(datacalc{q}.ttv.real(:,:));
    end
    
    
    %% ERSP and ITC    
    datacalc{q}.raw.ersp(:,:) = nanmean(datacat{q}.pow,3);
    
    datacalc{q}.ersp.pseudo(:,:) = (nanmean(datacat{q}.pow(:,poststim_pseudo,:),3)...
        -nanmean(nanmean(datacat{q}.pow(:,prestim_pseudo,:),3),2));
    datacalc{q}.ersp.real(:,:) = (nanmean(datacat{q}.pow(:,poststim_real,:),3)...
        -nanmean(nanmean(datacat{q}.pow(:,prestim_real,:),3),2));
    
    switch settings.units
        case 'prcchange'
            datacalc{q}.ersp.pseudo(:,:) = 100*datacalc{q}.ersp.pseudo(:,:)./nanmean(nanmean(datacat{q}.pow(:,prestim_pseudo,:),3),2);
            datacalc{q}.ersp.real(:,:) = 100*datacalc{q}.ersp.real(:,:)./nanmean(nanmean(datacat{q}.pow(:,prestim_real,:),3),2);
        case 'zscore'
            datacalc{q}.ersp.pseudo(:,:) = zscore(datacalc{q}.ersp.pseudo(:,:),0,2);
            datacalc{q}.ersp.real(:,:) = zscore(datacalc{q}.ersp.real(:,:),0,2);
        case 'log'
            datacalc{q}.ersp.pseudo(:,:) = 10*log10(datacalc{q}.ersp.pseudo(:,:));
            datacalc{q}.ersp.real(:,:) = 10*log10(datacalc{q}.ersp.real(:,:));
    end
    
    split_real = squeeze(nanmean(datacat{q}.pow(:,prestim_real,:),2));
    split_pseudo = squeeze(nanmean(datacat{q}.pow(:,prestim_pseudo,:),2));
    
    for c = 1:nbchan
        splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
        
        datacalc{q}.naddersp.raw.pseudo(c,:,1) = nanmean(datacat{q}.pow(c,:,find(~splitindex)),3);
        datacalc{q}.naddersp.raw.pseudo(c,:,2) = nanmean(datacat{q}.pow(c,:,find(splitindex)),3);
        
        datacalc{q}.naddersp.pseudo(c,:,1) = (nanmean(datacat{q}.pow(c,poststim_pseudo,find(~splitindex)),3)...
            -nanmean(nanmean(datacat{q}.pow(c,prestim_pseudo,find(~splitindex)),3),2));
        datacalc{q}.naddersp.pseudo(c,:,2) = (nanmean(datacat{q}.pow(c,poststim_pseudo,find(splitindex)),3)...
            -nanmean(nanmean(datacat{q}.pow(c,prestim_pseudo,find(splitindex)),3),2));
        
        switch settings.units
            case 'prcchange'
                datacalc{q}.naddersp.pseudo(c,:,1) = 100*datacalc{q}.naddersp.pseudo(c,:,1)./nanmean(nanmean(datacat{q}.pow(c,prestim_pseudo,:),3),2);
                datacalc{q}.naddersp.pseudo(c,:,2) = 100*datacalc{q}.naddersp.pseudo(c,:,2)./nanmean(nanmean(datacat{q}.pow(c,prestim_pseudo,:),3),2);
            case 'zscore'
                datacalc{q}.naddersp.pseudo(c,:,1) = zscore(datacalc{q}.naddersp.pseudo(c,:,1),0,2);
                datacalc{q}.naddersp.pseudo(c,:,2) = zscore(datacalc{q}.naddersp.pseudo(c,:,2),0,2);
            case 'log'
                datacalc{q}.naddersp.pseudo(c,:,1) = 10*log10(datacalc{q}.naddersp.pseudo(c,:,1));
                datacalc{q}.naddersp.pseudo(c,:,2) = 10*log10(datacalc{q}.naddersp.pseudo(c,:,2));
        end
        
        splitindex = split_real(c,:) > median(split_real(c,:));
        
        datacalc{q}.naddersp.raw.real(c,:,1) = nanmean(datacat{q}.pow(c,:,find(~splitindex)),3);
        datacalc{q}.naddersp.raw.real(c,:,2) = nanmean(datacat{q}.pow(c,:,find(splitindex)),3);
        
        datacalc{q}.naddersp.real(c,:,1) = (nanmean(datacat{q}.pow(c,poststim_real,find(~splitindex)),3)...
            -nanmean(nanmean(datacat{q}.pow(c,prestim_real,find(~splitindex)),3),2));
        datacalc{q}.naddersp.real(c,:,2) = (nanmean(datacat{q}.pow(c,poststim_real,find(splitindex)),3)...
            -nanmean(nanmean(datacat{q}.pow(c,prestim_real,find(splitindex)),3),2));
       
        switch settings.units
            case 'prcchange'
                datacalc{q}.naddersp.real(c,:,1) = 100*datacalc{q}.naddersp.real(c,:,1)./nanmean(nanmean(datacat{q}.pow(c,prestim_real,:),3),2);
                datacalc{q}.naddersp.real(c,:,2) = 100*datacalc{q}.naddersp.real(c,:,2)./nanmean(nanmean(datacat{q}.pow(c,prestim_real,:),3),2);
            case 'zscore'
                datacalc{q}.naddersp.real(c,:,1) = zscore(datacalc{q}.naddersp.real(c,:,1),0,2);
                datacalc{q}.naddersp.real(c,:,2) = zscore(datacalc{q}.naddersp.real(c,:,2),0,2);
            case 'log'
                datacalc{q}.naddersp.real(c,:,1) = 10*log10(datacalc{q}.naddersp.real(c,:,1));
                datacalc{q}.naddersp.real(c,:,2) = 10*log10(datacalc{q}.naddersp.real(c,:,2));
        end
        
        datacalc{q}.naddersp.diff = datacalc{q}.naddersp.real - datacalc{q}.naddersp.pseudo;
    end
    %% TTV of ERSP
    
    datacalc{q}.raw.ttversp(:,:) = nanstd(datacat{q}.pow,[],3);
    datacalc{q}.ttversp.pseudo(:,:) = (nanstd(datacat{q}.pow(:,poststim_pseudo,:),[],3)...
        -nanmean(nanstd(datacat{q}.pow(:,prestim_pseudo,:),[],3),2));
    datacalc{q}.ttversp.real(:,:) = (std(datacat{q}.pow(:,poststim_real,:),[],3)...
        -nanmean(nanstd(datacat{q}.pow(:,prestim_real,:),[],3),2));
    
    switch settings.units
        case 'prcchange'
            datacalc{q}.ttversp.real(:,:) = 100*datacalc{q}.ttversp.real(:,:)./nanmean(datacalc{q}.raw.ttversp(:,prestim_real),2);
            datacalc{q}.ttversp.pseudo(:,:) = 100*datacalc{q}.ttversp.pseudo(:,:)./nanmean(datacalc{q}.raw.ttversp(:,prestim_pseudo),2);
        case 'zscore'
            datacalc{q}.ttversp.real(:,:) = zscore(datacalc{q}.ttversp.real(:,:),0,2);
            datacalc{q}.ttversp.pseudo(:,:) = zscore(datacalc{q}.ttversp.pseudo(:,:),0,2);
        case 'log'
            datacalc{q}.ttversp.real(:,:) = 10*log10(datacalc{q}.ttversp.real(:,:));
            datacalc{q}.ttversp.pseudo(:,:) = 10*log10(datacalc{q}.ttversp.pseudo(:,:));
    end
    
    %% ERP nonadditivity
    datacalc{q}.erp(:,:) = nanmean(datacat{q}.real,3);
    
    split_real = squeeze(nanmean(datacat{q}.real(:,prestim_real,:),2));
    split_pseudo = squeeze(nanmean(datacat{q}.real(:,prestim_pseudo,:),2));
    
    
    for c = 1:nbchan
        splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
        
        datacalc{q}.nadderp.raw.pseudo(c,:,1) = nanmean(datacat{q}.real(c,:,find(~splitindex)),3);
        datacalc{q}.nadderp.raw.pseudo(c,:,2) = nanmean(datacat{q}.real(c,:,find(splitindex)),3);
        
        
        datacalc{q}.nadderp.pseudo(c,:,1) = (nanmean(datacat{q}.real(c,poststim_pseudo,find(~splitindex)),3)...
            -nanmean(nanmean(datacat{q}.real(c,prestim_pseudo,find(~splitindex)),3),2));
        datacalc{q}.nadderp.pseudo(c,:,2) = (nanmean(datacat{q}.real(c,poststim_pseudo,find(splitindex)),3)...
            -nanmean(nanmean(datacat{q}.real(c,prestim_pseudo,find(splitindex)),3),2));
        
        splitindex = split_real(c,:) > median(split_real(c,:));
        
        datacalc{q}.nadderp.raw.real(c,:,1) = nanmean(datacat{q}.real(c,:,find(~splitindex)),3);
        datacalc{q}.nadderp.raw.real(c,:,2) = nanmean(datacat{q}.real(c,:,find(splitindex)),3);
        
        
        datacalc{q}.nadderp.real(c,:,1) = (nanmean(datacat{q}.real(c,poststim_real,find(~splitindex)),3)...
            -nanmean(nanmean(datacat{q}.real(c,prestim_real,find(~splitindex)),3),2));
        datacalc{q}.nadderp.real(c,:,2) = (nanmean(datacat{q}.real(c,poststim_real,find(splitindex)),3)...
            -nanmean(nanmean(datacat{q}.real(c,prestim_real,find(splitindex)),3),2));
        
        datacalc{q}.nadderp.diff = datacalc{q}.nadderp.real-datacalc{q}.nadderp.pseudo;
    end
    
    datacalc{q}.itc = datacat{q}.itc;
    datacalc{q}.raw.itc = datacat{q}.raw.itc;
end


%% Calculating indices
for q = 1:numbands
    try
        datacalc{q}.ttvindex = squeeze(trapz(datacalc{q}.ttv.real(:,aucindex,:),2));% - squeeze(trapz(datacalc{q}.ttv.pseudo(:,aucindex,:),2));
        datacalc{q}.erspindex = squeeze(trapz(datacalc{q}.ersp.real(:,aucindex,:),2));% - squeeze(trapz(datacalc{q}.ersp.pseudo(:,aucindex,:),2));
        datacalc{q}.itcindex = squeeze(trapz(datacalc{q}.itc.real(:,aucindex,:),2));% - squeeze(trapz(datacalc{q}.itc.pseudo(:,aucindex,:),2));
        datacalc{q}.ttverspindex = squeeze(trapz(datacalc{q}.ttversp.real(:,aucindex,:),2));% - squeeze(trapz(datacalc{q}.ttversp.pseudo(:,aucindex,:),2));
    catch errormsg
        save(fullfile(settings.outputdir,'errorinfo.mat'),'-v7.3')
        error('Error saved')
    end
    %datacalc{q}.ttvsig =
    
    % abs((prestim high real - pseudo)) - abs((prestim low real - pseudo))
    %datacalc{q}.nattvindex = abs((squeeze(trapz(datacalc{q}.naddttv.real(:,aucindex,2,:),2)) - squeeze(trapz(datacalc{q}.naddttv.pseudo(:,aucindex,2,:),2)))) - ...
    %    abs((squeeze(trapz(datacalc{q}.naddttv.real(:,aucindex,1,:),2)) - squeeze(trapz(datacalc{q}.naddttv.pseudo(:,aucindex,1,:),2))));
    datacalc{q}.naerspindex = abs((squeeze(trapz(datacalc{q}.naddersp.real(:,aucindex,2,:),2)) - squeeze(trapz(datacalc{q}.naddersp.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(datacalc{q}.naddersp.real(:,aucindex,1,:),2)) - squeeze(trapz(datacalc{q}.naddersp.pseudo(:,aucindex,1,:),2))));
    
    datacalc{q}.naerpindex = abs((squeeze(trapz(datacalc{q}.nadderp.real(:,aucindex,2,:),2)) - squeeze(trapz(datacalc{q}.nadderp.pseudo(:,aucindex,2,:),2)))) - ...
        abs((squeeze(trapz(datacalc{q}.nadderp.real(:,aucindex,1,:),2)) - squeeze(trapz(datacalc{q}.nadderp.pseudo(:,aucindex,1,:),2))));
end

if strcmpi(settings.datatype,'ECoG')
    for q = 1:numbands
        if strcmpi(settings.ecog.method,'mean') % take the mean of all electrodes
            fields = fieldnames_recurse(datacalc{q});
            fields = cell_unpack(fields);
            tmp = struct;
            for c = 1:length(fields)
                tmp = assignfield_nest(tmp,fields{c},nanmean(getfield_nest(datacalc{q},fields{c}),1));
            end
            datacalc{q} = tmp;
            
        elseif strcmpi(settings.ecog.method,'median')
            fields = fieldnames_recurse(datacalc{q});
            fields = cell_unpack(fields);
            tmp = struct;
            for c = 1:length(fields)
                tmp = assignfield_nest(tmp,fields{c},nanmedian(getfield_nest(datacalc{q},fields{c}),1));
            end
            datacalc{q} = tmp;
            
        elseif strcmpi(settings.ecog.method,'roi') % organize data into ROIs
            labels = ft_getlabels(timefreq_data{q},settings.datasetinfo.atlas);
            alllabels = settings.datasetinfo.atlas.tissuelabel; % designed for AAL atlas
            tmpdata = struct;
            fields = fieldnames_recurse(datacalc{q});
            fields = cell_unpack(fields);
            for c = 1:length(fields)
                tmp2 = [];
                for cc = 1:length(alllabels)
                    tmp = getfield_nest(datacalc{q},fields{c});
                    if ~isempty(find(strcmpi(labels,alllabels{cc})))
                        tmp = nanmedian(tmp(find(strcmpi(labels,alllabels{cc})),:,:,:),1); % for each field, take the mean across all electrodes with the same label
                    else
                        tmp = NaN(size(tmp)); % fill with NaNs for regions not represented
                        tmp = mean(tmp,1);
                    end
                    tmp2 = cat(1,tmp2,tmp);
                    %tmp2(cc,:,:,:) = tmp;
                end
                tmpdata = assignfield_nest(tmpdata,fields{c},tmp2);
            end
            datacalc{q} = tmpdata;
        end
    end
end

settings.pseudo.prestim = prestim_pseudo;
settings.real.prestim = prestim_real;
settings.pseudo.poststim = poststim_pseudo;
settings.real.poststim = poststim_real;

% Convert to single precision to save space
fields = fieldnames_recurse(datacalc{1});
fields = cell_unpack(fields);

for c = 1:length(fields)
    for cc = 1:length(datacalc)
        tmp = getfield_nest(datacalc{cc},fields{c});
        datacalc{cc} = assignfield_nest(datacalc{cc},fields{c},single(tmp));
    end
end

%% Saving the file
save(fullfile(settings.outputdir,[settings.datasetname '_' filename '_calc.mat']),'datacalc','-v7.3')

end
%
function Calc_naddPLE(settings,freqdata,filename)

if strcmpi(settings.tfparams.method,'hilbert') && ~strcmpi(settings.tfparams.method,'irasa')
    prestim_pseudo = settings.pseudo.prestim;
    prestim_real = settings.real.prestim;
    poststim_pseudo = settings.pseudo.poststim;
    poststim_real = settings.real.poststim;
else
    prestim_pseudo = settings.pseudo.prestim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    prestim_real = settings.real.prestim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    poststim_pseudo = settings.pseudo.poststim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
    poststim_real = settings.real.poststim - settings.pseudo.prestim(1)+1+ceil(settings.srate/5);
end

pledata = cell(1,size(freqdata.fourierspctrm,1));

for i = 1:size(freqdata.fourierspctrm,1)
    for ii = 1:size(freqdata.fourierspctrm,2)
        for iii = 1:size(freqdata.fourierspctrm,4)
            tmp = log10(freqdata.freq);
            pow = log10(squeeze(freqdata.fourierspctrm(i,ii,:,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            pledata{i}(1,ii,iii) = -p(1);
        end
    end
end

pledata = cat(1,pledata{:});

pledata = permute(pledata,[2 3 1]);


split_real = squeeze(mean(pledata(:,prestim_real,:),2));
split_pseudo = squeeze(mean(pledata(:,prestim_pseudo,:),2));

for c = 1:settings.nbchan
    splitindex = split_pseudo(c,:) > median(split_pseudo(c,:));
    
    naple.raw.pseudo(c,:,1) = mean(pledata(c,:,find(~splitindex)),3);
    naple.raw.pseudo(c,:,2) = mean(pledata(c,:,find(splitindex)),3);
    
    naple.pseudo(c,:,1) = (mean(pledata(c,poststim_pseudo,find(~splitindex)),3)...
        -mean(mean(pledata(c,prestim_pseudo,find(~splitindex)),3),2));
    naple.pseudo(c,:,2) = (mean(pledata(c,poststim_pseudo,find(splitindex)),3)...
        -mean(mean(pledata(c,prestim_pseudo,find(splitindex)),3),2));
    
    %         tmp = abs(datacat(c,poststim_pseudo,find(~splitindex)))...
    %             -mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),2);
    %         tmppseudo{1} = squeeze(trapz(tmp,2));
    %
    %         tmp = abs(datacat(c,poststim_pseudo,find(splitindex)))...
    %             -mean(abs(datacat(c,prestim_pseudo,find(splitindex))),2);
    %         tmppseudo{2} = squeeze(trapz(tmp,2));
    
    switch settings.units
        case 'prcchange'
            naple.pseudo(c,:,1) = 100*naple.pseudo(c,:,1)./mean(mean(pledata(c,prestim_pseudo,:),3),2);
            naple.pseudo(c,:,2) = 100*naple.pseudo(c,:,2)./mean(mean(pledata(c,prestim_pseudo,:),3),2);
        case 'zscore'
            naple.pseudo(c,:,1) = zscore(naple.pseudo(c,:,1),0,2);
            naple.pseudo(c,:,2) = zscore(naple.pseudo(c,:,2),0,2);
        case 'log'
            naple.pseudo(c,:,1) = 10*log10(naple.pseudo(c,:,1));
            naple.pseudo(c,:,2) = 10*log10(naple.pseudo(c,:,2));
    end
    
    splitindex = split_real(c,:) > median(split_real(c,:));
    
    naple.raw.real(c,:,1) = mean(pledata(c,:,find(~splitindex)),3);
    naple.raw.real(c,:,2) = mean(pledata(c,:,find(splitindex)),3);
    
    naple.real(c,:,1) = (mean(pledata(c,poststim_real,find(~splitindex)),3)...
        -mean(mean(pledata(c,prestim_real,find(~splitindex)),3),2));
    naple.real(c,:,2) = (mean(pledata(c,poststim_real,find(splitindex)),3)...
        -mean(mean(pledata(c,prestim_real,find(splitindex)),3),2));
    
    %         tmp = abs(datacat(c,poststim_pseudo,find(~splitindex)))...
    %             -mean(abs(datacat(c,prestim_pseudo,find(~splitindex))),2);
    %         tmpreal{1}(:) = squeeze(trapz(tmp,2));
    %
    %         tmp = abs(datacat(c,poststim_pseudo,find(splitindex)))...
    %             -mean(abs(datacat(c,prestim_pseudo,find(splitindex))),2);
    %         tmpreal{2}(:) = squeeze(trapz(tmp,2));
    
    %        [~,~,~,stats] = ttest2(tmpreal{2}-tmppseudo{2},tmpreal{1}-tmppseudo{1});
    %        naple{q}.t(c) = stats.t;
    
    switch settings.units
        case 'prcchange'
            naple.real(c,:,1) = 100*naple.real(c,:,1)./mean(mean(pledata(c,prestim_real,:),3),2);
            naple.real(c,:,2) = 100*naple.real(c,:,2)./mean(mean(pledata(c,prestim_real,:),3),2);
        case 'zscore'
            naple.real(c,:,1) = zscore(naple.real(c,:,1),0,2);
            naple.real(c,:,2) = zscore(naple.real(c,:,2),0,2);
        case 'log'
            naple.real(c,:,1) = 10*log10(naple.real(c,:,1));
            naple.real(c,:,2) = 10*log10(naple.real(c,:,2));
    end
    
    naple.diff = naple.real - naple.pseudo;
end


save(fullfile(settings.outputdir,[settings.datasetname '_' filename '_naddPLE.mat']),'naple','-v7.3')

end
%
