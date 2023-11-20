function [data] = mbian_preproc(cfg,data)
% mbian_preproc implements a standard preprocessing pipeline for the Mind,
% Brain Imaging and Neuroethics group at the Royal
%
% Input arguments:
%      cfg: a config structure with the following fields
%           datatype: 'eeg' or 'meg' (default = 'eeg')
%           format: 'eeglab' or 'fieldtrip' (default = 'eeglab' for eeg,
%              'fieldtrip for meg)
%           bandpass: [low_freq high_freq] (default = [0.5 50] for eeg,
%              [0.05 150] for meg
%           resample: number specifying the sample rate for the
%              preprocessed recording (default = 500 Hz or the original
%              sampling rate, whichever is lower)
%
%           Datatype-specific options:
%
%           EEG:
%
%           chanlocs: either 'lookup' to look up channel locations in a
%              standard file, a string specifying the file to look up
%              locations in, 'none', or an EEGLAB chanlocs structure
%              (default = 'lookup')
%           line: structure with the following fields:
%              method: Can be 'cleanline', 'none', or 'notch' (default =
%                 'cleanline')
%              freq: an n x m matrix of line noise frequencies. For method
%                 'cleanline', input only the actual line frequencies (m =
%                 1). For method 'notch', input bands you want to notch
%                 filter (m = 2).
%           asr: clean data using artifact subspace reconstruction - 'yes'
%              or 'no' (default = 'yes')
%           reference: referencing scheme - 'default','average', or 'REST'
%              (default = 'average')
%           ica: decompose data using ICA - 'yes' or 'no' (default = 'yes')
%           mara: automatic ICA component rejection using MARA - 'yes' or
%              'no' (default = 'yes')
%
%           MEG:
%
%           not implemented for now
%
%      data: an EEGLAB or Fieldtrip data structure
%
% Outputs:
%      data: the preprocessed data structure
%
%

%% Set up defaults
cfg = setdefault(cfg,'datatype','eeg');
if ~cfgcheck(cfg,'format')
    if cfgcheck(cfg,'datatype','eeg')
        cfg.format = 'eeglab';
    elseif cfgcheck(cfg,'datatype','meg')
        cfg.format = 'fieldtrip';
    end
end

if cfgcheck(cfg,'datatype','eeg')
    cfg = setdefault(cfg,'bandpass',[0.5 50]);
elseif cfgcheck(cfg,'datatype','meg')
    cfg = setdefault(cfg,'bandpass',[0.05 150]);
end

if cfgcheck(cfg,'datatype','eeg') && data.srate > 500
    cfg = setdefault(cfg,'resample',500);
elseif cfgcheck(cfg,'datatype','meg') && data.fsample > 500
    cfg = setdefault(cfg,'resample',500);
end

if cfgcheck(cfg,'datatype','eeg')
    cfg = setdefault(cfg,'chanlocs','lookup');
    cfg = setdefault(cfg,'line','cleanline');
    cfg = setdefault(cfg,'asr','yes');
    cfg = setdefault(cfg,'reference','average');
    cfg = setdefault(cfg,'ica','yes');
    cfg = setdefault(cfg,'mara','yes');
elseif cfgcheck(cfg,'datatype','meg')
    % add later
end

ft_defaults

%% EEG pipeline
if cfgcheck(cfg,'datatype','eeg')
    %eeglab rebuild
    
    eegdir = extractBefore(which('eeglab'),'eeglab.m');
    
    EEG = data;
    data = [];
    
    % Read in channel locations
    if isstr(cfg.chanlocs)
        if cfgcheck(cfg,'chanlocs','lookup')
            %EEG=pop_chanedit(EEG, {'lookup','/Users/Soren/Documents/MATLAB/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc'});
            EEG = pop_chanedit(EEG,'lookup',fullfile(char(eegdir),'plugins','dipfit','standard_BEM','elec','standard_1005.elc'),'eval','chans = pop_chancenter( chans, [],[]);');
        elseif ~cfgcheck(cfg,'chanlocs','no')
            EEG = pop_chanedit(EEG,'lookup',cfg.chanlocs,'eval','chans = pop_chancenter( chans, [],[]);');
        end
    else
        EEG.chanlocs = cfg.chanlocs;
    end
    
    EEG  = eeg_checkset(EEG);
    
    % Resample the data
    
    if ~strcmpi(cfg.resample,'no') && (EEG.srate ~= cfg.resample)
        EEG = pop_resample(EEG,cfg.resample);
    end
    
    
    % Filter the data
    % filter in fieldtrip so you don't have to manually set the filter order
    chanlocs = EEG.chanlocs; %ft2eeglab doesn't handle chanlocs well, so save these
    
    if ~ischar(cfg.bandpass) && ~strcmpi(cfg.bandpass,'no')
        EEG = pop_eegfiltnew(EEG, cfg.bandpass(1),cfg.bandpass(2),[],0,[],0);
        %event = EEG.event;
        %data = eeglab2fieldtrip(EEG,'preprocessing','none');
        %tmpcfg = [] ; tmpcfg.bpfilter = 'yes'; tmpcfg.bpfreq = cfg.bandpass; tmpcfg.bpfilttype = 'fir';
        %data = ft_preprocessing(tmpcfg,data);
        %EEG = ft2eeglab(data);
        %EEG.chanlocs = chanlocs;
        %EEG.event = event;
    end
    
    EEG  = eeg_checkset(EEG);
    
    
    % Clean line noise
    switch cfg.line.method
        case 'cleanline'
            EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan,'computepower',1,'linefreqs',cfg.line.freq ,'normSpectrum',1,'p',0.05,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
            delete(gcf)
        case 'notch'
            for q = 1:size(cfg.line.freq,1)
                EEG = pop_eegfiltnew(EEG, cfg.line.freq(q,1),cfg.line.freq(q,2),[],1,[],0);
            end
            %chanlocs = EEG.chanlocs; %ft2eeglab doesn't handle chanlocs well, so save these
            %event = EEG.event;
            %data = eeglab2fieldtrip(EEG,'preprocessing','none');
            %tmpcfg = [] ; tmpcfg.bsfilter = 'yes'; tmpcfg.bsfreq = cfg.line.freq; tmpcfg.bpfilttype = 'fir';
            %data = ft_preprocessing(tmpcfg,data);
            %EEG = ft2eeglab(data);
            %EEG.chanlocs = chanlocs;
            %EEG.event = event;
    end
    EEG  = eeg_checkset(EEG);
    
    
    % Apply ASR
    if cfgcheck(cfg,'asr','yes')
        orig_chanlocs = EEG.chanlocs;
        
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.9,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
        %EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
        %EEG = clean_rawdata(EEG, 5, [-1], 0.85, 4, 20,-1,'availableRAM_GB',4);
        
        EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
        
        EEG  = eeg_checkset(EEG);
        
    end
    
    % Rereference
    switch cfg.reference
        case 'average'
            EEG = pop_reref( EEG, []);
        case 'REST'
            % add REST reference
    end
    EEG  = eeg_checkset(EEG);
    
    
    % Run ICA
    if cfgcheck(cfg,'ica','yes')
        EEG = pop_runica(EEG, 'interupt','off');
        EEG  = eeg_checkset(EEG);
    end
    
    % Reject components with MARA
    if cfgcheck(cfg,'mara','yes')
        [artcomps,info] = MARA(EEG);
        EEG  = pop_subcomp(EEG,artcomps,0);
        EEG  = eeg_checkset(EEG);
        EEG.etc.marainfo = info;
    end
    
    data = EEG;
    EEG = [];
    
elseif cfgcheck(cfg,'datatype','meg')
    %% MEG pipeline
end


