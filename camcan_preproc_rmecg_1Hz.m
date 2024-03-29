function [data,cont_data,comp_class,iteration] = camcan_preproc(subid,filename,cont_data,event)

% Assuming ComputeCanada
basedir = extractBefore(filename,['sub-' subid]);
basedir = char(basedir);
do_epoch = 0;
if contains(basedir,'task')
    do_epoch = 1;
elseif contains(basedir,'Rest')
    do_epoch = 0;
elseif exist('cont_data','var')
    do_epoch = 1;
end
%basedir = '/scratch/sorenwt/camcan/cc700/mri/pipeline/release004/BIDSsep/megmax_task';

% Option to input the continuous data in case something goes wrong and you
% only want to do the epoched preprocessing
if ~exist('cont_data','var')
    
    %% First do movement compensation in MNE somehow
    
    %% Load in the file
    hdr = ft_read_header(filename);
    
    
    cfg = []; cfg.dataset = filename;
    cfg.hdr = filename;
    cont_data = ft_preprocessing(cfg);
    
    
    
    %% Resample to 500 Hz
    
    cfg = []; cfg.resamplefs = 500;
    cont_data = ft_resampledata(cfg,cont_data);
    
    %% Cut data length to 5 mins
    
    cfg = []; cfg.latency = [60 360];
    cont_data = ft_selectdata(cfg,cont_data);
    
    %% Filter
    
    cfg = []; cfg.lpfilter = 'yes'; cfg.lpfreq = 150;
    cfg.bsfilter = 'yes'; cfg.bsfreq = [49 51; 99 101; 149 151];
    cont_data = ft_preprocessing(cfg,cont_data);

    cfg = []; cfg.hpfilter = 'yes'; cfg.hpfreq = 1;
    cont_data = ft_preprocessing(cfg,cont_data);
    
    %% Consider only magnetometers (and EOG) for now
    
    cfg = []; cfg.channel = {'megmag','eog','-eog062'};
    cont_data = ft_selectdata(cfg,cont_data);
    
    
    %% Use autoreject to find threshold to exclude from ICA
    
    % Filter at 1 Hz first for ICA later
   % cfg = []; cfg.hpfilter = 'yes'; cfg.hpfreq = 1;
   % data = ft_preprocessing(cfg,cont_data);
    data = cont_data;
    
    cfg = []; cfg.event = 1:2*cont_data.fsample:floor(length(cont_data.time{1}));
    cfg.event(end) = []; cfg.epoch = [0 (2*cont_data.fsample)-1];
    data = ft_epoch(cfg,data);
    data.trialinfo = ones(length(data.sampleinfo),1);
    
    cfg = []; cfg.channel = {'EOG'};
    refdata = ft_selectdata(cfg,data);
    
    cfg = []; cfg.channel = {'megmag'};
    data = ft_selectdata(cfg,data);
    
    save(fullfile(basedir,['sub-' subid],[subid '_cont_epochs.mat']),'data')
    
    pyscript = fopen([subid '_pyscript.py'],'w');
    fprintf(pyscript,'import sys \n')
    fprintf(pyscript,'sys.path.insert(0, ''/home/soren/Documents/MATLAB/Functions'') \n')
    fprintf(pyscript,'sys.path.insert(0, ''/home/sorenwt/projects/def-gnorthof/sorenwt/MATLAB/Functions'') \n')
    fprintf(pyscript,'from mne_preproc import autoreject_log \n')
    fprintf(pyscript,['autoreject_log(''' fullfile(basedir,['sub-' subid],[subid '_cont_epochs.mat'])...
        ''',''' fullfile(basedir,['sub-' subid],[subid '_badsegs.json']) ''')'])
    system(['python ' subid '_pyscript.py'])
    system(['rm ' subid '_pyscript.py'])
    
    bads = jsonread(fullfile(basedir,['sub-' subid],[subid '_badsegs.json']));
    
    %% Remove bad segments and re-concatenate
    
    cfg = []; cfg.trials = ~bads;
    data = ft_selectdata(cfg,data);
    
    cfg = []; cfg.trials = ~bads;
    refdata = ft_selectdata(cfg,refdata);
    
    cont_data_clean = ft_concat(data);
    cont_data_clean = rmfield(cont_data_clean,'trialinfo');
    cont_data_clean = rmfield(cont_data_clean,'sampleinfo');
    
        
    refdata = ft_concat(refdata);
    refdata = rmfield(refdata,'trialinfo');
    refdata = rmfield(refdata,'sampleinfo');
    
    %cfg = []; cfg.hpfilter = 'yes'; cfg.hpfreq = 1;
    %cont_data_clean = ft_preprocessing(cfg,cont_data_clean);
    
    %cont_data_clean = ft_concat(cont_data_clean);
    
    %%% Rescale planar gradiometers
    
    %     mag_cutoff = 62;
    %     plan_cutoff = 62;
    %
    %     chan_inds = ft_channelselection({'MEG'},cont_data_clean.grad);
    %     chan_inds = Subset_index(cont_data_clean.label,chan_inds);
    %     norm_vec = ones(numel(chan_inds),1);
    %     icadata = cont_data_clean.trial{1}(chan_inds,:);
    %
    %     if any(strcmp(cont_data_clean.grad.chantype,'megmag')) && any(strcmp(cont_data_clean.grad.chantype,'megplanar'))
    %         mag_min_eig = svd(cov(icadata(strcmp(cont_data_clean.grad.chantype(chan_inds),'megmag'),:)'));
    %         mag_min_eig = mean(mag_min_eig(mag_cutoff-2:mag_cutoff));
    %
    %         plan_min_eig = svd(cov(icadata(strcmp(cont_data_clean.grad.chantype(chan_inds),'megplanar'),:)'));
    %         plan_min_eig = mean(plan_min_eig(plan_cutoff-2:plan_cutoff));
    %
    %         norm_vec(strcmp(cont_data_clean.grad.chantype(chan_inds),'megmag'))    = mag_min_eig;
    %         norm_vec(strcmp(cont_data_clean.grad.chantype(chan_inds),'megplanar')) = plan_min_eig;
    %     else
    %         norm_vec = norm_vec*min(svd(cov(icadata(:,:)')));
    %     end
    %     norm_vec = sqrt(norm_vec);
    %
    %     icadata = icadata ./ repmat(norm_vec,1,size(icadata,2));
    %
    %     cont_data_clean.trial{1} = icadata;
    %     icadata = [];
    %
    %     cfg = []; cfg.channel = chan_inds;
    %     cont_data_clean_meg
    
    %     %% ICA using OSL-AFRICA
    %
    %     D = spm_eeg_ft2spm(cont_data_clean,fullfile(basedir,['sub-' subid],'tmp'));
    %     D = osl_africa(D,'used_maxfilter',true);
    %     system(['rm ' fullfile(basedir,['sub-' subid],'tmp*')])
    % %     cont_data_clean = spm2fieldtrip(D);
    % %     cont_data_clean.hdr = hdr; cont_data_clean.fsample = hdr.Fs;
    %     cont_data_clean = [];
    
    %% ICA using HCP scripts
    
    options.bandpass = [1 200];
    options.bandstop = [49 51; 99 101; 149 151; 199 201];
    [iteration,~] = hcp_ICA_unmix(cont_data_clean,{'channel','MEG','ica_iterations',20,'numIC',62});
    comp_class = hcp_ICA_RMEG_classification(refdata,options,iteration,cont_data_clean);
    cont_data_clean = [];
    
    %% Apply ICA matrix back to original data
    
    cfg = []; cfg.channel = {'megmag'};
    cont_data = ft_selectdata(cfg,cont_data);
    
    cfg = []; cfg.unmixing = comp_class.unmixing; cfg.topolabel = cont_data.label;
    comp = ft_componentanalysis(cfg,cont_data);
    
    cfg = []; cfg.component = except(1:62,comp_class.class.brain_ic);
    cont_data = ft_rejectcomponent(cfg,comp,cont_data);
    %    comp_class = []; comp = [];
    
end

%% Epoch data around stimuli
if do_epoch
    if ~exist('event','var')
    event = ft_read_event(filename);
    end
    types = extractfield(event,'type');
    tpoints = extractfield(event,'sample');
    latencies = tpoints(find(strcmpi(types,'Trigger')));
    allvalues = extractfield(event,'value');
    values = allvalues(find(strcmpi(types,'Trigger')));
    latencies(find(values < 0)) = [];
    
    rtevents = find(strcmpi(types,'Trigger'));
    rtevents(find(values < 0)) = [];
    rtlatencies = tpoints(rtevents+1);
    rts = rtlatencies-latencies;
    latencies = round(latencies/2);
    
    cfg = []; cfg.event = latencies; cfg.epoch = [-2*cont_data.fsample 1.5*cont_data.fsample]; %Epochs large in order to have trial padding
    data = ft_epoch(cfg,cont_data);
    data.trialinfo(:,1) = vert(values(find(values >= 0)));
    data.trialinfo(:,2) = vert(rts);
    data.grad = data.hdr.grad;
    
    sampinfo = data.sampleinfo;
    trialinfo = data.trialinfo;
else
    cfg = []; cfg.event = 1:2*cont_data.fsample:floor(length(cont_data.time{1}));
    cfg.event(end) = []; cfg.epoch = [0 (2*cont_data.fsample)-1];
    data = ft_epoch(cfg,cont_data);
    data.trialinfo = ones(length(data.sampleinfo),1);
    data.grad = cont_data.hdr.grad;
end
%% Use autoreject to remove or interpolate bad epochs
save(fullfile(basedir,['sub-' subid],[subid '_ftdata.mat']),'data');
data = [];

pyscript = fopen([subid '_pyscript.py'],'w');
fprintf(pyscript,'import sys \n')
fprintf(pyscript,'sys.path.insert(0, ''/home/soren/Documents/MATLAB/Functions'') \n')
fprintf(pyscript,'sys.path.insert(0, ''/home/sorenwt/projects/def-gnorthof/sorenwt/MATLAB/Functions'') \n')
fprintf(pyscript,'from mne_preproc import autoreject_epochs \n')
% fprintf(pyscript,'from mne_preproc import autoreject_log \n')
% if do_epoch
%     fprintf(pyscript,['autoreject_log(''' fullfile(basedir,['sub-' subid],[subid '_ftdata.mat'])...
%         ''',''' fullfile(basedir,['sub-' subid],[subid '_badsegs.json']) ''') \n'])
% end
fprintf(pyscript,['autoreject_epochs(''' fullfile(basedir,['sub-' subid],[subid '_ftdata.mat'])...
    ''',''' fullfile(basedir,['sub-' subid],[subid '_mne.fif']) ''','''...
    fullfile(basedir,['sub-' subid],[subid '_badsegs.json']) ''')'])
system(['python ' subid '_pyscript.py'])
system(['rm ' subid '_pyscript.py'])

if do_epoch
    bads = jsonread(fullfile(basedir,['sub-' subid],[subid '_badsegs.json']));
end

cfg = []; cfg.dataset = fullfile(basedir,['sub-' subid],[subid '_mne.fif']);
cfg.hdr = cfg.dataset;
data = ft_preprocessing(cfg);

if do_epoch
    data.trialinfo = trialinfo(find(~bads),:);
    data.sampleinfo = sampinfo(find(~bads),:);
end
