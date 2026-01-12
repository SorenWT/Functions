function [EEG_clean,thresh] = ft_autoreject_thresh(cfg,EEG,thresh)
% wrapper for autoreject threshold - applies to each channel
% independently

if EEG.nbchan > 1
    cfg = setdefault(cfg,'epochlen',2);
    cfg = setdefault(cfg,'removerej','yes');
    
    for i = 1:EEG.nbchan
        tmpEEG = pop_select(EEG,'channel',i);
        if nargin < 3
            [EEG_out{i},thresh(i)] = ft_autoreject_thresh(cfg,tmpEEG);
        else
            tmpthresh = thresh(i);
            [EEG_out{i}] = ft_autoreject_thresh(cfg,tmpEEG,tmpthresh);
        end
    end
    
    if strcmpi(cfg.removerej,'no')
        EEG_clean = EEG_out{1};
        EEG_clean.data = getfield_list(EEG_out,'data');
        EEG_clean.data = cat(1,EEG_clean.data{:});
        EEG_clean.chanlocs = getfield_list(EEG_out,'chanlocs'); EEG_clean.chanlocs = cat(1,EEG_clean.chanlocs{:});
        EEG_clean.nbchan = EEG.nbchan;
        EEG_clean = eeg_checkset(EEG_clean);
    else
        EEG_out = cat(1,EEG_out{:});
        EEG_clean = pop_mergeset(EEG_out,1:length(EEG_out));
    end
else
    
    
    data = eeglab2fieldtrip(EEG,'preprocessing','none');
    
    cfg = setdefault(cfg,'epochlen',2);
    cfg = setdefault(cfg,'removerej','yes');
    
    tmpcfg = []; tmpcfg.event = 1:cfg.epochlen*data.fsample:floor(length(data.time{1}));
    tmpcfg.event(end) = []; tmpcfg.epoch = [0 (cfg.epochlen*data.fsample)-1];
    data = ft_epoch(tmpcfg,data);
    data.trialinfo = ones(length(data.sampleinfo),1);
    
    EEG_epoch = eeg_regepochs(EEG,'recurrence',cfg.epochlen,'limits',[0 cfg.epochlen],'rmbase',[NaN])
    EEG_epoch.chanlocs.X = rand; EEG_epoch.chanlocs.Y = rand; EEG_epoch.chanlocs.Z = rand;
    EEG_epoch = eeg_checkset(EEG_epoch);
    
    
    if nargin < 3
            thisrand = randi(1000000);

    %save(['data_cont_epochs' num2str(thisrand) '.mat'],'data')
    pop_saveset(EEG_epoch,'filename',['cont_epochs_' num2str(thisrand) '.set'],'filepath',pwd);
    
    setenv('PATH','/opt/anaconda3/envs/mne/bin:/opt/anaconda3/condabin:/usr/local/fsl/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/Library/Apple/usr/bin:/Users/Soren/Library/Python/2.7/bin:/Users/Soren/abin:/Applications/workbench/bin_macosx64:/Users/Soren/.local/bin')
    
    pyscript = fopen(['tmp_pyscript_' num2str(thisrand) '.py'],'w');
    fprintf(pyscript,'import sys \n')
    fprintf(pyscript,'sys.path.insert(0, ''/Users/Soren/Documents/MATLAB/Functions'') \n')
    fprintf(pyscript,'from mne_preproc import autoreject_threshold \n')
    fprintf(pyscript,['autoreject_threshold(''' fullfile(pwd,['cont_epochs_' num2str(thisrand) '.set'])...
        ''',''' fullfile(pwd,['badsegs_' num2str(thisrand) '.json']) ''')'])
    fclose(pyscript)
    %system(['conda init'])
    %system(['conda activate mne'])
    system(['python tmp_pyscript_' num2str(thisrand) '.py'])
    system(['rm tmp_pyscript_' num2str(thisrand) '.py'])
    
    thresh = jsonread(fullfile(pwd,['badsegs_' num2str(thisrand) '.json']));
    thresh = thresh.eeg;
    end
    
    data.trial = data.trial(1:size(EEG_epoch.data,3));
    
    for i = 1:size(EEG_epoch.data,3)
        data.p2p(:,i) = max(data.trial{i},[],2)-min(data.trial{i},[],2);
        EEG_epoch.p2p(:,i) = max(EEG_epoch.data(:,:,i),[],2)-min(EEG_epoch.data(:,:,i),[],2);
    end
    bads = mean(EEG_epoch.p2p,1)>(thresh*1e6);
    EEG_epoch.etc.autoreject_thresh = thresh*1e6;
    EEG_epoch.etc.nrejected = sum(bads);
    tmp = bwconncomp(~bads);
    if length(tmp.PixelIdxList)>0
    EEG_epoch.etc.longestseg = cfg.epochlen*max(cellfun(@length,tmp.PixelIdxList));
    else
        EEG_epoch.etc.longestseg = 0;
    end
    
    
    if strcmpi(cfg.removerej,'yes')
        if all(bads)
            warning('All trials rejected - returning an empty dataset')
            EEG_clean = EEG;
            EEG_clean.data = [];
            EEG_clean = eeg_checkset(EEG_clean);
            return
        end
        
        EEG_epoch = pop_rejepoch(EEG_epoch,bads,0);
        EEG_epoch = eeg_checkset(EEG_epoch);
        
        EEG_clean = eeg_epoch2continuous(EEG_epoch);
        %EEG_clean.etc = EEG_epoch.etc;
    else
        % fill with NaN to keep different channels together
        
        for i = 1:size(EEG_epoch.data,3)
            if bads(i)
                EEG_epoch.data(:,:,i) = NaN;
                %data.trial{i} = NaN.*data.trial{i};
            end
        end
        
        EEG_clean = eeg_epoch2continuous(EEG_epoch);
        
%         data_clean = ft_concat(data);
%         data_clean = rmfield(data_clean,'trialinfo');
%         data_clean = rmfield(data_clean,'sampleinfo');
        
%         EEG_clean = EEG;
%         EEG_clean.data = double(data_clean.trial{1});
    end
end

end