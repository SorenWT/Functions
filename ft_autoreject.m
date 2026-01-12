function [data_clean] = ft_autoreject(cfg,data)
    % wrapper for autoreject taking data in Fieldtrip format
    
    cfg = setdefault(cfg,'epochlen',2);
    
    tmpcfg = []; tmpcfg.event = 1:cfg.epochlen*data.fsample:floor(length(data.time{1}));
    tmpcfg.event(end) = []; tmpcfg.epoch = [0 (cfg.epochlen*data.fsample)-1];
    data = ft_epoch(tmpcfg,data);
    data.trialinfo = ones(length(data.sampleinfo),1);
    
    save('data_cont_epochs.mat','data')
    
    setenv('PATH','/opt/anaconda3/envs/mne/bin:/opt/anaconda3/condabin:/usr/local/fsl/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/Library/Apple/usr/bin:/Users/Soren/Library/Python/2.7/bin:/Users/Soren/abin:/Applications/workbench/bin_macosx64:/Users/Soren/.local/bin')
    
    pyscript = fopen(['tmp_pyscript.py'],'w');
    fprintf(pyscript,'import sys \n')
    fprintf(pyscript,'sys.path.insert(0, ''/Users/Soren/Documents/MATLAB/Functions'') \n')
    fprintf(pyscript,'from mne_preproc import autoreject_log \n')
    fprintf(pyscript,['autoreject_log(''' fullfile(pwd,'data_cont_epochs.mat')...
        ''',''' fullfile(pwd,'badsegs.json') ''')'])
    system(['conda init'])
    system(['conda activate mne'])
    system(['python tmp_pyscript.py'])
    system(['rm tmp_pyscript.py'])
    
    bads = jsonread(fullfile(pwd,'badsegs.json'));
    
    
    tmpcfg = []; tmpcfg.trials = ~bads;
    data = ft_selectdata(tmpcfg,data);
    
    data_clean = ft_concat(data);
    data_clean = rmfield(data_clean,'trialinfo');
    data_clean = rmfield(data_clean,'sampleinfo');
    
end