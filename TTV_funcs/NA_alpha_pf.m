function settings = NA_alpha_pf(settings)

cd(settings.inputdir)

files = dir(settings.files);

pool = gcp('nocreate');

if isempty(pool) || pool.NumWorkers ~= settings.pool
    delete(gcp('nocreate'))
    parpool(settings.pool)
end

fbands = cell(1,length(files));
pf = cell(1,length(files));

parfor i = 1:length(files)
    if strcmpi(settings.datatype,'EEG')
        EEG = pop_loadset('filename',files(i).name,'filepath',settings.inputdir);
        data = eeglab2fieldtrip(EEG,'preprocessing','none');
        data = ft_struct2single(data);
    else
        data = parload(fullfile(files(i).folder,files(i).name),'data');
    end
        
   [fbands{i},pf{i}] = NA_convert_alpha_pf(settings,ft_concat(data));
end

fbands = cat(1,fbands{:});
pf = cat(1,pf{:});

settings.tfparams.fbands = fbands;
settings.alpha_pf = pf;
