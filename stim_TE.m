function tevals = stim_TE(braints,musicts,fs,labels,filename)

if size(braints,1) > size(braints,2)
   braints = braints'; 
end

if isstruct(labels)
    labels = extractfield(labels,'labels');
end

labels = vert(labels);

data = struct; 
data.trial{1} = [braints; horz(musicts)]; 
data.label = [labels; {'Music'}];
data.fsample = fs;
data.time{1} = linspace(0,length(braints)/fs,length(braints));

cfg = [];
cfg.event = 1:(fs*4):ceil(length(data.trial{1})-fs*4);
cfg.epoch = [0 fs*4];
cfg.unit = 'samples';
trl = ft_epoch(cfg,data);

cfg = [];
cfg.ensemblemethod = 'no';
cfg.sgncmb = cat(2,repmat({'Music'},length(labels),1),vert(labels));
cfg.toi = [0 4];
cfg.ragtaurange = [0.2 0.5]; cfg.ragdim = 2:12; cfg.actthrvalue = fs;
cfg.minnrtrials = 12; cfg.maxlag = fs; cfg.repPred = fs*0.5;
cfg.predicttime_u = 20;
cfg.flagNei = 'Mass'; cfg.sizeNei = 4;

dataprep = TEprepare(cfg,trl);

cfg = []; cfg.optdimusage = 'maxdim';
cfg.fileidout = filename;
cfg.surrogatetype = 'trialshuffling';
cfg.shifttest = 'no';
%cfg.numpermutation = 10000;

teresults = TEsurrogatestats(cfg,dataprep);



