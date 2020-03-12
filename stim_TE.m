function teresults = stim_TE(braints,musicts,fs,labels,lags,stats)
% Inputs: 
%     braints: a 2D or 3D matrix of the measure time series in the form
%        channels x time points x music pieces (can do only one music piece
%        if you want)
%     musicts: a row vector or 2D matrix of music measure time series in
%        the form music pieces x time points (should have same order as the 
%        3rd dimension of braints)
%     fs: the sampling frequency of your data (should be equal to 1 over 
%        the step size for measure time series)
%     labels: either a cell array containing the channel names for your
%        braints data, or an EEGLAB chanlocs structure (the channel names
%        will be extracted from the latter)
%     lags: a scalar or vector of lags (in samples) at which to compute
%        transfer entropy (default = 1:ceil(fs))
%     stats: 'yes' or 'no', do statistics on the data (only put 'yes' if
%        you are inputting a 3D matrix with multiple music pieces; default
%        = 'no')
%
% Outputs: 
%     teresults: transfer entropy results. Important fields are:
%         TEmat: the actual transfer entropy values. The format is channels
%         x music pieces
%         MImat: mutual information values. Format is the same as TEmat
%         (channels x music pieces)


if size(braints,1) > size(braints,2)
   braints = braints'; 
end

if isstruct(labels)
    labels = extractfield(labels,'labels');
end

if ~exist('lags','var')
   lags = 1:ceil(fs); 
end

if ~exist('stats','var')
   stats = 'no'; 
end

labels = vert(labels);

data = struct; 
%data.trial{1} = [braints; horz(musicts)]; 
for i = 1:size(braints, 3)
   data.trial{i} = [braints(:,:,i); musicts(i,:)]; 
   data.time{i} = linspace(0,length(braints)/fs,length(braints));
end
data.label = [labels; {'Music'}];
data.fsample = fs;

% cfg = [];
% cfg.event = 1:(fs*4):ceil(length(data.trial{1})-fs*4);
% cfg.epoch = [0 fs*4];
% cfg.unit = 'samples';
% trl = ft_epoch(cfg,data);

%trl = data;

for i = 1:length(lags)
    cfg = [];
    cfg.ensemblemethod = 'no';
    cfg.sgncmb = cat(2,repmat({'Music'},length(labels),1),vert(labels));
    cfg.toi = [-Inf Inf];
    cfg.ragtaurange = [0.2 0.5]; cfg.ragdim = 2:12; cfg.actthrvalue = length(musicts);
    cfg.minnrtrials = length(data.trial); cfg.maxlag = 0.5*length(musicts); cfg.repPred = ceil(0.2*length(musicts));
    cfg.predicttime_u = lags(i);
    cfg.flagNei = 'Mass'; cfg.sizeNei = 4;
    
    dataprep = TEprepare(cfg,data);
    
    if strcmpi(stats,'yes')
        cfg = []; cfg.optdimusage = 'maxdim';
        cfg.fileidout = 'test';
        cfg.surrogatetype = 'trialshuffling';
        cfg.shifttest = 'no';
        cfg.numpermutation = 1000;
        
        teresults{i} = TEsurrogatestats(cfg,dataprep);
    else
        cfg.tau(1:size(dataprep.TEprepare.channelcombi,1)) = dataprep.TEprepare.opttau;
        cfg.kth_neighbors = 4; cfg.TheilerT = 'ACT'; cfg.extracond = 'none';
        cfg.dim(1:size(dataprep.TEprepare.optdimmat,1),1) = dataprep.TEprepare.optdim;
        cfg.calctime = 'no'; cfg.shuffle = 'no';
        
        teresults{i} = transferentropy(cfg,dataprep);
    end
end



