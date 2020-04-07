function teresults = stim_TE_wentao(braints,musicts,fs,labels,lags,shiftstats)
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
%     shiftstats: 'yes' or 'no', normalize single-trial values by
%        time-shifted surrogates (default = 'yes' if only one trial input -
%        otherwise trial-shuffling is implemented instead)
%
% Outputs: 
%     teresults: transfer entropy results. Important fields are:
%         TEmat: the actual transfer entropy values. The format is channels
%         x music pieces, or channels x shifts for single-trial data
%         MImat: mutual information values. Format is the same as TEmat
%         (channels x music pieces/shifts)


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
    if size(braints,3) == 1
        shiftstats = 'yes';
    else
        shiftstats = 'no';
    end
end

labels = vert(labels);

data = struct; 
%data.trial{1} = [braints; horz(musicts)]; 
for i = 1:size(braints, 3)
   data.trial{i} = [braints(:,:,i); musicts(i,:)]; 
   data.time{i} = linspace(0,length(braints)/fs,length(braints));
end

if size(braints,3) == 1
   shifttimes = [(2*max(lags)):(2*max(lags)+5*length(lags)-1)];
   for i = 1:length(shifttimes)
        data.trial{i+1} = [braints; musicts([(shifttimes(i)+1):end 1:(shifttimes(i))])];
        data.time{i+1} = linspace(0,length(braints)/fs,length(braints));
   end
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
    cfg.ragtaurange = [0.2 0.4]; cfg.ragdim = 2:8; cfg.actthrvalue = length(musicts)/3;
    cfg.minnrtrials = length(data.trial); cfg.maxlag = ceil(0.5*length(musicts)); cfg.repPred = ceil(0.2*length(musicts));
    cfg.predicttime_u = lags(i);
    cfg.flagNei = 'Mass'; cfg.sizeNei = 4;
    
    dataprep = TEprepare(cfg,data);
    
    if strcmpi(shiftstats,'no')
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
        teresults{i}.te_z = zscore(teresults{i}.TEmat,[],2);
        teresults{i}.mi_z = zscore(teresults{i}.MImat,[],2);
    end
end



