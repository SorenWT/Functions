function [sources,roidata,voxeldata] = SourceEst_HCP(data,subid,atlas,parflag,sources)

headmodel = parload(['/data/hcp/meg/' subid '/anatomy/' subid '_MEG_anatomy_headmodel.mat'],'headmodel');
sourcemodel = parload(['/data/hcp/meg/' subid '/anatomy/' subid '_MEG_anatomy_sourcemodel_2d.mat'],'sourcemodel2d');
loretaflag = 1;
if ~loretaflag
% %defining trials for noise data
cfg = [];
cfg.dataset = ['/data/hcp/meg/' subid '/unprocessed/1-Rnoise/4D/c,rfDC'];
cfg.trialdef.trialDuration = 2;
trl = trialfun_Restin(cfg);

%reading noise data, band stop filtering for line noise
cfg = []; cfg.dataset = ['/data/hcp/meg/' subid '/unprocessed/1-Rnoise/4D/c,rfDC'];
cfg.bsfilter = 'yes'; cfg.bsfreq = [59 119 179 239 299;61 121 181 241 301];
cfg.trl = trl;

%resampling noise data
noisedata = ft_preprocessing(cfg);
cfg = []; cfg.resamplefs = 508.6275; cfg.detrend = 'no';
noisedata = ft_resampledata(cfg,noisedata);

% Select only the MEG channels
cfg = []; cfg.channel = {'MEG'};
data = ft_selectdata(cfg,data);

%calculating noise covariance
cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-Inf Inf];
cfg.channel = data.label;
noise_avg = ft_timelockanalysis(cfg,noisedata);
else
    noise_avg.cov = ones(length(data.label));
end

%[~,ftpath] = ft_version;
%ftpath = '/group/northoff/share/fieldtrip-master';

% % Convert the data to Fieldtrip format
% data = eeglab2fieldtrip(EEG,'preprocessing','none');

% Load headmodel, sourcemodel, atlas
% 
% if isstr(headmodel)
%     load(headmodel) %variable is "vol"
% end
% 
% if isstr(sourcemodel)
%     sourcemodel = ft_read_headshape(sourcemodel);
% end

% if isstr(atlas)
%     atlas = ft_read_atlas(atlas);
% end

% % Interpolate template surface on atlas
% cfg = [];
% cfg.interpmethod = 'nearest';
% cfg.parameter = 'tissue';
% sourcemodel = ft_sourceinterpolate(cfg, atlas, sourcemodel);

if ~exist('sources','var')
data.grad = ft_convert_units(data.grad,'mm');
% Create forward model
cfg = []; cfg.grad = data.grad; cfg.channel = {'MEG'};
sourcemodel = ft_convert_units(sourcemodel,'mm');
headmodel = ft_convert_units(headmodel,'mm');
cfg.grid.pos = sourcemodel.pos; cfg.grid.inside = 1:size(sourcemodel.pos,1);
cfg.headmodel = headmodel;
leadfield = ft_prepare_leadfield(cfg);

% clear data

% % Epoch into arbitrary 2-second segments
% cfg = []; cfg.event = 1:2*data.fsample:length(data.time{1}); cfg.epoch = [0 (2*data.fsample)-1];
% cfg.event(end) = [];
% data = ft_epoch(cfg,data);


cfg = []; cfg.keeptrials = 'yes';
tlock = ft_timelockanalysis(cfg,data);
tlock.cov = permute(repmat(noise_avg.cov,1,1,length(data.trial)),[3 1 2]);

% Estimate sources
cfg = [];
cfg.method = 'eloreta';
cfg.grid = leadfield;
cfg.eloreta.keepfilter = 'yes';
cfg.eloreta.normalize = 'yes';
cfg.headmodel = headmodel;
sources = ft_sourceanalysis(cfg,tlock);
end




if nargout > 1
cfg = []; cfg.channel = {'MEG'};
data = ft_selectdata(cfg,data);

% Get voxel time courses
voxeldata = struct;
datacat = cat(2,data.trial{:});
source_datacat = zeros(size(sources.pos,1),size(datacat,2));
if parflag
    filter = sources.avg.filter;
    parfor c = 1:size(sources.pos,1)
        tmp = filter{c}*datacat;
        [u,s,v] = svd(tmp,'econ');
        source_datacat(c,:) = u(:,1)'*filter{c}*datacat;
    end
else
    for c = 1:size(sources.pos,1)
        tmp = sources.avg.filter{c}*datacat;
        [u,s,v] = svd(tmp,'econ');
        source_datacat(c,:) = u(:,1)'*sources.avg.filter{c}*datacat;
    end
end
%voxeldata.trial{1} = source_datacat;
voxeldata.time{1} = linspace(0,size(source_datacat,2)/data.fsample-1/data.fsample,size(source_datacat,2));
voxeldata.label = cellstr(num2str([1:size(sources.pos,1)]'));
voxeldata.fsample = data.fsample;
%source_datacat = []; %saving memory

% Get ROI time courses
roidata = voxeldata;
roidata.label = atlas.parcellationlabel;
roidata.trial = [];
trialdat = zeros(length(roidata.label),size(source_datacat,2));
for cc = 1:length(roidata.label)
    trialdat(cc,:) = mean(source_datacat(atlas.parcellation == cc,:),1);
end
roidata.trial{1} = trialdat;
voxeldata.trial{1} = source_datacat;
source_datacat = []; trialdat = []; 
end
end
