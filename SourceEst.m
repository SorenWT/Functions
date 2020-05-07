function [roidata,voxeldata,sources] = SourceEst(data,headmodel,sourcemodel,atlas,opts)
% SourceEst does virtual-channel source estimation for one subject, giving
% a time-series of source-level activity for each point in the sourcemodel
% input, or for each region in the atlas input
%
% Inputs: 
%    data: the data for source estimation, in Fieldtrip format
%    headmodel: the headmodel for the subject
%    sourcemodel: the sourcemodel for the subject
%    atlas: an atlas dividing the sourcemodel into regions
%    opts: opts can have the following fields
%       atlasparam: the field of the atlas containing the parcellation
%         information (default = 'parcellation')
%       atlaslabel: the field of the atlas containing the labels for each
%         region (default = 'parcellationlabel')
%       interp: interpolate the sourcemodel onto the atlas ('yes' or 'no' -
%         default = 'no', meaning it is assumed the atlas and sourcemodel 
%         are already aligned)
%       datatype: 'EEG' or 'MEG' - for channel selection (default = 'MEG')
%       noisecov: noise covariance matrix, if desired
%       method: source imaging method: currently supports 'eloreta' or 'mne'
%         (default = 'eloreta')
%
% Outputs: 
%    roidata: a Fieldtrip-format data structure containing the source-level
%       time series of each region in atlas
%    voxeldata: a Fieldtrip-format data structure containing the
%       source-level time series of each voxel or vertex in sourcemodel
%    sources: the original Fieldtrip source-estimation structure. This
%       contains the spatial filter used to compute the source-level time 
%       series 
% 
% Created by Soren Wainio-Theberge 

if ~exist('atlas','var')
   atlas = []; 
end

if ~exist('opts','var')
    opts = struct;
end

opts = setdefault(opts,'atlasparam','parcellation');
opts = setdefault(opts,'atlaslabel','parcellationlabel');
opts = setdefault(opts,'interp','no');
opts = setdefault(opts,'datatype','MEG');
opts = setdefault(opts,'noisecov',eye(length(data.label)));
opts = setdefault(opts,'method','eloreta');

% Load headmodel, sourcemodel, atlas if they're strings

if isstr(headmodel)
    load(headmodel) %variable is "vol"
end

if isstr(sourcemodel)
    sourcemodel = ft_read_headshape(sourcemodel);
end

if isstr(atlas)
    atlas = ft_read_atlas(atlas);
end

% % Interpolate template surface on atlas
if strcmpi(opts.interp,'yes') && ~isempty(atlas)
    cfg = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter = opts.atlasparam;
    sourcemodel = ft_sourceinterpolate(cfg, atlas, sourcemodel);
    atlas.(opts.atlasparam) = sourcemodel.(opts.atlasparam);
end

% Create forward model
cfg = [];
if strcmpi(opts.datatype,'MEG')
    cfg.grad = data.grad;
elseif strcmpi(opts.datatype,'EEG')
    cfg.elec = data.elec;
end
cfg.channel = {opts.datatype};
sourcemodel = ft_convert_units(sourcemodel,'mm');
headmodel = ft_convert_units(headmodel,'mm');
cfg.grid.pos = sourcemodel.pos; cfg.grid.inside = 1:size(sourcemodel.pos,1);
cfg.headmodel = headmodel;
leadfield = ft_prepare_leadfield(cfg);

% % Select only the relevant channels
cfg = []; cfg.channel = {opts.datatype};
data = ft_selectdata(cfg,data);

% If continuous, epoch into arbitrary 2-second segments
if length(data.trial) == 1
    cfg = []; cfg.event = 1:2*data.fsample:length(data.time{1}); cfg.epoch = [0 (2*data.fsample)-1];
    cfg.event(end) = [];
    data = ft_epoch(cfg,data);
end

cfg = []; cfg.keeptrials = 'yes';
tlock = ft_timelockanalysis(cfg,data);
tlock.cov = permute(repmat(opts.noisecov,1,1,length(data.trial)),[3 1 2]);

% Estimate sources
if strcmpi(opts.method,'eloreta')
    cfg = [];
    cfg.method = 'eloreta';
    cfg.grid = leadfield;
    cfg.eloreta.keepfilter = 'yes';
    cfg.eloreta.normalize = 'yes';
    cfg.headmodel = headmodel;
elseif strcmpi(opts.method,'mne')
    cfg = [];
    cfg.method = 'mne';
    cfg.grid = leadfield;
    cfg.headmodel = headmodel;
    cfg.mne.prewhiten = 'yes';
    cfg.mne.lambda    = 3;
    cfg.mne.keepfilter = 'yes';
    cfg.mne.scalesourcecov = 'yes';
end
sources = ft_sourceanalysis(cfg,tlock);

% Get voxel time courses
voxeldata = struct;
datacat = cat(2,data.trial{:});
for c = 1:size(sources.pos,1)
    tmp = sources.avg.filter{c}*datacat;
    u = svd(tmp,'econ');
    tmp = u(:,1)'*sources.avg.filter{c}*datacat;
    source_datacat(c,:) = tmp;
end
voxeldata.trial{1} = source_datacat;
voxeldata.time{1} = linspace(0,length(voxeldata.trial{1}/data.fsample),size(voxeldata.trial{1},2));
voxeldata.label = cellstr(num2str([1:size(sources.pos,1)]'));
voxeldata.fsample = data.fsample;
source_datacat = []; %saving memory

% Get ROI time courses
roidata = voxeldata;
roidata.label = atlas.(opts.atlaslabel);
roidata.trial = [];
roidata.trial{1} = NaN(length(roidata.label),length(voxeldata.trial{1}));
for cc = 1:length(roidata.label)
    roidata.trial{1}(cc,:) = mean(voxeldata.trial{1}(find(atlas.(opts.atlasparam) == cc),:),1);
end
%roidata.trial{1}(find(isnan(roidata.trial{1}(:,1))),:) = [];

end
