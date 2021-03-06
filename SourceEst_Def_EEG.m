function [roidata,voxeldata,sources,sourcemodel] = SourceEst_Def_EEG(EEG,sourcemodel,atlas,parflag)

[~,ftpath] = ft_version;
%ftpath = '/group/northoff/share/fieldtrip-master';

% Convert the data to Fieldtrip format
data = eeglab2fieldtrip(EEG,'preprocessing','none');

% Load templates
load('standard_bem.mat') %variable is "vol"
if ~exist('sourcemodel','var') || isempty(sourcemodel)
	sourcemodel = parload('conte69_fs_LR_8k.mat','conte69');
%    sourcemodel = ft_read_headshape('conte69_fs_LR_8k.mat');
end
if ~exist('atlas','var') || isempty(atlas)
    load('brod_balsa_fs_LR_8k.mat');
    atlas = brod;
end
%sourcemodel = ft_read_headshape('cortex_5124.surf.gii');
%atlas = ft_read_atlas(fullfile(ftpath,'template','atlas','aal','ROI_MNI_V4.nii'));

% % Interpolate template surface on atlas
% cfg = [];
% cfg.interpmethod = 'nearest';
% cfg.parameter = 'tissue';
% sourcemodel2 = ft_sourceinterpolate(cfg, atlas, sourcemodel);

% Align electrodes to scalp surface
data.elec = ft_convert_units(data.elec,'mm');
sourcemodel = ft_convert_units(sourcemodel,'mm');
vol = ft_convert_units(vol,'mm');
cfg = []; cfg.method = 'project'; cfg.headshape = vol.bnd(1);
elec = ft_electroderealign(cfg,data.elec);
data.elec = elec;

% Create forward model
cfg = []; cfg.elec = data.elec; cfg.channel = {'EEG'};
cfg.grid.pos = sourcemodel.pos; cfg.grid.inside = 1:size(sourcemodel.pos,1);
cfg.headmodel = vol;
leadfield = ft_prepare_leadfield(cfg);

% Select only the EEG channels
cfg = []; cfg.channel = {'EEG'};
eegdata = ft_selectdata(cfg,data);
clear data

% Epoch into arbitrary 2-second segments
cfg = []; cfg.event = 1:2*eegdata.fsample:length(eegdata.time{1}); cfg.epoch = [0 (2*eegdata.fsample)-1];
cfg.event(end) = [];
eegdata = ft_epoch(cfg,eegdata);

% Estimate sources
cfg = [];
cfg.method = 'eloreta';
cfg.grid = leadfield;
cfg.eloreta.keepfilter = 'yes';
cfg.eloreta.normalize = 'yes';
cfg.headmodel = vol;
sources = ft_sourceanalysis(cfg,eegdata);

% Get voxel time courses
roidata = struct;
voxeldata = struct;
datacat = cat(2,eegdata.trial{:});
source_datacat = zeros(size(sources.pos,1),size(datacat,2));
if parflag
    parfor c = 1:size(sources.pos,1)
        tmp = sources.avg.filter{c}*datacat;
        u = svd(tmp,'econ');
        source_datacat(c,:) = u(:,1)'*sources.avg.filter{c}*datacat;
    end
else
    for c = 1:size(sources.pos,1)
        tmp = sources.avg.filter{c}*datacat;
        u = svd(tmp,'econ');
        tmp = u(:,1)'*sources.avg.filter{c}*datacat;
        source_datacat(c,:) = u(:,1)'*sources.avg.filter{c}*datacat;
    end
end
voxeldata.trial{1} = source_datacat;
voxeldata.time{1} = linspace(0,length(voxeldata.trial{1})/eegdata.fsample-1/eegdata.fsample,size(voxeldata.trial{1},2));
voxeldata.label = cellstr(num2str([1:size(sources.pos,1)]'));
voxeldata.fsample = eegdata.fsample;
source_datacat = []; %saving memory

% Get ROI time courses
roidata = voxeldata;
roidata.label = atlas.BAlabel;
roidata.trial = [];
roidata.trial{1} = NaN(length(roidata.label),length(voxeldata.trial{1}));
for cc = 1:length(roidata.label)
    roidata.trial{1}(cc,:) = mean(voxeldata.trial{1}(find(atlas.BA == cc),:),1);
end
%roidata.trial{1}(find(isnan(roidata.trial{1}(:,1))),:) = [];

end
