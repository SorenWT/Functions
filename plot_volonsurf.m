function [pnl,clim,ax,cbar] = plot_volonsurf(vect,vatlas,info,surf,varargin)

% to do: fix bad spots with interpolation

V = vect2vol(vect,vatlas);

argsin = varargin;

niftiwrite(V,'tmp.nii',info)

voldat = ft_read_mri('tmp.nii');

if isstr(surf)
    if strcmpi(surf,'conte69')
        load('conte69_fs_LR_8k.mat')
        surf = conte69;
        load('Conte69_inflated_8k.mat')
        infl = 1;
    else
       load(surf) 
       infl = 0;
    end
else
    infl = 0;
end


voldat.inside = voldat.anatomy~=0;
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'anatomy';
plotmodel = ft_sourceinterpolate(cfg, voldat, surf);
plotmodel.brainstructure = surf.brainstructure;
zeroindx = find(plotmodel.anatomy==0); replacezero = zeros(size(zeroindx));
for i = 1:length(zeroindx)
    replacezero(i) = neighb_replacezero(zeroindx(i),plotmodel.tri,plotmodel.anatomy,0);
end
plotmodel.anatomy(zeroindx) = replacezero;

if CheckInput(argsin,'mask')
    mask = EasyParse(argsin,'mask');
    for i = 1:length(plotmodel.anatomy)
        masksurf(i) = mask(find(abs((vect-plotmodel.anatomy(i))/plotmodel.anatomy(i)) < 0.001,1));
    end
    indx = find(strcmpi(argsin,'mask'));
    argsin(indx+1) = {masksurf};
end


if infl
    [pnl,clim,ax,cbar] = brains4(plotmodel.anatomy,inflated,[],argsin{:});%,'mask',plotmodel.anatomy~=0);
else
    [pnl,clim,ax,cbar] = brains4(plotmodel.anatomy,surf,[],argsin{:});%,'mask',plotmodel.anatomy~=0);
end

system('rm tmp.nii')
end

function replacezero = neighb_replacezero(indx,tri,vals,returnnan)
neighbs = tri(any((tri-indx) == 0,2),:);
neighbs = unique(neighbs); neighbs(neighbs==indx) = []; 
neighbvals = vals(neighbs); neighbvals(neighbvals==0) = [];
replacezero = mode(neighbvals);
if isnan(replacezero) && ~returnnan
    for i = 1:length(neighbs) % get the neighbors of the neighbors if all neighbors are zero initially 
        rz(i) = neighb_replacezero(neighbs(i),tri,vals,1);
    end
    replacezero = mode(rz(~isnan(rz)));
end
end