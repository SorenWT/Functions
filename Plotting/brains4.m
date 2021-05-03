function [pnl,clim,ax,cbar] = brains4(datain,sourcemodel,atlas,varargin)

if ~exist('atlas','var')
    atlas = [];
end

if CheckInput(varargin,'mask')
    roimask = EasyParse(varargin,'mask');
end

if ~CheckInput(varargin,'material')
    mtrl = [0.45 0.55 0 20 1];
else
    mtrl = EasyParse(varargin,'material');
end

if ~isempty(atlas)
    for c = 1:length(sourcemodel.pos)
        if isfield(atlas,'parcellation')
            plotdata(c) = datain(atlas.parcellation(c));
        elseif isfield(atlas,'parcels')
            plotdata(c) = datain(atlas.parcels(c));
        elseif isfield(atlas,'tissue')
            if atlas.tissue(c) > 0
                plotdata(c) = datain(atlas.tissue(c));
            else
                plotdata(c) = NaN;
            end
        end
        
        
        if CheckInput(varargin,'mask')
            if isfield(atlas,'parcellation')
                plotmask(c) = roimask(atlas.parcellation(c));
            else
                plotmask(c) = roimask(atlas.parcels(c));
            end
        end
        
    end
else
    plotdata = datain;
    if CheckInput(varargin,'mask')
        plotmask = roimask;
    end
end

if ~CheckInput(varargin,'mask')
    plotmask = ones(size(plotdata));
end

cort_size = length(find(sourcemodel.brainstructure == 1));

%clr = [0.87 0.87 0.87];
clr = [0.7804 0.7609 0.6627].*1.1;


if CheckInput(varargin,'panel')
    panelindx = EasyParse(varargin,'panelindx');
    pnl = EasyParse(varargin,'panel');
    pnl(panelindx{:}).pack(2,2)
    pnl(panelindx{:}).de.margin = [3 3 3 3];
    pnl(panelindx{:},1,1).select();
else
    subplot(2,2,1)
end
ax(1) = gca;
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 1),:);
trindx = find(max(sourcemodel.tri,[],2) <= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);

ft_plot_mesh(bnd,'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==1),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==1)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==1))), 'maskstyle', 'opacity','edgealpha',0);
material(mtrl)
try
    set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1))) max(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1)))]);
catch
end


%ft_plot_mesh(bnd,'vertexcolor',plotdata(sourcemodel.brainstructure == 1)','edgealpha',0);
view(-90,0)
lighting gouraud
camlight
axis tight

cbar = colorbar('Location','EastOutside');
if CheckInput(varargin,'panel')
    pnl(panelindx{:},1,2).select();
else
    subplot(2,2,2)
end
ax(2) = gca;
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 2),:);
trindx = find(min(sourcemodel.tri,[],2) >= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);
bnd.tri = bnd.tri - cort_size;
view(90,0)
ft_plot_mesh(bnd,'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==2),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==2)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==2))), 'maskstyle', 'opacity','edgealpha',0);
material(mtrl)
try
    set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2))) max(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2)))]);
catch
end

lighting gouraud
camlight
axis tight

if CheckInput(varargin,'panel')
    pnl(panelindx{:},2,1).select();
else
    subplot(2,2,3)
end
ax(3) = gca;
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 1),:);
trindx = find(max(sourcemodel.tri,[],2) <= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);

ft_plot_mesh(bnd,'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==1),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==1)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==1))), 'maskstyle', 'opacity','edgealpha',0);
material(mtrl)
try
    set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1))) max(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1)))]);
catch
end

view(90,0)

lighting gouraud
camlight
axis tight

if CheckInput(varargin,'panel')
    pnl(panelindx{:},2,2).select();
else
    subplot(2,2,4)
end
ax(4) = gca;
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 2),:);
trindx = find(min(sourcemodel.tri,[],2) >= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);
bnd.tri = bnd.tri - cort_size;

ft_plot_mesh(bnd,'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==2),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==2)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==2))), 'maskstyle', 'opacity','edgealpha',0);
material(mtrl)
try
    set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2))) max(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2)))]);
catch
end

view(-90,0)

lighting gouraud
camlight
axis tight

set(gcf,'color','w')

if ~CheckInput(varargin,'panel')
    clim = Normalize_Clim(gcf);
    ax = findobj('parent',gcf,'type','axes');
    ax = ax(end:-1:1);
else
    clim = Normalize_Clim(ax);
end

if CheckInput(varargin,'CLim')
    Set_Clim(ax,EasyParse(varargin,'CLim'))
end

if CheckInput(varargin,'colormap')
    for i = 1:length(ax)
        colormap(ax(i),EasyParse(varargin,'colormap'))
    end
end

if CheckInput(varargin,'frame')
    %replaces first output with frame - can't be used with panel-based plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % make figure as big as possible for high res image
    for i = 1:length(ax)
        pnl(i) = getframe(ax(i));
    end
    delete(gcf)
end

if ~exist('pnl','var')
    pnl = [];
end

