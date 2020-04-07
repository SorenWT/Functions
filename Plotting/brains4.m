function [pnl,clim,ax,cbar] = brains4(datain,sourcemodel,atlas,varargin)

if CheckInput(varargin,'mask')
    roimask = EasyParse(varargin,'mask');
end

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

if ~CheckInput(varargin,'mask')
   plotmask = ones(size(plotdata)); 
end

cort_size = length(find(sourcemodel.brainstructure == 1));

clr = [0.7804 0.7609 0.6627].*1.1;


if CheckInput(varargin,'panel')
    panelindx = EasyParse(varargin,'panelindx');
    pnl = EasyParse(varargin,'panel');
    pnl(panelindx{:}).pack(2,2)
    pnl(panelindx{:},1,1).select();
    ax(1) = gca;
else
    figure
    subplot(2,2,1)
end
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 1),:);
trindx = find(max(sourcemodel.tri,[],2) <= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);

ft_plot_mesh(bnd,'facealpha',vert(1-plotmask(sourcemodel.brainstructure==1)),'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==1),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==1)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==1))), 'maskstyle', 'opacity','edgealpha',0);
set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1))) max(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1)))]);

%ft_plot_mesh(bnd,'vertexcolor',plotdata(sourcemodel.brainstructure == 1)','edgealpha',0);
view(-180,0)
lighting gouraud 
camlight

cbar = colorbar('Location','EastOutside');
if CheckInput(varargin,'panel')
    pnl(panelindx{:},1,2).select();
    ax(2) = gca;
else
    subplot(2,2,2)
end
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 2),:);
trindx = find(min(sourcemodel.tri,[],2) >= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);
bnd.tri = bnd.tri - cort_size;
view(0,0)
ft_plot_mesh(bnd,'facealpha',vert(1-plotmask(sourcemodel.brainstructure==2)),'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==2),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==2)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==2))), 'maskstyle', 'opacity','edgealpha',0);
set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2))) max(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2)))]);

lighting gouraud 
camlight

if CheckInput(varargin,'panel')
    pnl(panelindx{:},2,1).select();
    ax(3) = gca;
else
    subplot(2,2,3)
end
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 1),:);
trindx = find(max(sourcemodel.tri,[],2) <= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);

ft_plot_mesh(bnd,'facealpha',vert(1-plotmask(sourcemodel.brainstructure==1)),'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==1),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==1)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==1))), 'maskstyle', 'opacity','edgealpha',0);
set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1))) max(vert(plotdata(sourcemodel.brainstructure==1)).*vert(plotmask(sourcemodel.brainstructure==1)))]);

view(0,0)

lighting gouraud 
camlight

if CheckInput(varargin,'panel')
    pnl(panelindx{:},2,2).select();
    ax(4) = gca;
else
    subplot(2,2,4)
end
bnd.pnt = sourcemodel.pos(find(sourcemodel.brainstructure == 2),:);
trindx = find(min(sourcemodel.tri,[],2) >= cort_size);
bnd.tri = sourcemodel.tri(trindx,:);
bnd.tri = bnd.tri - cort_size;

ft_plot_mesh(bnd,'facealpha',vert(1-plotmask(sourcemodel.brainstructure==2)),'edgealpha',0,'vertexcolor',repmat(clr,sum(sourcemodel.brainstructure==2),1));
hold on
ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata(sourcemodel.brainstructure==2)), 'facealpha', double(vert(plotmask(sourcemodel.brainstructure==2))), 'maskstyle', 'opacity','edgealpha',0);
set(gca,'CLim',[min(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2))) max(vert(plotdata(sourcemodel.brainstructure==2)).*vert(plotmask(sourcemodel.brainstructure==2)))]);

view(180,0)

lighting gouraud 
camlight

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
   colormap(EasyParse(varargin,'colormap')) 
end

if CheckInput(varargin,'frame')
    %replaces first output with frame - can't be used with panel-based plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % make figure as big as possible for high res image
    for i = 1:length(ax)
        pnl(i) = getframe(ax(i));
    end
    delete(gcf)
end
