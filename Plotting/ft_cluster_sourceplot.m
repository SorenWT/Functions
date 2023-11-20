function [p,clim,ax,cbar] = ft_cluster_sourceplot(datain,sourcemodel,atlas,roimask,varargin)

for c = 1:length(sourcemodel.pos)
    if isfield(atlas,'parcellation')
        plotdata(c) = datain(atlas.parcellation(c));
        plotmask(c) = roimask(atlas.parcellation(c));
    elseif isfield(atlas,'parcels')
        plotdata(c) = datain(atlas.parcels(c));
        plotmask(c) = roimask(atlas.parcels(c));
    end
end
plotmask = double(plotmask);

argsin = varargin;

argsin = setdefault(argsin,'method','brains4');

if EasyParse(argsin,'method','wholebrain')
    bnd.pnt = sourcemodel.pos;
    bnd.tri = sourcemodel.tri;
    
    mtrl = [0.45 0.55 0 20 1];
    
    if ~CheckInput(argsin,'panel')
        p = panel('no-manage-font');
        pindx = {};
    else
        p = EasyParse(argsin,'panel');
        pindx = EasyParse(argsin,'panelindx');
    end
    
    
    p(pindx{:}).pack()
    p(pindx{:},1).pack({[0 0 1 1]}) %make 2 separate panels so you can have a title that doesn't rotate
    p(pindx{:},1).pack({[0 0 1 1]})
    
    
    p(pindx{:},1,1).select()
    %p(pindx{:}).select()
    ft_plot_mesh(bnd,'facealpha',vert(1-plotmask));
    ft_plot_mesh(bnd, 'vertexcolor', vert(plotdata), 'facealpha', vert(plotmask), 'maskstyle', 'opacity');
    
    lighting gouraud
    %camlight
    material(mtrl)
    
    % automatically set the view to look at the most interesting part
    % (highest value)
    
    peakcoords = bnd.pnt(find(abs(plotdata)==max(abs(plotdata.*plotmask))),:);
    peakcoords = nanmean(peakcoords,1);
    
    view(peakcoords)
    camlight
    
    p(pindx{:},1,2).select()
    set(p(pindx{:},1,2).axis,'visible','off')
    set(findall(p(pindx{:},1,2).axis,'type','text'),'visible','on')
    %p(pindx{:},1,2).title('test')
    
    %a = p(pindx{:},1,2).axis;
    %a.Visible = 'off';
    
elseif EasyParse(argsin,'method','brains4')
    [p,clim,ax,cbar] = brains4(datain,sourcemodel,atlas,'mask',roimask,argsin{:});
end





