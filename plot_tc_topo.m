function plot_tc_topo(time,datamat,settings,p,pindx,varargin)
% datamat should be in the form channels x time points x subjects x
% conditions
% time should be in ms
%settings only needs to have a layout and the field 'datatype' (EEG or MEG)

if CheckInput(varargin,'avgfunc')
    avgfunc = EasyParse(varargin,'avgfunc');
else
    avgfunc = @nanmean;
end

p(pindx{:}).pack();
for cc = 1:4
    p(pindx{:}).pack({[0.25*(cc-1) 0 0.25 0.15]})
end
p(pindx{:},1).select()

hold on

if CheckInput(varargin,'color')
    clr = EasyParse(varargin,'color');
else
    clr = {'b' 'r'};
end

stdshade(time,squeeze(avgfunc(datamat(:,:,:,1),1)),clr{1},0.15,1,'sem');

if size(datamat,4) > 1
    stdshade(time,squeeze(avgfunc(datamat(:,:,:,2),1)),clr{2},0.15,1,'sem');
end


xlabel('Time (s)')

FixAxes(gca)
set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)

plotindx = linspace(0,length(time),5);
plotindx = movmean(plotindx,2);
plotindx(1) = [];
plotindx = round(plotindx);

%tindx =
for cc = 1:4
    p(pindx{:},cc+1).select()
    %axes(p(2,c,cc+1).axis)
    if size(datamat,4) > 1
        plotdata = avgfunc(squeeze(datamat(:,plotindx(cc),:,2))...
            - squeeze(datamat(:,plotindx(cc),:,1)),2);
    else
        plotdata = avgfunc(squeeze(datamat(:,plotindx(cc),:)),2);
    end
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
            ones(size(datamat,1)),zeros(size(datamat,1)));
    else
        cluster_topoplot(plotdata,settings.layout,...
            ones(size(datamat,1)),zeros(size(datamat,1)));
    end
    if cc == 4
        cbar = colorbar('EastOutside');
    end
    title([num2str(time(plotindx(cc))) ' ms'],'FontSize',10)
    ax(cc) = gca;
    %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
end
Normalize_Clim(ax,1);
ax2 = p(pindx{:},1).axis;
cbar.Position = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) cbar.Position(3) 0.15*ax2.Position(4)];

