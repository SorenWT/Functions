function plot_tc_topo(time,datamat,settings,p,pindx,varargin)
% time should be in ms
% datamat should be in the form channels x time points x subjects x
% conditions
% settings is a struct with fields 'layout' (a Fieldtrip layout) and 'datatype' (EEG or MEG)
% p is a panel object
% pindx is the index of the subpanel you want to plot into, supplied as a
% cell array (i.e. if you want to plot into p(1,2), supply {1 2}). If you
% don't want to plot into a subpanel, put {}


if CheckInput(varargin,'avgfunc')
    avgfunc = EasyParse(varargin,'avgfunc');
else
    avgfunc = @nanmean;
end

if CheckInput(varargin,'topolocation')
    topoloc = EasyParse(varargin,'topolocation');
else
    topoloc = 'bot_outside';
end

% doesn't do anything right now, but later should add the ability to choose
% how many topoplots to use
if CheckInput(varargin,'ntopos')
    ntopos = EasyParse(varargin,'ntopos');
else
    ntopos = 4;
end

if nargin < 4
    p = panel('no-manage-font');
    pindx = {};
end

switch topoloc
    case 'bot_inside'
        p(pindx{:}).pack();
        for cc = 1:4
            p(pindx{:}).pack({[0.25*(cc-1) 0 0.25 0.15]})
        end
        mainindx = {1};
        topoindx = {{2} {3} {4} {5}};
    case 'bot_outside'
        p(pindx{:}).pack('v',{50 50})
        p(pindx{:},2).pack(2,2)
        mainindx = {1};
        topoindx = {{2 1 1},{2 1 2},{2 2 1},{2 2 2}};
    case 'right_outside'
        p(pindx{:}).pack('h',{60 40})
        p(pindx{:},2).pack('v',repmat({1/ntopos},1,ntopos))
        mainindx = {1};
        topoindx = {{2 1} {2 2} {2 3} {2 4}};
end


p(pindx{:},mainindx{:}).select()
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
    %axes(p(2,c,cc+1).axis)
    if size(datamat,4) > 1
        plotdata = avgfunc(squeeze(datamat(:,plotindx(cc),:,2))...
            - squeeze(datamat(:,plotindx(cc),:,1)),2);
    else
        plotdata = avgfunc(squeeze(datamat(:,plotindx(cc),:)),2);
    end
    if strcmpi(settings.datatype,'MEG')
        p(pindx{:},topoindx{cc}{:}).select()
        
        ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
            ones(size(datamat,1)),zeros(size(datamat,1)));
        p(pindx{:},topoindx{cc}{:}).title([num2str(time(plotindx(cc))) ' ms']);
        ax(cc) = gca;
    elseif strcmpi(settings.datatype,'EEG')
        p(pindx{:},topoindx{cc}{:}).select()
        
        cluster_topoplot(plotdata,settings.layout,...
            ones(size(datamat,1)),zeros(size(datamat,1)));
        t = p(pindx{:},topoindx{cc}{:}).title([num2str(time(plotindx(cc))) ' ms']);
        ax(cc) = gca;
    elseif strcmpi(settings.datatype,'source')
        p = ft_cluster_sourceplot(plotdata,settings.layout,settings.layout,ones(size(plotdata)),'method','wholebrain','panel',p,'panelindx',{pindx{:} topoindx{cc}{:}});
        t = p(pindx{:},topoindx{cc}{:},1,2).title([num2str(time(plotindx(cc))) ' ms']);
        ax(cc) = p(pindx{:},topoindx{cc}{:},1,1).axis;
                axis(ax(cc),'tight')
        set(findall(p(pindx{:},topoindx{cc}{:},1,2).axis,'type','text'),'visible','on')
        
    end
    if cc == 4
        cbar = colorbar(ax(cc),'EastOutside');
    end
    t.FontSize = 10; t.FontWeight = 'bold';
    %ax(cc) = gca;
    %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
end
Normalize_Clim(ax,1);
ax2 = p(pindx{:},1).axis;
cbar.Position = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) cbar.Position(3) 0.15*ax2.Position(4)];

