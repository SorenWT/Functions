function plot_erp_clusts(clust,datamat,settings,p,pindx,varargin)
% specialized version of this function for plotting ERP stuff, specifically
% with cluster-based permutation tests
% time should be in ms
% datamat should be in the form channels x time points x subjects x
% conditions
% settings is a struct with fields:
%     'layout' (a Fieldtrip layout or EEGLAB chanlocs structure)
%     'datatype' (EEG or MEG)
%     'nclusts' (the number of cluster pairs to plot - assumes that
%     clusters come in pairs, as for dipolar topographies). Defaults to 1
%     'condnames' (the names of the conditions, to put in the legend.
%     Defaults to {'Condition 1','Condition 2'}
% p is a panel object (not required)
% pindx is the index of the subpanel you want to plot into, supplied as a
% cell array (i.e. if you want to plot into p(1,2), supply {1 2}). If you
% don't want to plot into a subpanel, put {}


set(gcf,'color','w')
settings = setdefault(settings,'nclusts',1);
settings = setdefault(settings,'condnames',{'Condition 1','Condition 2'});

nclusts = settings.nclusts;

if CheckInput(varargin,'avgfunc')
    avgfunc = EasyParse(varargin,'avgfunc');
else
    avgfunc = @nanmean;
end

if CheckInput(varargin,'shading')
    shadeopts = EasyParse(varargin,'shading');
else
    shadeopts = 'sem';
end

if nargin < 5
    p = panel('no-manage-font');
    pindx = {};
end

p(pindx{:}).pack('v',repmat({1/nclusts},1,nclusts));
for i = 1:nclusts
    p(pindx{:},i).pack('h',{1/3 1/3 1/3})
    p(pindx{:},i,2).pack(3,3)
end
% p(pindx{:}).pack();
% for cc = 1:4
%     p(pindx{:}).pack({[0.25*(cc-1) 0 0.25 0.15]})
% end
% p(pindx{:},1).select()

%hold on

if CheckInput(varargin,'color')
    clr = EasyParse(varargin,'color');
else
    clr = {'b' 'r'};
end

for i = 1:nclusts
    if isfield(clust,'posclusterslabelmat')
        thesechans_pos = find(sum(clust.posclusterslabelmat==i,2)>0.25*max(sum(clust.posclusterslabelmat==i,2)));
    else
        thesechans_pos = NaN;
    end
    if isfield(clust,'negclusterslabelmat')
        thesechans_neg = find(sum(clust.negclusterslabelmat==1,2)>0.25*max(sum(clust.negclusterslabelmat==1,2)));
    else
        thesechans_neg = NaN;
    end
    
    if ~all(isnan(thesechans_pos))
        p(pindx{:},i,1).select()
        stdshade(clust.time,squeeze(avgfunc(datamat(thesechans_pos,:,:,1),1)),clr{1},0.15,1,shadeopts);
        hold on
        
        if size(datamat,4) > 1
            stdshade(clust.time,squeeze(avgfunc(datamat(thesechans_pos,:,:,2),1)),clr{2},0.15,1,shadeopts);
        end
        
        xlabel('Time (s)')
        if strcmpi(settings.datatype,'EEG')
            ylabel('Voltage (\muV)')
        elseif strcmpi(settings.datatype,'MEG')
            ylabel('Field strength (T)');
        end
        
        FixAxes(gca)
        set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
        title(['Positive cluster ' num2str(i)])
        manlegend(settings.condnames,clr)
    end
    
    plotindx = linspace(0,length(clust.time),10);
    plotindx = movmean(plotindx,2);
    plotindx(1) = [];
    plotindx = round(plotindx);
    
    %plothelper = [1 1; 1 2; 2 1; 2 2];
    plothelper = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];
    
    for cc = 1:9
        p(pindx{:},i,2,plothelper(cc,1),plothelper(cc,2)).select()
        %axes(p(2,c,cc+1).axis)
        if size(datamat,4) > 1
            plotdata = avgfunc(squeeze(datamat(:,plotindx(cc),:,2))...
                - squeeze(datamat(:,plotindx(cc),:,1)),2);
        else
            plotdata = avgfunc(squeeze(datamat(:,plotindx(cc),:)),2);
        end
        
        topochans = ones(size(clust.mask,1),1);
        if isfield(clust,'posclusterslabelmat')
            topochans(clust.posclusterslabelmat(:,plotindx(cc))==i) = 0;
        end
        if isfield(clust,'negclusterslabelmat')
            topochans(clust.negclusterslabelmat(:,plotindx(cc))==i) = 0;
        end
        
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                topochans,clust.mask(:,plotindx(cc)).*topochans);
        else
            cluster_topoplot(plotdata,settings.layout,...
                topochans,clust.mask(:,plotindx(cc)).*topochans);
        end
        if cc == 4
            cbar = colorbar('EastOutside');
            if strcmpi(settings.datatype,'EEG')
                cbar.Label.String = 'Voltage (\muV)';
            elseif strcmpi(settings.datatype,'MEG')
                cbar.Label.String = 'Field strength (T)';
            end
        end
        title([num2str(clust.time(plotindx(cc))) ' ms'],'FontSize',10)
        ax(cc) = gca;
        %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
    Normalize_Clim(ax,1);
    %ax2 = p(pindx{:},1).axis;
    %cbar.Position = [ax2.Position(1)+ax2.Position(3) ax2.Position(2) cbar.Position(3) 0.15*ax2.Position(4)];
    
    if ~all(isnan(thesechans_neg))
        p(pindx{:},i,3).select()
        stdshade(clust.time,squeeze(avgfunc(datamat(thesechans_neg,:,:,1),1)),clr{1},0.15,1,shadeopts);
        hold on
        
        if size(datamat,4) > 1
            stdshade(clust.time,squeeze(avgfunc(datamat(thesechans_neg,:,:,2),1)),clr{2},0.15,1,shadeopts);
        end
        
        xlabel('Time (s)')
        if strcmpi(settings.datatype,'EEG')
            ylabel('Voltage (\muV)')
        elseif strcmpi(settings.datatype,'MEG')
            ylabel('Field strength (T)');
        end
        
        FixAxes(gca)
        set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
        title(['Negative cluster ' num2str(i)])
        manlegend(settings.condnames,clr)
    end
end
