function [p,ax] = tf_clustplot(stats,layout,whichcluster,plotdata,p,pindx,varargin)

argsin = varargin;

argsin = setdefault(argsin,'fsizemultiplier',1);
argsin = setdefault(argsin,'selmethod','peak');
argsin = setdefault(argsin,'sourceplotmethod','brains4');

selmethod = EasyParse(argsin,'selmethod');
fsizmult = EasyParse(argsin,'fsizemultiplier');
sourceplotmethod = EasyParse(argsin,'sourceplotmethod');



if nargin < 4 || isempty(plotdata)
    plotdata = stats.stat;
end
stats.plotdata = plotdata;


if ~isfield(stats,'negclusters')
    stats.negclusters = struct();
end

if ~isfield(stats,'posclusters')
    stats.posclusters = struct();
end

if nargin < 3
    whichcluster = -(length(stats.negclusters)):length(stats.posclusters);
    whichcluster(whichcluster==0) = [];
    whichcluster(abs(whichcluster)>5) = []; % don't plot more than 5 clusters of each type
end


if length(whichcluster)>1
    figure
    set(gcf,'units','normalized','position',[0 0 1 1])
    
    p = panel('no-manage-font');
    %p = panel();%'no-manage-font');
    %p.fontsize = 18*fsizmult;
    
    p.pack('h',{1/2 1/2})
    p(1).pack('v',repmat({1/max(abs(whichcluster))},1,max(abs(whichcluster))))
    p(2).pack('v',repmat({1/max(abs(whichcluster))},1,max(abs(whichcluster))))
    
    
    for i = 1:length(whichcluster)
        %if i < 0
        %    p(1,abs(i)).select()
        %else
        %    p(2,abs(i)).select()
        %end
        [p,ax{i}]=tf_clustplot(stats,layout,whichcluster(i),stats.plotdata,p,{double(whichcluster(i)>0)+1 abs(whichcluster(i))});
    end
    
    ax = cat(2,ax{:});
    
    Normalize_Clim(ax,1);
    
    p(2).marginleft = 50;
    p.margintop = 12;
    p.marginbottom = 25;
    p.marginright = 25;
    
    
else
    if ~exist('p','var')
        p = panel('no-manage-font');
    end
    if ~exist('pindx','var')
        pindx = {};
    end
    
    p(pindx{:}).pack('h',{1/2 1/2})
    p.fontsize = 18*fsizmult;
    p.fontweight = 'normal';
    
    %ax(1) = gca;
    
    i=1;
    
    % find the most representative time/frequency point
    
    if whichcluster(i) < 0
        clustlabs = reshape(stats.negclusterslabelmat,size(stats.negclusterslabelmat,1),[]);
    else
        clustlabs = reshape(stats.posclusterslabelmat,size(stats.posclusterslabelmat,1),[]);
    end
    
    maps = reshape(stats.plotdata,size(stats.plotdata,1),[]);
    maps = maps(:,any(clustlabs==abs(whichcluster(i))));
    
    reprmask = clustlabs(:,any(clustlabs==abs(whichcluster(i))));
    
    
    switch selmethod
        case 'representative'
            [~,~,~,d] = kmeans(maps',1);
            choiceindx = find(d == min(d));
        case 'peak'
            tmp = find(abs(maps)==max(max(abs(maps.*(reprmask==abs(whichcluster(i)))))));
            [~,choiceindx] = ind2sub(size(maps),tmp);
    end
    
    
    
    reprmap = maps(:,choiceindx);
    
    reprmask = reprmask(:,choiceindx);
    reprmask = reprmask==abs(whichcluster(i));
    
    tmp = find(any(clustlabs==abs(whichcluster(i))));
    indxpoint = tmp(choiceindx);
    
    %     tmp = zeros(1,size(clustlabs,2));
    %     tmp(indxpoint) = 1;
    %     tfindxmat = reshape(tmp,size(stats.stat,2),size(stats.stat,3));
    %     indxpoint = find(tfindxmat);
    
    [findx,tindx] = ind2sub([size(stats.stat,2) size(stats.stat,3)],indxpoint);
    
    reprtfpoint = {[num2str(round(stats.freq(findx),3,'sig')) ' Hz, ' num2str(stats.time(tindx)) ' seconds']};
    
    if ~contains(ft_datatype(layout),'mesh')
        %p(pindx{:},1).select()
        f = figure;
        ft_cluster_topoplot(layout,reprmap,stats.label,stats.prob,reprmask);
        p(pindx{:},1).select(gca); close(f);
        ax{1} = gca;
    else
        if strcmpi(sourceplotmethod,'wholebrain')
            %p(pindx{:},1).select()
            ft_cluster_sourceplot(reprmap,layout,layout,reprmask,'panel',p,'panelindx',{pindx{:} 1},'method','wholebrain')
            ax{1} = p(pindx{:},1,1).axis;
        else
            [p,~,ax{1},cbar] = brains4(reprmap,layout,layout,'mask',reprmask,'panel',p,'panelindx',{pindx{:} 1});
            delete(cbar)
        end
    end
    t = p(pindx{:},1).title(reprtfpoint);%,'FontSize',18*fsizmult)
    t.FontWeight = 'bold'; t.FontSize = 18*fsizmult; t.Units = 'inches';
    
    % find the most representative sensor/region
    
    %p(pindx{:},2).select()
    ax{2} = gca;
    
    if whichcluster(i) < 0
        clustlabs = reshape(stats.negclusterslabelmat,size(stats.negclusterslabelmat,1),[]);
    else
        clustlabs = reshape(stats.posclusterslabelmat,size(stats.posclusterslabelmat,1),[]);
    end
    
    tf = reshape(stats.plotdata,size(stats.plotdata,1),[]);
    tf = tf(any(clustlabs==abs(whichcluster(i)),2),:);
    
    chanlab = stats.label(any(clustlabs==abs(whichcluster(i)),2));
    clustlabs = clustlabs(any(clustlabs==abs(whichcluster(i)),2),:);
    nonantf = tf;
    nonantf(:,any(isnan(tf),1)) = [];
    
    reprmask = clustlabs;
    %reprmask = clustlabs(any(clustlabs==abs(whichcluster(i))),:);
    
    
    switch selmethod
        case 'representative'
            [~,~,~,d] = kmeans(nonantf,1);
            choiceindx = find(d==min(d));
        case 'peak'
            tmp = find(abs(tf)==max(max(abs(tf.*(reprmask==abs(whichcluster(i)))))));
            [choiceindx,~] = ind2sub(size(tf),tmp);
    end
    
    
    reprtf = tf(choiceindx,:);
    reprtf = reshape(reprtf,size(stats.stat,2),size(stats.stat,3));
    
    reprmask = clustlabs(choiceindx,:);
    reprmask = reprmask==abs(whichcluster(i));
    %reprmask = clustlabs(choiceindx,:) == abs(whichcluster(i));
    reprmask = reshape(reprmask,size(stats.stat,2),size(stats.stat,3));
    
    reprchan = chanlab(choiceindx);
    f = figure;
    easy_freqplot(struct('toplot',reprtf,'time',stats.time,'freq',stats.freq,'mask',reprmask))
    p(pindx{:},2).select(gca); close(f)
    
    FixAxes(gca,16*fsizmult)
    set(gca,'YScale','log')
    if whichcluster < 0
        thisp = stats.negclusters(abs(whichcluster(i))).prob;
    else
        thisp = stats.posclusters(abs(whichcluster(i))).prob;
        
    end
    t = p(pindx{:},2).title([reprchan ': p = ' num2str(round(thisp,2,'significant'))]);
    t.FontWeight = 'bold'; t.FontSize = 18*fsizmult;
    
    
    p(pindx{:}).de.marginleft = 30;
    p(pindx{:}).marginbottom = 40;
    
    ax = cat(2,ax{:});
    
    
end

