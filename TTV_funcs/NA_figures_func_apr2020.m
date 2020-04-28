function NA_figures_func(settings)

fbands = settings.tfparams.fbandnames;
load('lkcmap2.mat')

load([settings.outputdir '/' settings.datasetname '_allmeas.mat'])
if exist([settings.outputdir '/' settings.datasetname '_results.mat'],'file')
    load([settings.outputdir '/' settings.datasetname '_results_FDR.mat'])
    if isfield(settings,'rest')
        load([settings.outputdir '/' settings.datasetname '_restmeas_FDR.mat'])
    end
else
    load([settings.outputdir '/' settings.datasetname '_results.mat'])
    if isfield(settings,'rest')
        load([settings.outputdir '/' settings.datasetname '_restmeas.mat'])
    end
end

if strcmpi(settings.datatype,'EEG')
    %eeglab
end

cd([settings.outputdir '/' settings.datasetname '_figures'])

if strcmpi(settings.tfparams.method,'hilbert') || ~isempty(find(contains(settings.steps,'tf_filter')))
    prestim_pseudo = settings.pseudo.prestim;
    prestim_real = settings.real.prestim;
    poststim_pseudo = settings.pseudo.poststim;
    poststim_real = settings.real.poststim;
else
    prestim_pseudo = settings.pseudo.prestim - settings.pseudo.prestim(1)+1+settings.srate/5;
    prestim_real = settings.real.prestim - settings.pseudo.prestim(1)+1+settings.srate/5;
    poststim_pseudo = settings.pseudo.poststim - settings.pseudo.prestim(1)+1+settings.srate/5;
    poststim_real = settings.real.poststim - settings.pseudo.prestim(1)+1+settings.srate/5;
end


%% Figure 2a: Nonadditivity of ERSP in different frequency bands

figure

p = panel('no-manage-font');

set(gcf,'units','normalized','position',[0 0.5 1 0.5],'Color','w');

p.pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
for c = 1:settings.nfreqs
    p(c).pack();
    for cc = 1:4
        p(c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(c,1).select();    
    
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    hold on
    stdshade(t,squeeze(nanmean(allmeas{c}.naddersp.diff(:,:,1,:),1)),'b',0.15,1,'sem');
    
    stdshade(t,squeeze(nanmean(allmeas{c}.naddersp.diff(:,:,2,:),1)),'r',0.15,1,'sem');
    title(fbands{c})
    if c == settings.nfreqs
        legend({'Corrected prestim low','Corrected prestim high'})
    end
    xlabel('Time (s)')
    if c == 1
        ylabelunits(settings)
    end
    FixAxes(gca)
    set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
    
    plotindx = linspace(0,max(settings.aucindex),5);
    plotindx = movmean(plotindx,2);
    plotindx(1) = [];
    plotindx = round(plotindx);

    %tindx =
    for cc = 1:4
        p(c,cc+1).select()
        %axes(p(2,c,cc+1).axis)
        plotdata = nanmean(squeeze(allmeas{c}.naddersp.diff(:,plotindx(cc),2,:))...
            - squeeze(allmeas{c}.naddersp.diff(:,plotindx(cc),1,:)),2);
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                ones(size(alloutputs.ersp.pt.stats{c}.mask(:,plotindx(cc)))),0.*alloutputs.ersp.pt.stats{c}.mask(:,plotindx(cc)));
        else
            cluster_topoplot(plotdata,settings.layout,...
                ~(0.*alloutputs.ersp.pt.stats{c}.mask(:,plotindx(cc))),0.*alloutputs.ersp.pt.stats{c}.mask(:,plotindx(cc)));
        end
        colormap(lkcmap2)
        if cc == 4
            cbars(c) = colorbar;
        end
        ax(cc) = p(c,cc+1).axis;
        title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
        %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
    Normalize_Clim(ax,1);
    
end


p.de.margin = [5 5 5 5];
p.de.marginright = 20;
p.marginright = 10;
% fix margins here

set(gcf,'Color','w')

for c = 1:settings.nfreqs
    ax2(c) = p(c,1).axis;
    cbars(c).Position = [ax2(c).Position(1)+ax2(c).Position(3) ax2(c).Position(2) cbars(c).Position(3) 0.15*ax2(c).Position(4)];
end
Normalize_Ylim(ax2)

for c = 1:settings.nfreqs
    p(c,1).select()
    %Plot_sigmask(p(2,c,1).axis,alloutputs.ersp.pt.stats{c}.prob < 0.05,'cmapline','LineWidth',5)
    H = sigstar({get(p(c,1).axis,'XLim')},2*min(min(alloutputs.ersp.pt.stats{c}.prob)),0,18)
    delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
end


savefig('Fig2a.fig')
export_fig('Fig2a.png','-m5')
save('Panel2a.mat','p')
close


%% Figure 2b: TTV of ERSP

figure

p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'units','normalized','position',[0 0.5 1 0.5],'Color','w');

p.pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
for c = 1:settings.nfreqs
    p(c).pack();
    for cc = 1:4
        p(c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(c,1).select();
    plotband = c;
    
    t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
    %t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
    %prestimdata = 100*(allmeas{plotband}.raw.ttversp(:,prestim_real,:)-nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2))./...
    %    nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
    %prestimdata = nanmean(nanmean(prestimdata,3),1);
    poststimdata = (allmeas{plotband}.ttversp.real);
    poststimdata = nanmean(nanmean(poststimdata,3),1);
    %plotdata = [prestimdata poststimdata];
    plotdata = poststimdata;
    %FillBetween(t((length(settings.real.prestim)+1):end),poststimdata,...
    %    zeros(1,length(poststimdata)));
    hold on
    stdshade(t,squeeze(nanmean(allmeas{plotband}.ttversp.real,1)),'b',0.15,1,'sem')
    %plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
    
    title(fbands{c})
    xlabel('Time (s)')
    if c == 1
    ylabel('% change of TTV of ERSP')
    end
    %ylim = get(gca,'YLim');
    %line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
    %set(gca,'YLim',ylim)
    FixAxes(gca,14)
    set(gca,'XLim',[0 max(t)])
    %set(gca,'FontSize',16)
    
    plotindx = linspace(0,max(settings.aucindex),5);
    plotindx = movmean(plotindx,2);
    plotindx = round(plotindx);
    plotindx(1) = [];
    for cc = 1:4
        p(c,cc+1).select()
        plotdata = mean(allmeas{c}.ttversp.real(:,plotindx(cc),:),3);
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                ones(size(alloutputs.ersp.ttv.stats{c}.mask(:,plotindx(cc)))),0.*alloutputs.ersp.ttv.stats{c}.mask(:,plotindx(cc)));
        else
            cluster_topoplot(plotdata,settings.layout,...
                1-(0.*alloutputs.ersp.ttv.stats{c}.mask(:,plotindx(cc))),0.*alloutputs.ersp.ttv.stats{c}.mask(:,plotindx(cc)));
        end
        colormap(lkcmap2)
        if cc == 4
            cbars1(c) = colorbar;
        end
        ax(cc) = p(c,cc+1).axis;
        title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
        %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
    Normalize_Clim(ax,1);
end
clear ax

for c = 1:settings.nfreqs
    ax1(c) = p(c,1).axis;
    cbars1(c).Position = [ax1(c).Position(1)+ax1(c).Position(3) ax1(c).Position(2) cbars1(c).Position(3) 0.15*ax1(c).Position(4)];
end
Normalize_Ylim(ax1)
%Set_Ylim(ax1,[-80 60])

for c = 1:settings.nfreqs    
    p(c,1).select()
    %Plot_sigmask(p(c,1).axis,alloutputs.ersp.ttv.stats{c}.prob < 0.05,'cmapline','LineWidth',5)
    H = sigstar({get(p(c,1).axis,'XLim')},2*min(min(alloutputs.ersp.ttv.stats{c}.prob)),0,18)
    delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
end


p.de.margin = [5 5 5 5];
p.marginleft = 18;
p.marginbottom = 18;
p.margintop = 8;
p.marginright = 10;
p.de.marginleft = 18;

set(gcf,'Color','w')

savefig('Fig2b.fig')
export_fig('Fig2b.png','-m5')
save('Panel2b.mat','p')

%% Figure 2c: Method correlation

figure

p = panel('no-manage-font')

pos = get(gcf,'position');
set(gcf,'units','normalized','position',[0 0.5 1 0.5],'Color','w');

p.pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
alloutputs.ersp.corr.r = real(alloutputs.ersp.corr.r);
for c = 1:settings.nfreqs
    p(c).pack()
    p(c).pack({[0 0.7 0.5 0.3]})
    if ~isempty(find(~isnan(alloutputs.ersp.corr.r(:,c)))) && ~isempty(find(alloutputs.ersp.corr.r(:,c)))
    p(c,1).select()
    nicecorrplot(nanmean(allmeas{c}.naerspindex,1),nanmean(allmeas{c}.ttverspindex,1),...
        {['Pseudotrial-based' newline 'ERSP nonadditivity'],'TTV-based ERSP nonadditivity'},'Plot','r');
    FixAxes(gca,14)
    
    title(fbands{c})
    p(c,2).select()
    if strcmpi(settings.datatype,'MEG')
        if isempty(find(~isnan(alloutputs.ersp.corr.r(:,c))))
            alloutputs.ersp.corr.r(:,c) = zeros(size(alloutputs.ersp.corr.r(:,c)));
        end
        ft_cluster_topoplot(settings.layout,real(alloutputs.ersp.corr.r(:,c)),settings.datasetinfo.label,...
            alloutputs.ersp.corr.p(:,c)',alloutputs.ersp.corr.stats{c}.mask);
    else
        if isempty(find(~isnan(alloutputs.ersp.corr.r(:,c))))
            alloutputs.ersp.corr.r(:,c) = zeros(size(alloutputs.ersp.corr.r(:,c)));
        end
        cluster_topoplot(real(alloutputs.ersp.corr.r(:,c)),settings.layout,...
            alloutputs.ersp.corr.p(:,c)',alloutputs.ersp.corr.stats{c}.mask);
    end
    colormap(lkcmap2)
    cbars2(c) = colorbar('EastOutside');
    FixAxes(gca,14)
    ax(c) = p(c,2).axis;
    ax2(c) =p(c,1).axis;
    end
    %cbar(c).Position = [ax2(c).Position(1)+ax2(c).Position(3)-cbar(c).Position(3) ax2(c).Position(2) cbar(c).Position(3) cbar(c).Position(4)];
end
Normalize_Clim(ax,1);


p.de.margin = [5 5 5 5];
p.marginleft = 18;
p.marginbottom = 25;
p.margintop = 8;
p.marginright = 10;
p.de.marginleft = 14;
% fix margins here

for c = 1:settings.nfreqs
    ax2(c) = p(c,1).axis;
    cbars2(c).Position = [ax2(c).Position(1)+0.8*ax2(c).Position(3) ax2(c).Position(2)+0.7*ax2(c).Position(4) 0.07*ax2(c).Position(3) 0.28*ax2(c).Position(4)];
end
Normalize_Ylim(ax2)

set(gcf,'Color','w')

savefig('Fig2c.fig')
export_fig('Fig2c.png','-m4')
save('Panel2c.mat','p')
close



%% Figure 3: ERP nonadditivity

figure

pos = get(gcf,'position');
set(gcf,'units','normalized','position',[0 0 1 1],'Color','w');

p = panel('no-manage-font');
p.pack('v',{50 50});
p(1).pack('h',{50 50})
p(2).pack('h',{25 50 25});

p(1,1).pack()
for cc = 1:4
    p(1,1).pack({[0.25*(cc-1) 0 0.25 0.15]})
end
p(1,1,1).select();
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.real(:,:,1,:),1))-...
    squeeze(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),1)),'b',0.15,1,'sem')
stdshade(t,squeeze(nanmean(allmeas{1}.nadderp.real(:,:,2,:),1))-...
    squeeze(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),1)),'r',0.15,1,'sem')
%Plot_sigmask(gca,alloutputs.erp.pt.stats{1}.mask,'cmapline','LineWidth',5)
H = sigstar({get(p(1,1,1).axis,'XLim')},2*min(min(alloutputs.erp.pt.stats{1}.prob)),0,18)
delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])

%FillBetween(t,nanmean(nanmean(allmeas{1}.nadderp.real(:,:,1,:),4),1)-...
%    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,1,:),4),1),nanmean(nanmean(allmeas{2}.nadderp.real(:,:,2,:),4),1)-...
%    nanmean(nanmean(allmeas{1}.nadderp.pseudo(:,:,2,:),4),1));
legend({'Corrected prestim low','Corrected prestim high'},'EdgeColor','none')
xlabel('Time (s)')
ylabel('Voltage (uV)')
title('Pseudotrial-based nonadditivity')
%ylabelunits(settings)
FixAxes(gca,16)
%axes('position',[0.75 0.135 0.15 0.2])
plotindx = linspace(0,max(settings.aucindex),5);
plotindx(1) = [];
plotindx = plotindx - settings.srate/10;
for cc = 1:4
    p(1,1,cc+1).select()
    plotdata = nanmean(squeeze(allmeas{1}.nadderp.diff(:,plotindx(cc),2,:)-allmeas{1}.nadderp.diff(:,plotindx(cc),1,:)),2);
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
            0.*alloutputs.erp.pt.stats{1}.mask(:,plotindx(cc)),0.*alloutputs.erp.pt.stats{1}.mask(:,plotindx(cc)));
    else
        cluster_topoplot(plotdata,settings.layout,...
            1-(0.*alloutputs.erp.pt.stats{1}.mask(:,plotindx(cc))),0.*alloutputs.erp.pt.stats{1}.mask(:,plotindx(cc)));
    end
    colormap(lkcmap2)
    if cc == 4
        colorbar
    end
    ax(cc) = p(1,1,cc+1).axis;
    title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
    Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
end
Normalize_Clim(ax,1)


p(1,2).pack();
for cc = 1:4
    p(1,2).pack({[0.25*(cc-1) 0 0.25 0.15]})
end
p(1,2,1).select();
t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
hold on
stdshade(t,squeeze(nanmean(allmeas{1}.ttv.real,1)),'k',0.15,1,'sem')
%Plot_sigmask(gca,alloutputs.erp.ttv.stats{1}.mask,'cmapline','LineWidth',5)
H = sigstar({get(p(1,2,1).axis,'XLim')},2*min(min(alloutputs.erp.ttv.stats{1}.prob)),0,18)
delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
xlabel('Time (s)')

ylabelunits(settings)
title('TTV-based nonadditivity')
%ylabelunits(settings)
FixAxes(gca,16)
%axes('position',[0.75 0.135 0.15 0.2])
plotindx = linspace(0,max(settings.aucindex),5);
plotindx(1) = [];
plotindx = plotindx - settings.srate/10;
for cc = 1:4
    p(1,2,cc+1).select()
    plotdata = nanmean(allmeas{1}.ttv.real(:,plotindx(cc),:),3);
    if strcmpi(settings.datatype,'MEG')
        ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
            0.*alloutputs.erp.ttv.stats{1}.mask(:,plotindx(cc)),0.*alloutputs.erp.ttv.stats{1}.mask(:,plotindx(cc)));
    else
        cluster_topoplot(plotdata,settings.layout,...
            1-(0.*alloutputs.erp.pt.stats{1}.mask(:,plotindx(cc))),0.*alloutputs.erp.pt.stats{1}.mask(:,plotindx(cc)));
    end
    colormap(lkcmap2)
    if cc == 4
        colorbar
    end
    ax(cc) = p(1,2,cc+1).axis;
    title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
    Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
end
Normalize_Clim(ax,1)

p(2,2).pack();
p(2,2).pack({[0.7 0 0.3 0.3]});
p(2,2,1).select();
nicecorrplot(nanmean(allmeas{1}.naerpindex,1),nanmean(allmeas{1}.ttvindex,1),{'Pseudotrial-based ERP nonadditivity','TTV-based ERP nonadditivity'});
FixAxes(gca,16)
title('Correlation of pseudotrial and TTV methods')
p(2,2,2).select()
plotdata = alloutputs.erp.corr.r(:,1);
if isempty(find(~isnan(plotdata)))
    plotdata = zeros(size(plotdata));
end

if strcmpi(settings.datatype,'MEG')
    ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
        alloutputs.erp.corr.p(:,1),alloutputs.erp.corr.stats{1}.mask);
else
    cluster_topoplot(plotdata,settings.layout,...
        1-(0.*alloutputs.erp.corr.p(:,1)),alloutputs.erp.corr.stats{1}.mask);
end
Normalize_Clim(gca,1)
colorbar('WestOutside')


p.margin = [20 20 5 8];
p.de.margin = [5 5 5 5];
p(1).marginright = 25;
p(1,2).marginleft = 28;
p(1).marginbottom = 32;

colormap(lkcmap2)
%Plot_sigmask(gca,alloutputs.erp.ttv.stats{1}.mask,'cmapline','LineWidth',5)

set(gcf,'Color','w')

%set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
savefig('Fig3.fig')
export_fig('Fig3.png','-m5')
save('Panel3.mat','p')
close

%% Figure 5b: IRASA osci and frac power time course on same plot


figure

pris = prism;

p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*3.5 pos(4)],'Color','w');

p.pack('h',repmat({1/settings_osci.nfreqs},settings_osci.nfreqs,1)')
for c = 1:settings_osci.nfreqs
    p(c).pack();
    for cc = 1:4
        p(c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(c,1).select();
    plotband = c;
    
    t = linspace(0,length(settings_osci.real.poststim)*(1/settings_osci.srate),length(settings_osci.real.poststim));
    %t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
    %prestimdata = 100*(allmeas{plotband}.raw.ttversp(:,prestim_real,:)-nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2))./...
    %    nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
    %prestimdata = nanmean(nanmean(prestimdata,3),1);
    poststimdata = (allmeas_osci{plotband}.ersp.real);
    poststimdata = nanmean(nanmean(poststimdata,3),1);
    %plotdata = [prestimdata poststimdata];
    plotdata = poststimdata;
    %FillBetween(t((length(settings.real.prestim)+1):end),poststimdata,...
    %    zeros(1,length(poststimdata)));
    hold on
    stdshade(t,squeeze(nanmean(allmeas_osci{plotband}.ersp.real,1)),pris(2,:),0.15,1,'sem')
    %plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
    hold on
    stdshade(t,squeeze(nanmean(allmeas_frac{plotband}.ersp.real,1)),pris(4,:),0.15,1,'sem')
    
    title(fbands{c})
    xlabel('Time (s)')
    if c ==1 
    ylabel('% change of ERSP')
    end
    %ylim = get(gca,'YLim');
    %line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
    %set(gca,'YLim',ylim)
    FixAxes(gca,14)
    set(gca,'XLim',[0 max(t)])
    %set(gca,'FontSize',16)
    
end
clear ax

for c = 1:settings_osci.nfreqs
   ax1(c) = p(c,1).axis; 
end
Normalize_Ylim(ax1)

for c = 1:settings_osci.nfreqs
    p(c,1).select();
    %Plot_sigmask(p(c,1).axis,compstats.erspdiff{c}.prob < 0.05,'cmapline','LineWidth',5)
    H = sigstar({get(p(c,1).axis,'XLim')},2*min(min(compstats.erspdiff{c}.prob)),0,18)
    delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)]);
end

p(settings_osci.nfreqs,1).select()
legend({'Oscillatory power','Fractal power'})

p.margintop = 8;
p.marginleft = 18;
p.de.marginleft = 12;

savefig('Fig5b.fig')
export_fig('Fig5b.png','-m5')
save('Panel5b.mat','p')

% %% Figure ??: IRASA osci and frac power time course
% % Need allmeas_osci and allmeas_frac, alloutputs_osci and alloutputs_frac
% 
% figure
% 
% pris = prism;
% 
% p = panel('no-manage-font');
% 
% pos = get(gcf,'position');
% set(gcf,'position',[pos(1:2) pos(3)*3.5 pos(4)*3.5],'Color','w');
% 
% p.pack('v',{50 50})
% p(1).pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
% for c = 1:settings.nfreqs
%     p(1,c).pack();
%     for cc = 1:4
%         p(1,c).pack({[0.25*(cc-1) 0 0.25 0.15]})
%     end
%     p(1,c,1).select();
%     plotband = c;
%     
%     t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
%     %t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
%     %prestimdata = 100*(allmeas{plotband}.raw.ttversp(:,prestim_real,:)-nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2))./...
%     %    nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
%     %prestimdata = nanmean(nanmean(prestimdata,3),1);
%     poststimdata = (allmeas_osci{plotband}.ersp.real);
%     poststimdata = nanmean(nanmean(poststimdata,3),1);
%     %plotdata = [prestimdata poststimdata];
%     plotdata = poststimdata;
%     %FillBetween(t((length(settings.real.prestim)+1):end),poststimdata,...
%     %    zeros(1,length(poststimdata)));
%     hold on
%     stdshade(t,squeeze(nanmean(allmeas_osci{plotband}.ersp.real,1)),'k',0.15,1,'sem')
%     %plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
%     
%     title(fbands{c})
%     xlabel('Time (s)')
%     if c == 1
%     ylabel('% change of Oscillatory ERSP')
%     end
%     %ylim = get(gca,'YLim');
%     %line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
%     %set(gca,'YLim',ylim)
%     FixAxes(gca,14)
%     set(gca,'XLim',[0 max(t)])
%     %set(gca,'FontSize',16)
%     
%     plotindx = linspace(0,max(settings.aucindex),5);
%     plotindx = round(plotindx);
%     plotindx(1) = [];
%     for cc = 1:4
%         p(1,c,cc+1).select()
%         plotdata = mean(allmeas_osci{c}.ersp.real(:,plotindx(cc),:),3);
%         if strcmpi(settings.datatype,'MEG')
%             ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
%                 0.*alloutputs_osci.ersp.ttv.stats{c}.mask(:,plotindx(cc)),0.*alloutputs_osci.ersp.ttv.stats{c}.mask(:,plotindx(cc)));
%         else
%             cluster_topoplot(plotdata,settings.layout,...
%                 1-(0.*alloutputs_osci.ersp.ttv.stats{c}.mask(:,plotindx(cc))),0.*alloutputs_osci.ersp.ttv.stats{c}.mask(:,plotindx(cc)));
%         end
%         colormap(lkcmap2)
%         if cc == 4
%             cbars1(c) = colorbar;
%         end
%         ax(cc) = p(1,c,cc+1).axis;
%         title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
%         Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
%     end
%     Normalize_Clim(ax,1);
% end
% clear ax
% 
% for c = 1:settings.nfreqs
%     ax1(c) = p(1,c,1).axis;
%     cbars1(c).Position = [ax1(c).Position(1)+ax1(c).Position(3) ax1(c).Position(2) cbars1(c).Position(3) 0.15*ax1(c).Position(4)];
% end
% %Normalize_Ylim(ax1)
% 
% p(2).pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
% 
% for c = 1:settings.nfreqs
%     p(2,c).pack();
%     for cc = 1:4
%         p(2,c).pack({[0.25*(cc-1) 0 0.25 0.15]})
%     end
%     p(2,c,1).select();
%     plotband = c;
%     
%     t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
%     %t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
%     %prestimdata = 100*(allmeas{plotband}.raw.ttversp(:,prestim_real,:)-nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2))./...
%     %    nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
%     %prestimdata = nanmean(nanmean(prestimdata,3),1);
%     poststimdata = (allmeas_frac{plotband}.ersp.real);
%     poststimdata = nanmean(nanmean(poststimdata,3),1);
%     %plotdata = [prestimdata poststimdata];
%     plotdata = poststimdata;
%     %FillBetween(t((length(settings.real.prestim)+1):end),poststimdata,...
%     %    zeros(1,length(poststimdata)));
%     hold on
%     stdshade(t,squeeze(nanmean(allmeas_frac{plotband}.ersp.real,1)),'k',0.15,1,'sem')
%     %plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
%     
%     title(fbands{c})
%     xlabel('Time (s)')
%     if c == 1
%     ylabel('% change of Fractal ERSP')
%     end
%     %ylim = get(gca,'YLim');
%     %line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
%     %set(gca,'YLim',ylim)
%     FixAxes(gca,14)
%     set(gca,'XLim',[0 max(t)])
%     %set(gca,'FontSize',16)
%     
%     plotindx = linspace(0,max(settings.aucindex),5);
%     plotindx = round(plotindx);
%     plotindx(1) = [];
%     for cc = 1:4
%         p(2,c,cc+1).select()
%         plotdata = mean(allmeas_frac{c}.ersp.real(:,plotindx(cc),:),3);
%         if strcmpi(settings.datatype,'MEG')
%             ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
%                 0.*alloutputs_frac.ersp.ttv.stats{c}.mask(:,plotindx(cc)),0.*alloutputs_frac.ersp.ttv.stats{c}.mask(:,plotindx(cc)));
%         else
%             cluster_topoplot(plotdata,settings.layout,...
%                 1-(0.*alloutputs_frac.ersp.ttv.stats{c}.mask(:,plotindx(cc))),0.*alloutputs_frac.ersp.ttv.stats{c}.mask(:,plotindx(cc)));
%         end
%         colormap(lkcmap2)
%         if cc == 4
%             cbars2(c) = colorbar;
%         end
%         ax(cc) = p(2,c,cc+1).axis;
%         title([num2str(plotindx(cc)*(1000/settings.srate)) ' ms'],'FontSize',10)
%         Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
%     end
%     Normalize_Clim(ax,1);
% end
% clear ax
% 
% for c = 1:settings.nfreqs
%     ax2(c) = p(2,c,1).axis;
%     cbars2(c).Position = [ax2(c).Position(1)+ax2(c).Position(3) ax2(c).Position(2) cbars2(c).Position(3) 0.15*ax2(c).Position(4)];
% 
%     %p(1,c,1).select()
%     %Plot_sigmask(p(1,c,1).axis,alloutputs.ersp.ttv.stats{c}.prob < 0.05,'cmapline','LineWidth',5)
% end
% Normalize_Ylim(cat(2,ax1,ax2))
% 
% for c = 1:settings.nfreqs
%     p(1,c,1).select()
%     %Plot_sigmask(p(1,c,1).axis,compstats.erspdiff{c}.mask,'cmapline','LineWidth',5)
%     H = sigstar({get(p(1,c,1).axis,'XLim')},2*min(compstats.erspdiff{c}.prob),0,18)
%     delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
%     
%     p(2,c,1).select()
%     %Plot_sigmask(p(2,c,1).axis,compstats.erspdiff{c}.mask,'cmapline','LineWidth',5)
%     H = sigstar({get(p(2,c,1).axis,'XLim')},2*min(compstats.erspdiff{c}.prob),0,18)
%     delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
% end
% 
% p.de.margin = [5 5 5 5];
% p.marginleft = 24;
% p.marginbottom = 25;
% p.margintop = 8;
% p.marginright = 10;
% p.de.marginleft = 15;
% p(1).marginbottom = 24;
% 
% AddFigureLabel(p(1,1,1).axis,'A')
% AddFigureLabel(p(2,1,1).axis,'B')
% 
% savefig('Fig6_of_tc.fig')
% export_fig('Fig6_of_tc.png','-m4')
% save('Panel6_of_tc.mat','p')

%% Figure 6a&b: Nonadditivity for IRASA oscillatory vs fractal

figure

p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'units','normalized','position',[0 0 1 1],'Color','w');

p.pack('v',{50 50})
p(1).pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
for c = 1:settings.nfreqs
        p(1,c).pack();
    for cc = 1:4
        p(1,c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(1,c,1).select();
    t = linspace(0,length(settings_osci.real.poststim)*(1/settings_osci.srate),length(settings_osci.real.poststim));
    hold on
    stdshade(t,squeeze(nanmean(allmeas_osci{c}.naddersp.diff(:,:,1,:),1)),'b',0.15,1,'sem');
    
    stdshade(t,squeeze(nanmean(allmeas_osci{c}.naddersp.diff(:,:,2,:),1)),'r',0.15,1,'sem');
    title(fbands{c})
    if c == settings.nfreqs
        legend({'Corrected prestim low','Corrected prestim high'})
    end
    xlabel('Time (s)')
    if c == 1
        ylabel('% change in Oscillatory ERSP')
    end
    FixAxes(gca,14)
    set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
        set(gca,'XLim',[0 max(t)])

    
        plotindx = linspace(0,max(settings_osci.aucindex),5);
        plotindx = movmean(plotindx,2);
    plotindx = round(plotindx);
    plotindx(1) = [];
    %tindx =
    for cc = 1:4
        p(1,c,cc+1).select()
        %axes(p(2,c,cc+1).axis)
        plotdata = nanmean(squeeze(allmeas_osci{c}.naddersp.diff(:,plotindx(cc),2,:))...
            - squeeze(allmeas_osci{c}.naddersp.diff(:,plotindx(cc),1,:)),2);
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                ones(size(alloutputs_osci.ersp.pt.stats{c}.mask(:,plotindx(cc)))),0.*alloutputs_osci.ersp.pt.stats{c}.mask(:,plotindx(cc)));
        else
            cluster_topoplot(plotdata,settings.layout,...
                ~(0.*alloutputs_osci.ersp.pt.stats{c}.mask(:,plotindx(cc))),0.*alloutputs_osci.ersp.pt.stats{c}.mask(:,plotindx(cc)));
        end
        colormap(lkcmap2)
        if cc == 4
            cbars1(c) = colorbar;
        end
        ax(cc) = p(1,c,cc+1).axis;
        title([num2str(plotindx(cc)*(1000/settings_osci.srate)) ' ms'],'FontSize',10)
        %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
    Normalize_Clim(ax,1);
end

for c = 1:settings.nfreqs
    ax1(c) = p(1,c,1).axis;
    cbars1(c).Position = [ax1(c).Position(1)+ax1(c).Position(3) ax1(c).Position(2) cbars1(c).Position(3) 0.15*ax1(c).Position(4)];
end
Normalize_Ylim(ax1)
%Set_Ylim(ax1,[-60 70])


for c = 1:settings.nfreqs
    p(1,c,1).select()
    %Plot_sigmask(p(1,c,1).axis,alloutputs_osci.ersp.pt.stats{c}.mask,'cmapline','LineWidth',5)
    H = sigstar({get(p(1,c,1).axis,'XLim')},2*min(min(alloutputs_osci.ersp.pt.stats{c}.prob)),0,18)
    delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
end

p(2).pack('h',repmat({1/settings.nfreqs},settings.nfreqs,1)')
for c = 1:settings.nfreqs
    p(2,c).pack();
    for cc = 1:4
        p(2,c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(2,c,1).select();
    t = linspace(0,length(settings_osci.real.poststim)*(1/settings_osci.srate),length(settings_osci.real.poststim));
    hold on
    stdshade(t,squeeze(nanmean(allmeas_frac{c}.naddersp.diff(:,:,1,:),1)),'b',0.15,1,'sem');
    
    stdshade(t,squeeze(nanmean(allmeas_frac{c}.naddersp.diff(:,:,2,:),1)),'r',0.15,1,'sem');
    title(fbands{c})
    if c == settings.nfreqs
        legend({'Corrected prestim low','Corrected prestim high'})
    end
    xlabel('Time (s)')
    if c == 1
        ylabel('% change in Fractal ERSP')
    end
    FixAxes(gca,14)
    set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
        set(gca,'XLim',[0 max(t)])

    
        plotindx = linspace(0,max(settings_osci.aucindex),5);
        plotindx = movmean(plotindx,2);
    plotindx = round(plotindx);
    plotindx(1) = [];
    %tindx =
    for cc = 1:4
        p(2,c,cc+1).select()
        %axes(p(2,c,cc+1).axis)
        plotdata = nanmean(squeeze(allmeas_frac{c}.naddersp.diff(:,plotindx(cc),2,:))...
            - squeeze(allmeas_frac{c}.naddersp.diff(:,plotindx(cc),1,:)),2);
        if strcmpi(settings.datatype,'MEG')
            ft_cluster_topoplot(settings.layout,plotdata,settings.datasetinfo.label,...
                ones(size(alloutputs_frac.ersp.pt.stats{c}.mask(:,plotindx(cc)))),0.*alloutputs_frac.ersp.pt.stats{c}.mask(:,plotindx(cc)));
        else
            cluster_topoplot(plotdata,settings.layout,...
                ~(0.*alloutputs_frac.ersp.pt.stats{c}.mask(:,plotindx(cc))),0.*alloutputs_frac.ersp.pt.stats{c}.mask(:,plotindx(cc)));
        end
        colormap(lkcmap2)
        if cc == 4
            cbars2(c) = colorbar;
        end
        ax(cc) = p(2,c,cc+1).axis;
        title([num2str(plotindx(cc)*(1000/settings_osci.srate)) ' ms'],'FontSize',10)
        %Set_Clim(ax(cc),[prctile(plotdata,20) prctile(plotdata,80)]);
    end
    Normalize_Clim(ax,1);
end

for c = 1:settings.nfreqs
    ax2(c) = p(2,c,1).axis;
    cbars2(c).Position = [ax2(c).Position(1)+ax2(c).Position(3) ax2(c).Position(2) cbars2(c).Position(3) 0.15*ax2(c).Position(4)];
end
Normalize_Ylim(ax2)
%Set_Ylim(ax2,[-25 40])

for c = 1:settings.nfreqs
    p(2,c,1).select()
    %Plot_sigmask(p(2,c,1).axis,alloutputs_frac.ersp.pt.stats{c}.mask,'cmapline','LineWidth',5)
    H = sigstar({get(p(2,c,1).axis,'XLim')},2*min(min(alloutputs_frac.ersp.pt.stats{c}.prob)),0,18);
    delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)]);
end

p.de.margin = [5 5 5 5];
p.marginleft = 18;
p.marginbottom = 18;
p.margintop = 8;
p.marginright = 12;
p.de.marginleft = 18;
p(1).marginbottom = 20;
p(1).de.marginleft = 18;
p(2).de.marginleft = 18;


AddFigureLabel(p(1,1,1).axis,'A');
AddFigureLabel(p(2,1,1).axis,'B');

savefig('Fig6ab.fig')
export_fig('Fig6ab.png','-m4')
save('Panel6ab.mat','p')

%% Figure 6c: Oscillatory-fractal nonadditivity time course

figure

pris = prism;

p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'units','normalized','position',[0 0.5 1 0.5],'Color','w');

p.pack('h',repmat({1/settings_osci.nfreqs},settings_osci.nfreqs,1)')
for c = 1:settings_osci.nfreqs
    p(c).pack();
    for cc = 1:4
        p(c).pack({[0.25*(cc-1) 0 0.25 0.15]})
    end
    p(c,1).select();
    plotband = c;
    
    t = linspace(0,length(settings_osci.real.poststim)*(1/settings_osci.srate),length(settings_osci.real.poststim));
    %t = linspace(-length(settings.real.prestim)*(1/settings.srate),length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim)+length(settings.real.prestim));
    %prestimdata = 100*(allmeas{plotband}.raw.ttversp(:,prestim_real,:)-nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2))./...
    %    nanmean(allmeas{plotband}.raw.ttversp(:,prestim_real,:),2); %assuming percent change units
    %prestimdata = nanmean(nanmean(prestimdata,3),1);
    poststimdata = squeeze(diff(allmeas_osci{plotband}.naddersp.diff,1,3));
    poststimdata = nanmean(nanmean(poststimdata,3),1);
    %plotdata = [prestimdata poststimdata];
    plotdata = poststimdata;
    %FillBetween(t((length(settings.real.prestim)+1):end),poststimdata,...
    %    zeros(1,length(poststimdata)));
    hold on
    stdshade(t,squeeze(nanmean(diff(allmeas_osci{plotband}.naddersp.diff,1,3),1)),pris(2,:),0.15,1,'sem')
    %plot(t,zeros(1,length(plotdata)),'k--','LineWidth',1.5)
    hold on
    stdshade(t,squeeze(nanmean(diff(allmeas_frac{plotband}.naddersp.diff,1,3),1)),pris(4,:),0.15,1,'sem')
    
    title(fbands{c})
    xlabel('Time (s)')
    if c ==1 
    ylabel('Nonadditive effect (%)')
    end
    %ylim = get(gca,'YLim');
    %line([0 0],ylim,'Color',[0.5 0.5 0.5],'LineWidth',2)
    %set(gca,'YLim',ylim)
    FixAxes(gca,14)
    set(gca,'XLim',[0 max(t)])
    %set(gca,'FontSize',16)
    
end
clear ax

for c = 1:settings_osci.nfreqs
   ax1(c) = p(c,1).axis; 
end
Normalize_Ylim(ax1)

for c = 1:settings_osci.nfreqs
    p(c,1).select();
    %Plot_sigmask(p(c,1).axis,compstats.erspdiff{c}.prob < 0.05,'cmapline','LineWidth',5)
    H = sigstar({get(p(c,1).axis,'XLim')},2*min(min(compstats.ptdiff{c}.prob)),0,18)
    delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)]);
end

p(settings_osci.nfreqs,1).select()
legend({'Oscillatory nonadditive effect','Fractal nonadditive effect'},'FontSize',10)

p.margintop = 8;
p.marginleft = 18;
p.de.marginleft = 12;

savefig('Fig6c.fig')
export_fig('Fig6c.png','-m4')
save('Panel6c.mat','p')


%% Figure xxx (supplement) - pseudotrial-based and ttv-based time course differences?

%Pseudotrial based first
p = panel('no-manage-font');

set(gcf,'position',[pos(1:2) pos(3)*3 pos(4)*3],'Color','w');

p.pack(settings.nfreqs,settings.nfreqs);

for q = 1:settings.nfreqs-1
    for qq = 1:settings.nfreqs-1
        if qq > q %below the diagonal = pseudotrial
            p(q,qq).select();
            t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
            hold on
            stdshade(t,squeeze(nanmean(allmeas{q}.naddersp.diff(:,:,2,:),1))-...
                squeeze(nanmean(allmeas{q}.naddersp.diff(:,:,1,:))),'b',0.15,1,'std');
            stdshade(t,squeeze(nanmean(allmeas{qq}.naddersp.diff(:,:,2,:),1))-...
                squeeze(nanmean(allmeas{qq}.naddersp.diff(:,:,1,:))),'r',0.15,1,'std');
            %Plot_sigmask(p(q,qq).axis,alloutputs.ersp.pt.tcoursestats{q,qq}.prob < 0.05,'cmapline','LineWidth',5)
            H = sigstar({get(p(q,qq).axis,'XLim')},2*min(alloutputs.ersp.pt.tcoursestats{c}.prob),0,18)
            delete(H(1)); pos = get(H(2),'position'); yl = get(gca,'YLim'); set(H(2),'position',[pos(1) yl(2)-0.05*diff(yl) pos(3)])
            
            if q == 1
                title(fbands{qq})
            end
            %legend({'Corrected prestim low','Corrected prestim high'})
            xlabel('Time (s)')
            ylabelunits(settings)
            FixAxes(gca)
            set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
            
            
        elseif qq < q %above the diagonal, plot TTV
            p(q,qq).select()
            t = linspace(0,length(settings.real.poststim)*(1/settings.srate),length(settings.real.poststim));
            hold on
            stdshade(t,squeeze(nanmean(allmeas{q}.ttversp.real,1)),'b',0.15,1,'std');
            stdshade(t,squeeze(nanmean(allmeas{qq}.ttversp.real,1)),'r',0.15,1,'std');
            
            Plot_sigmask(p(q,qq).axis,2*min(alloutputs.ersp.pt.tcoursestats{qq,q}.prob),'cmapline','LineWidth',5)
            
            if q == 1
                title(fbands{qq})
            end
            %legend({'Corrected prestim low','Corrected prestim high'})
            xlabel('Time (s)')
            ylabelunits(settings)
            FixAxes(gca)
            set(gca,'FontSize',11,'TitleFontSizeMultiplier',1.1)
            
        elseif qq == q % on the diagonal, just plot a line
            p(q,qq).select()
            % plot(linspace(0,1,100),linspace(1,0,100),'k','LineWidth',5)
            axis(p(q,qq).axis,'off')
            if q == 1
                title(fbands{q},'FontSize',11)
            end
            
        end
    end
end

% fix margins, add text boxes etc

savefig('FigS1.fig')
export_fig('FigS1.png','-m4')
save('PanelS1.mat','p')
close
end

function ylabelunits(settings)
switch settings.units
    case 'prcchange'
        ylabel('% change from prestim')
    case 'log'
        ylabel('10*log10 unit change')
    case 'raw'
        ylabel('Change from prestim (respective units)')
    case 'zscore'
        ylabel('Normalized change from prestim')
end
end

