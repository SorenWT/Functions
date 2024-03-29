function [pvals] = Plot_figure_oscifrac_final(spec,fbands,opts)
% Plots a standardized figure for the OsciFrac paper
%
% Input arguments:
%    spec: a structure array. Each field should be titled according to the
%       condition, and should contain an array of IRASA spectra for each
%       subject
%    fbands: an nx2 array containing the definitions of the frequency bands
%       of interest
%    opts: a structure with the following fields:
%       paired: 'yes' or 'no' for paired statistics. Currently not used
%          (default = 'no')
%       statfields: fields to do statistics on (default = all the fields in
%          the specs structure). The order of this field also specifies the
%          reference category - the last field is taken as the reference
%       fbandnames: names of the frequency bands you're using (default =
%          {'Delta','Theta','Alpha','Beta','Gamma'})
%       frange_ple: a 1x2 array containing the range over which the PLE
%          should be calculated (estimate this graphically based on when
%          the fractal power "drops off"; default = [lowest frequency in
%          fbands highest frequency in fbands])
%       stats_mixd: can input the stats from a previous call to this
%          function in order to do the plot faster
%       stats_oscifrac: can input the stats from a previous call to this
%          function in order to do the plot faster
%       cmpare: can input the stats from a previous call to this
%          function in order to do the plot faster

%% Setting up options
fields = fieldnames(spec);

if nargin < 3
    opts = struct;
end

opts = setdefault(opts,'paired','no');
opts = setdefault(opts,'fbandnames',{'Delta','Theta','Alpha','Beta','Gamma'});
opts = setdefault(opts,'statfields',fields);
opts = setdefault(opts,'frange_ple',[fbands(1,1) fbands(end,end)]);
opts = setdefault(opts,'powermode','abs');

frange = intersect(find(spec.(fields{1})(1).freq(:,1) > opts.frange_ple(1)),find(spec.(fields{1})(1).freq(:,1) < opts.frange_ple(1,2)));

%% Calculating power
for i = 1:length(fields)
    for c = 1:length(spec.(fields{i}))
        for cc = 1:size(fbands,1)
            if strcmpi(opts.powermode,'rel')
                relmixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),{'mixd',[fbands(1,1),fbands(end,end)]},'trapz');
            end
            absmixd.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'mixd',fbands(cc,:),0,'mean');
        end
        
        for cc = 1:size(fbands,1)
            if strcmpi(opts.powermode,'rel')
                relosci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),{'osci',[fbands(1,1),fbands(end,end)]},'trapz');
            end
            absosci.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',fbands(cc,:),0,'mean');
        end
        
        for cc = 1:size(fbands,1)
            absfrac.(fields{i})(c,:,cc) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'frac',fbands(cc,:),0,'mean');
        end
        
        tmp = amri_sig_plawfit(spec.(fields{i})(c),opts.frange_ple);
        ple.(fields{i})(c,:) = tmp.Beta;
        if strcmpi(opts.powermode,'rel')
            bb.(fields{i})(c,:) = 1./OsciFrac_EEG_wrapper(spec.(fields{i})(c),[fbands(1,1) fbands(end,end)]);
        else
            bb.(fields{i})(c,:) = IRASAPower_EEG_wrapper(spec.(fields{i})(c),'osci',[fbands(1,1) fbands(end,end)],0,'mean');
        end
    end
end

if strcmpi(opts.powermode,'rel')
    mixd = relmixd;
    osci = relosci;
else
    mixd = absmixd;
    osci = absosci;
end

% for c = 1:length(fields)
%     spec.(fields{c}) = mergestructs(spec.(fields{c}));
%     spec.(fields{c}) = structfun(@(data)nanmean(data,3),spec.(fields{c}),'UniformOutput',false);
% end

%% Logistic regression


data = [];
for i = 1:length(opts.statfields)
    tmp = [squeeze(mean(mixd.(opts.statfields{i}),2)) squeeze(mean(osci.(opts.statfields{i}),2)) ...
        mean(ple.(opts.statfields{i}),2) mean(bb.(opts.statfields{i}),2) ...
        ones(size(mean(ple.(opts.statfields{i}),2)))*i];
    data = cat(1,data,tmp);
end

if strcmpi(opts.paired,'yes') %adding subjects for random effects
    data = cat(2,data,repmat([1:size(tmp,1)]',length(fields),1));
end

if strcmpi(opts.paired,'yes')
    varnames = [cellcat('Mixd_',opts.fbandnames,'',0) ...
        cellcat('Osci_',opts.fbandnames,'',0) {'PLE' 'Frac_BB' 'Condition' 'Subject'}];
else
    varnames = [cellcat('Mixd_',opts.fbandnames,'',0) ...
        cellcat('Osci_',opts.fbandnames,'',0) {'PLE' 'Frac_BB' 'Condition'}];
end

data = array2table(data,'VariableNames',varnames);

mixnames = cellcat('+',cellcat('Mixd_',opts.fbandnames,'',0),'',1);
mixnames = [mixnames{:}];
mixnames(end) = [];

oscinames = cellcat('+',cellcat('Osci_',opts.fbandnames,'',0),'',1);
oscinames = [oscinames{:}];

frmla1 = ['Condition ~ ' mixnames];
frmla2 = ['Condition ~ ' oscinames 'PLE+Frac_BB'];

if strcmpi(opts.paired, 'yes')
    frmla1 = [frmla1 '+(1|Subject)-1'];
    frmla2 = [frmla2 '+(1|Subject)-1'];
end

if strcmpi(opts.paired,'yes')
    data{:,1:end-2} = zscore(data{:,1:end-2},[],1);
else
    data{:,1:end-1} = zscore(data{:,1:end-1},[],1);
end

if strcmpi(opts.paired,'yes')
    for c = 1:size(fbands,1)
        %         t = rm_anova2(cat(1,data.(['Mixd_' opts.fbandnames{c}]),data.(['Osci_' opts.fbandnames{c}])),...
        %            cat(1,data.Subject,data.Subject),cat(1,data.Condition,data.Condition),...
        %            Make_designVect([height(data) height(data)])',{'Condition','Oscifrac'});
        %                 pvals.p_dif_irasa(c) = t{4,6};
        
        t = Permtest_rmanova2(cat(1,data.(['Mixd_' opts.fbandnames{c}]),data.(['Osci_' opts.fbandnames{c}])),...
            cat(1,data.Subject,data.Subject),cat(1,data.Condition,data.Condition),...
            Make_designVect([height(data) height(data)])',{'Condition','Oscifrac'},1000);
        pvals.p_dif_irasa(c) = t(end);
        
        
        if length(fields) > 2
            tbl2 = [];
            for cc = 1:length(fields)
                tbl2 = [tbl2 mean(mixd.(fields{cc})(:,:,c),2)];
            end
            pvals.p_mixd(c) = friedman(tbl2,1,'off');
            
            pvals.p_mixd_mcompare(:,:,c) = mcompare_bfholm(tbl2,@signrank);
            
            tbl2 = [];
            for cc = 1:length(fields)
                tbl2 = [tbl2 mean(osci.(fields{cc})(:,:,c),2)];
            end
            pvals.p_osci(c) = friedman(tbl2,1,'off');
            
            pvals.p_osci_mcompare(:,:,c) = mcompare_bfholm(tbl2,@signrank);
            
            
        else
            pvals.p_mixd(c) = signrank(mean(mixd.(fields{1})(:,:,c),2),mean(mixd.(fields{2})(:,:,c),2));
            pvals.p_osci(c) = signrank(mean(osci.(fields{1})(:,:,c),2),mean(osci.(fields{2})(:,:,c),2));
        end
    end
    
    if length(fields) > 2
        tbl2 = [];
        for cc = 1:length(fields)
            tbl2 = [tbl2 mean(ple.(fields{cc}),2)];
        end
        pvals.p_ple = friedman(tbl2,1,'off');
        
        pvals.p_ple_mcompare = mcompare_bfholm(tbl2,@signrank);
        
        tbl2 = [];
        for cc = 1:length(fields)
            tbl2 = [tbl2 mean(bb.(fields{cc}),2)];
        end
        pvals.p_bb = friedman(tbl2,1,'off');
        
        pvals.p_bb_mcompare = mcompare_bfholm(tbl2,@signrank);
        
    else
        pvals.p_ple = signrank(mean(ple.(fields{1}),2),mean(ple.(fields{2}),2));
        pvals.p_bb = signrank(mean(bb.(fields{1}),2),mean(bb.(fields{2}),2));
    end
else
    for c = 1:size(fbands,1)
        %[~,t] = anovan(cat(1,data.(['Mixd_' opts.fbandnames{c}]),data.(['Osci_' opts.fbandnames{c}])),...
        %    {Make_designVect([height(data) height(data)])',cat(1,data.Condition,data.Condition)},'model','interaction');
        %pvals.p_dif_irasa(c) = t{4,7};
        t = Permtest_anovan(cat(1,data.(['Mixd_' opts.fbandnames{c}]),data.(['Osci_' opts.fbandnames{c}])),...
            {Make_designVect([height(data) height(data)])',cat(1,data.Condition,data.Condition)},1000);
        pvals.p_dif_irasa(c) = t(end);
        pvals.p_mixd(c) = ranksum(mean(mixd.(fields{1})(:,:,c),2),mean(mixd.(fields{2})(:,:,c),2));
        pvals.p_osci(c) = ranksum(mean(osci.(fields{1})(:,:,c),2),mean(osci.(fields{2})(:,:,c),2));
    end
    
    pvals.p_ple = ranksum(mean(ple.(fields{1}),2),mean(ple.(fields{2}),2));
    pvals.p_bb = ranksum(mean(bb.(fields{1}),2),mean(bb.(fields{2}),2));
end

% %% Do multiple comparison correction across the whole dataset
%
% pfields = fieldnames(pvals);
% for c = 1:length(pfields)
%
% end

%% Plotting the figure

for c = 1:length(fields)
    spec.(fields{c}) = mergestructs(spec.(fields{c}));
    spec.(fields{c}) = structfun(@(data)nanmean(data,3),spec.(fields{c}),'UniformOutput',false);
end

p = panel('no-manage-font');

pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*2 pos(4)*3],'color','w');

p.de.margin = [5 5 5 5];

%fix margins
p.marginleft = 18;
p.margintop = 8;

p.pack('h',{50 50})

p(1).pack('v',{1/2 1/2})
%p(1).pack('v',{40 30 30})
p(1,1).marginbottom = 24;
p(1,2).pack('h',{1/2 1/2})

cmap = lines;


p(1,1).select()
hold on
for c = 1:length(fields)
    if strcmpi(opts.powermode,'abs')
        plot(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).mixd(frange,:),2),'LineWidth',1.5)
    else
        plot(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).mixd(frange,:),2)./...
            nanmean(nanmean(spec.(fields{c}).mixd(frange,:))),'LineWidth',1.5)
    end
end
legend(fields)
FixAxes(gca,12)
title('Mixed power spectrum')
xlabel('Log Frequency (Hz)')
ylabel('Log Power')
set(gca,'XScale','log','YScale','log','XLim',opts.frange_ple)

p(1,2,1).select()
hold on
for c = 1:length(fields)
    if strcmpi(opts.powermode,'abs')
         loglog(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).osci(frange,:),2),'LineWidth',1.5)
    else
         loglog(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).osci(frange,:),2)./...
            nanmean(nanmean(spec.(fields{c}).mixd(frange,:))),'LineWidth',1.5) 
    end
end
FixAxes(gca,12)
title('Oscillatory power spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'XLim',opts.frange_ple)

p(1,2,2).select()
hold on
for c = 1:length(fields)
    if strcmpi(opts.powermode,'abs')
    loglog(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).frac(frange,:),2),'LineWidth',1.5)
    else
               loglog(spec.(fields{c}).freq(frange,1),nanmean(spec.(fields{c}).frac(frange,:),2)./...
            nanmean(nanmean(spec.(fields{c}).mixd(frange,:))),'LineWidth',1.5) 
    end
end
FixAxes(gca,12)
title('Fractal power spectrum')
xlabel('Log Frequency (Hz)')
ylabel('Log Power')
set(gca,'XScale','log','YScale','log','XLim',opts.frange_ple)
%Normalize_YLim(cat(1,p(1,2,1).axis,p(1,2,2).axis))


%p(2,1).pack(2,6)
%p(2,2).pack(2,6)
%p(2,2).pack('h',{1/6 1/6 1/6 1/6 1/6 1/6})

plotstruct = struct;
for c = 1:length(fields)
    plotstruct.(fields{c}) = [];
end

fbandnames = opts.fbandnames;
%ple_fbands = {'2-8 Hz','8-13 Hz','13-50 Hz'};
%ple_fbands = {'2-50 Hz'};

% if ~exist('plotstats','var')
%     plotstats = cell(2,6);
% end

% if ~strcmpi(statfun,'ft_statfun_friedman') && ~strcmpi(statfun,'ft_statfun_kruskal')
%     numstat = 2;
% else
%     numstat = length(fields);
% end

for c = 1:size(fbands,1)
    for cc = 1:length(fields)
        plotstack(c,cc,1) = nanmean(nanmean(absfrac.(fields{cc})(:,:,c),2),1);
        plotstack(c,cc,2) = nanmean(nanmean(absosci.(fields{cc})(:,:,c),2),1);
        %plotstack(c,cc,2) = nanmean(nanmean(linfrac.(fields{cc})(:,:,c),2),1);
        %plotstack(c,cc,3) = nanmean(nanmean(resfrac.(fields{cc})(:,:,c),2),1);
    end
end

for c = 1:size(fbands,1)
    for cc = 1:length(fields)
        plotstack(c,cc,1) = nanmean(nanmean(absfrac.(fields{cc})(:,:,c),2),1);
        plotstack(c,cc,2) = nanmean(nanmean(absosci.(fields{cc})(:,:,c),2),1);
        mixsum(c,cc) = nanmean(nanmean(absmixd.(fields{cc})(:,:,c),2),1);
    end
end

l = lines;
l = l(1:length(fields),:);
for c = 1:length(fields)
    palel(c,:) = palecol(l(c,:));
end
for c = 1:size(fbands,1)
    col{c} = cat(3,l,palel);
    col{c} = permute(col{c},[3 1 2]);
end

lgnd = [];
for c = 1:length(fields)
    lgnd = [lgnd {[fields{c} ' Fractal']} {[fields{c} ' Oscillatory']}];
end

p(2).pack(floor(size(fbands,1)/2),ceil(size(fbands,1)/2))
linindx = cat(2,Make_designVect(repmat(ceil(size(fbands,1)/2)',1,floor(size(fbands,1)/2)))',...
    repmat((1:ceil(size(fbands,1)/2))',floor(size(fbands,1)/2),1));

for c = 1:size(fbands,1)
    p(2,linindx(c,1),linindx(c,2)).select()
    b = bar([1:length(fields)],squeeze(plotstack(c,:,:)),'stacked','HandleVisibility','off');
    for cc = 1:2
        b(cc).FaceColor = 'flat';
    end
    b(1).CData = l;
    b(2).CData = palel;
    hold on;
    for cc = 1:length(fields)
        if size(mixd.(fields{cc}),1) > 1
            tmp = squeeze(nanstd(nanmean(absmixd.(fields{cc}),2),[],1)./sqrt(size(absmixd.(fields{cc}),1)));
            er = errorbar(cc,mixsum(c,cc),tmp(c),'HandleVisibility','off');
            er.Color = [0 0 0];
            er.LineStyle = 'none';
        end
    end
    if c == ceil(size(fbands,1)/2)
        tmp = cat(1,l,palel);
        indx = [];
        for cc = 1:length(fields)
            indx = cat(2,indx,[cc:length(fields):size(tmp,1)]);
        end
        tmp = tmp(indx,:);
       manlegend(lgnd,tmp)
    end
    FixAxes(gca,14)
ylabel('Power')
set(gca,'XTickLabel',fields)
title(fbandnames{c})
end

% 
% p(2).select()
% [h,barpos] = plotBarStackGroups(plotstack,opts.fbandnames,col);
% hold on;
% for cc = 1:length(fields)
%     if size(mixd.(fields{cc}),1) > 1
%         er = errorbar(barpos(cc,:),mixsum(:,cc),squeeze(nanstd(nanmean(absmixd.(fields{cc}),2),[],1)./sqrt(size(absmixd.(fields{cc}),1))));
%         er.Color = [0 0 0];
%         er.LineStyle = 'none';
%     end
% end
% legend(l)
% FixAxes(gca,14)
% ylabel('Power')
% set(gca,'YScale','log')
if isfield(opts,'filename')
    savefig([opts.filename '_a.fig'])
    export_fig([opts.filename '_a.png'],'-m4')
end

figure

p = panel('no-manage-font')
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) pos(3)*2 pos(4)*3],'color','w');


p.pack('v',{75 25})
p(1).pack('v',{50 50})
p(2).pack('h',{50 50})
p(2).de.margin = [10 10 5 5];

%p(2,2,1).pack('v',{50 50})

l = lines;

p(1,1).select()


%scalefact = 0.8/length(fields);
h = cell(1,length(fields));

if length(fields) == 2
    s = 2;
else
    s = 3;
end
for c = 1:length(fields)
    %     tmpstruct = mixd;
    %     for cc = 1:length(fields)
    %         tmpstruct.(fields{cc}) = tmpstruct.(fields{cc})(:,:,c);
    %     end
    %     violinplot(tmpstruct)
    plotdata = mean(mixd.(fields{c}),2);
    plotdata(find(plotdata <= 0)) = min(min(plotdata));
    h{c} = notBoxPlot(squeeze(log10(plotdata)),(1:s:s*size(fbands,1))+((c-1)*0.4));
    for cc = 1:size(fbands,1)
        set(h{c}(cc).sdPtch,'FaceColor',palel(c,:),'EdgeColor',palel(c,:).*0.75,'HandleVisibility','off')
        set(h{c}(cc).semPtch,'FaceColor',l(c,:),'EdgeColor',l(c,:).*0.75,'HandleVisibility','off')
        set(h{c}(cc).mu,'Color',[0 0 0],'HandleVisibility','off')
        set(h{c}(cc).data,'HandleVisibility','off')
    end
    set(h{c}(1).semPtch,'HandleVisibility','on')
end
xl = get(gca,'XLim');
set(gca,'XLim',[0.5 max(xl)])

if length(fields) == 2
    %hard coded for now
    barlocs = {[1 1.4],[3 3.4],[5 5.4],[7 7.4]};
    sigstar(barlocs(find(pvals.p_mixd < 0.05)),pvals.p_mixd(find(pvals.p_mixd < 0.05)),0,18)
else
    for c = 1:size(fbands,1)
        barpos = (1+s*(c-1)):0.4:(1+s*(c-1)+0.4*(length(fields)-1));
        barmat = cat(3,repmat(barpos,5,1),repmat(barpos',1,5));
        barcell = mat2cell(barmat,ones(5,1),ones(5,1),2);
        barcell = cellfun(@squeeze,barcell,'UniformOutput',false);
        tmpp = pvals.p_mixd_mcompare(:,:,c);
        tmpp = tmpp(find(belowDiag(ones(5))));
        barcell = barcell(find(belowDiag(ones(5))));
        barcell = barcell(find(tmpp < 0.05));
        tmpp = tmpp(find(tmpp < 0.05));
        sigstar(barcell,tmpp,0,18);
    end
end

ax = gca;
ax.XAxis.Visible = 'off';

ylabel('Log mixed power')
legend(fields)
FixAxes(gca,14)
ax(1) = gca;

p(1,2).select()
for c = 1:length(fields)
    %     tmpstruct = mixd;
    %     for cc = 1:length(fields)
    %         tmpstruct.(fields{cc}) = tmpstruct.(fields{cc})(:,:,c);
    %     end
    %     violinplot(tmpstruct)
    plotdata = mean(osci.(fields{c}),2);
    plotdata(find(plotdata <= 0)) = min(min(plotdata(find(plotdata > 0)))); % fix this later
    h{c} = notBoxPlot(squeeze(log10(plotdata)),(1:s:s*size(fbands,1))+((c-1)*0.4));
    for cc = 1:size(fbands,1)
        set(h{c}(cc).sdPtch,'FaceColor',palel(c,:),'EdgeColor',palel(c,:).*0.75,'HandleVisibility','off')
        set(h{c}(cc).semPtch,'FaceColor',l(c,:),'EdgeColor',l(c,:).*0.75,'HandleVisibility','off')
        set(h{c}(cc).mu,'Color',[0 0 0],'HandleVisibility','off')
        set(h{c}(cc).data,'HandleVisibility','off')
    end
    set(h{c}(1).semPtch,'HandleVisibility','on')
end
xl = get(gca,'XLim')
set(gca,'XLim',[0.5 max(xl)])

if length(fields) == 2
    %hard coded for now
    barlocs = {[1 1.4],[3 3.4],[5 5.4],[7 7.4]};
    sigstar(barlocs(find(pvals.p_osci < 0.05)),pvals.p_osci(find(pvals.p_osci < 0.05)),0,18)
else
    for c = 1:size(fbands,1)
        barpos = (1+s*(c-1)):0.4:(1+s*(c-1)+0.4*(length(fields)-1));
        barmat = cat(3,repmat(barpos,5,1),repmat(barpos',1,5))
        barcell = mat2cell(barmat,ones(5,1),ones(5,1),2);
        barcell = cellfun(@squeeze,barcell,'UniformOutput',false);
        tmpp = pvals.p_osci_mcompare(:,:,c);
        tmpp = tmpp(find(belowDiag(ones(5))));
        barcell = barcell(find(belowDiag(ones(5))));
        barcell = barcell(find(tmpp < 0.05));
        tmpp = tmpp(find(tmpp < 0.05));
        sigstar(barcell,tmpp,0,18);
    end
end

set(gca,'XTick',(1:s:s*size(fbands,1))+(0.4*(length(fields)-1)/2),'XTickLabel',opts.fbandnames)

ylabel('Log oscillatory power')
FixAxes(gca,14)
ax(2) = gca;

Normalize_Ylim(ax)

p(2,1).select()
if strcmpi(opts.paired,'yes')
    plotstruct = table;
end
for cc = 1:length(fields)
    plotstruct.(fields{cc}) = squeeze(nanmean(bb.(fields{cc}),2));
end
clear h
if strcmpi(opts.paired,'yes')
    h=notBoxPlot(plotstruct{:,:})
else
    h(1) = notBoxPlot(plotstruct.(fields{1}),1);
    h(2) = notBoxPlot(plotstruct.(fields{2}),2);
end
for c = 1:length(fields)
    set(h(c).sdPtch,'FaceColor',palel(c,:),'EdgeColor',palel(c,:).*0.75,'HandleVisibility','off')
    set(h(c).semPtch,'FaceColor',l(c,:),'EdgeColor',l(c,:).*0.75)
    set(h(c).mu,'Color',[0 0 0],'HandleVisibility','off')
end
if length(fields) == 2
    sigstar({[1 2]},pvals.p_bb,0,18)
else
    barpos = 1:length(fields);
    barmat = cat(3,repmat(barpos,5,1),repmat(barpos',1,5))
    barcell = mat2cell(barmat,ones(5,1),ones(5,1),2);
    barcell = cellfun(@squeeze,barcell,'UniformOutput',false);
    tmpp = pvals.p_bb_mcompare;
    tmpp = tmpp(find(belowDiag(ones(5))));
    barcell = barcell(find(belowDiag(ones(5))));
    barcell = barcell(find(tmpp < 0.05));
    tmpp = tmpp(find(tmpp < 0.05));
    sigstar(barcell,tmpp,0,18);
end
title('Fractal broadband power')
set(gca,'XTickLabel',fields)
ylabel('Fractal broadband power')
FixAxes(gca,14)


p(2,2).select()
if strcmpi(opts.paired,'yes')
    plotstruct = table;
end
for cc = 1:length(fields)
    plotstruct.(fields{cc}) = squeeze(nanmean(ple.(fields{cc}),2));
end
clear h
if strcmpi(opts.paired,'yes')
    h=notBoxPlot(plotstruct{:,:})
else
    h(1) = notBoxPlot(plotstruct.(fields{1}),1);
    h(2) = notBoxPlot(plotstruct.(fields{2}),2);
end
for c = 1:length(fields)
    set(h(c).sdPtch,'FaceColor',palel(c,:),'EdgeColor',palel(c,:).*0.75,'HandleVisibility','off')
    set(h(c).semPtch,'FaceColor',l(c,:),'EdgeColor',l(c,:).*0.75)
    set(h(c).mu,'Color',[0 0 0],'HandleVisibility','off')
end
if length(fields) == 2
    sigstar({[1 2]},pvals.p_ple,0,18)
else
    barpos = 1:length(fields);
    barmat = cat(3,repmat(barpos,5,1),repmat(barpos',1,5))
    barcell = mat2cell(barmat,ones(5,1),ones(5,1),2);
    barcell = cellfun(@squeeze,barcell,'UniformOutput',false);
    tmpp = pvals.p_ple_mcompare;
    tmpp = tmpp(find(belowDiag(ones(5))));
    barcell = barcell(find(belowDiag(ones(5))));
    barcell = barcell(find(tmpp < 0.05));
    tmpp = tmpp(find(tmpp < 0.05));
    sigstar(barcell,tmpp,0,18);
end
%violinplot(plotstruct)
title('Fractal PLE')
set(gca,'XTickLabel',fields)
ylabel('PLE')
FixAxes(gca,14)

set(gcf,'Color','w')

p.marginleft = 20;
p(1,1).marginbottom = 5;
p(1).marginbottom = 18;
p(2,1).marginright = 12;

if isfield(opts,'filename')
    savefig([opts.filename '_b.fig'])
    export_fig([opts.filename '_b.png'],'-m4')
end

% pvals = struct;
% pvals.pvals.p_irasa_dif = pvals.p_dif_irasa;
% pvals.pvals.p_mixd = pvals.p_mixd;
% pvals.pvals.p_osci = pvals.p_osci;
% pvals.pvals.p_ple = pvals.p_ple;
% pvals.pvals.p_bb = pvals.p_bb;
%
% if length(fields) > 2
%    pvals.pvals.p_mixd_mcompare = pvals.p_mixd_mcompare;
%    pvals.pvals.p_osci_mcompare = pvals.p_osci_mcompare;
%    pvals.pvals.p_ple_mcompare = pvals.p_ple_mcompare;
%    pvals.pvals.p_bb_mcompare = pvals.p_bb_mcompare;
% end


