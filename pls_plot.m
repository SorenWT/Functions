function p = pls_plot(plsmdl,whichcomps,datasetlabels,xlabels,ylabels,varargin)
% this function makes a standard plot for a pls model defined by plsregress_perm
%
% Inputs:
%     plsmdl: a model from the plsregress_perm function
%     ncomps: number of components to plot
%     plotlabels: labels for the X and Y datasets
%     xlabels: labels for the x variables
%     ylabels: labels for the y variables
% Optional inputs:
%     panel: plot into a particular panel
%     panelindx: plot into a particular part of a panel
%     coeffplots: the types of coefficient plots to be used. Default is
%     {{'Xloads','bar'},{'Yloads','bar'}}. Can add more types if desired



argsin = varargin;

if CheckInput(argsin,'panel')
    p = EasyParse(argsin,'panel');
else
    figure
    set(gcf,'units','normalized','position',[0 0.4-0.4*double(length(whichcomps)>1) 1 0.6+0.4*double(length(whichcomps)>1)]);
    p = panel('no-manage-font');
end

if CheckInput(argsin,'panelindx')
    pindx = EasyParse(argsin,'panelindx');
else
    pindx = {};
end

if CheckInput(argsin,'coeffplots')
    coeffplots = EasyParse(argsin,'coeffplots');
else
    coeffplots = {{'Xloads','bar'},{'Yloads','bar'}};
end

load('lkcmap2')
argsin = setdefault(argsin,'colormap',lkcmap2)

if ~exist('xlabels','var') || isempty(xlabels)
    xlabels = cellcat('Predictor ',cellstr(num2str([1:size(plsmdl.XL,1)]')),'',0);
end

if ~exist('ylabels','var') || isempty(ylabels)
    ylabels = cellcat('Response ',cellstr(num2str([1:size(plsmdl.YL,1)]')),'',0);
    
end

for i = 1:length(coeffplots)
    if strcmpi(coeffplots{i}{2},'posture')
        plotwidth(i) = 0.5;
    elseif strcmpi(coeffplots{i}{2},'embody') & size(plsmdl.(coeffplots{i}{1}),1) > length(unique(EasyParse(argsin,'embodyatlas')))
        plotwidth(i) = 2;
    elseif strcmpi(coeffplots{i}{2},'embody')
        plotwidth(i) = 0.5;
    elseif strcmpi(coeffplots{i}{2},'embody_top')
        plotwidth(i) = 1.5;
    else
        plotwidth(i) = 1;
    end
end
plotwidth = (100./sum(plotwidth))*plotwidth;
plotwidth = mat2cell(plotwidth,1,ones(1,length(plotwidth)));

ncomps = length(whichcomps);

if isempty(pindx)
    p.pack('v',repmat({1/ncomps},1,ncomps))
    for i = 1:ncomps
        p(i).pack('h',{1/5 4/5})
    end
else
    p(pindx{:}).pack('v',repmat({1/ncomps},1,ncomps))
    for i = 1:ncomps
        p(pindx{:},i).pack('h',{1/3 2/3})
    end
end


for q = whichcomps
    p(pindx{:},find(whichcomps==q),1).select()
    nicecorrplot(plsmdl.XS(:,q),plsmdl.YS(:,q),{['Latent component ' num2str(q) ' - ' datasetlabels{1}],...
        ['Latent component ' num2str(q) ' - ' datasetlabels{2}]},'type','pearson','externalp',plsmdl.pperm(q));
    FixAxes(gca,20)
    
    p(pindx{:},find(whichcomps==q),2).pack('h',plotwidth);
    
    for i = 1:length(coeffplots)
        if contains(coeffplots{i}{1},'X')
            labs = xlabels;
        else
            labs = ylabels;
        end
        switch coeffplots{i}{2}
            case 'bar'
                p(pindx{:},find(whichcomps==q),2,i).select()
                l = lines;
                b = bar(plsmdl.(coeffplots{i}{1})(:,q),'facecolor',palecol(l(i,:)));
                hold on
                b2 = bar(plsmdl.(coeffplots{i}{1})(:,q).*(plsmdl.([coeffplots{i}{1} '_bootp'])(:,q)<0.05),'facecolor',l(1,:));
                if isfield(plsmdl,[coeffplots{i}{1} '_boot'])
                    e = errorbar(plsmdl.(coeffplots{i}{1})(:,q),1.96*std(plsmdl.([coeffplots{i}{1} '_boot'])(:,q,:),[],3),...
                        'LineStyle','none','LineWidth',2,'Color','k');
                end
                set(gca,'XTick',[1:size(plsmdl.(coeffplots{i}{1}),1)],'XTickLabel',labs)
                ylabel('Loading');
                xtickangle(90)
                FixAxes(gca,14)
                
                %t = p(pindx{:},find(whichcomps==q),2,i).title(datasetlabels{i});
                %t.FontWeight = 'Bold'; t.FontSize = 20;
            case 'barh'
                p(pindx{:},find(whichcomps==q),2,i).select()
                l = lines;
                b = barh(plsmdl.(coeffplots{i}{1})(:,q),'facecolor',palecol(l(i,:)));
                hold on
                b2 = barh(plsmdl.(coeffplots{i}{1})(:,q).*(plsmdl.([coeffplots{i}{1} '_bootp'])(:,q)<0.05),'facecolor',l(i,:));
                if isfield(plsmdl,[coeffplots{i}{1} '_boot'])
                    e = errorbar(plsmdl.(coeffplots{i}{1})(:,q),1:length(plsmdl.(coeffplots{i}{1})(:,q)),1.96*std(plsmdl.([coeffplots{i}{1} '_boot'])(:,q,:),[],3),...
                        'horizontal','LineStyle','none','LineWidth',2,'Color','k');
                end
                set(gca,'YTick',[1:size(plsmdl.(coeffplots{i}{1}),1)],'YTickLabel',labs)
                xlabel('Loading');
                FixAxes(gca,14)
                
                %t = p(pindx{:},find(whichcomps==q),2,i).title(plotlabels{i});
                %t.FontWeight = 'Bold'; t.FontSize = 20;
            case 'wordcloud'
                f = figure;
                if isfield(plsmdl,([coeffplots{i}{1} '_bootp']))
                    bootp = plsmdl.([coeffplots{i}{1} '_bootp'])(:,q);
                else
                    bootp = zeros(size(plsmdl.(coeffplots{i}{1})(:,q)));
                end
                wc = wordcloud_bipolar(labs,plsmdl.(coeffplots{i}{1})(:,q),EasyParse(argsin,'colormap'),'specific',bootp<0.05);
                %ax = findobj('parent','f','type','axes');
                p(pindx{:},find(whichcomps==q),2,i).select(wc)
                close(f)
                
                %t = p(pindx{:},find(whichcomps==q),2,i).title(plotlabels{i});
                %t.FontWeight = 'Bold'; t.FontSize = 20;
            case 'embody'
                atlas = EasyParse(argsin,'embodyatlas');
                
                emos = [{'ANGER'}    {'CONTEMPT'}    {'DISGUST'}    {'ENVY'}    {'FEAR'}    {'HAPPINESS'}    {'PRIDE'}    {'SADNESS'}    {'SHAME'}    {'SURPRISE'}];
                
                if length(plsmdl.(coeffplots{i}{1})(:,q)) == length(unique(atlas))
                    p(pindx{:},find(whichcomps==q),2,i).select()
                    map = vect2vol(plsmdl.(coeffplots{i}{1})(:,q),atlas);
                    embody_plotmap(map)
                    set(gca,'YDir','reverse')
                    cbar = colorbar;
                    cbar.Label.String = 'Loading'; cbar.FontSize = 14;
                else
                    theseloads = reshape(plsmdl.(coeffplots{i}{1})(:,q),sum(unique(atlas)>0),[]);
                    if ncomps ==1
                        dims = numSubplots(size(theseloads,2),1.5); dims = fliplr(dims);
                    else
                        dims = numSubplots(size(theseloads,2),2.5); %dims = fliplr(dims);
                    end
                    
                    emos = repmat(emos,1,size(theseloads,2)/10);
                    p(pindx{:},find(whichcomps==q),2,i).pack(dims(1),dims(2))
                    for qq = 1:size(theseloads,2)
                        [i1,i2] = ind2sub(dims,qq);
                        p(pindx{:},find(whichcomps==q),2,i,i1,i2).select()
                        axs(qq) = gca;
                        map = vect2vol(theseloads(:,qq),atlas);
                        embody_plotmap(map)
                        set(gca,'YDir','reverse')
                        title(emos{qq});
                    end
                    Normalize_Clim(axs,1);
                    set(gca,'YDir','reverse')
                    cbar = colorbar;
                    cbar.Label.String = 'Loading'; cbar.FontSize = 14;
                end
                
                %t = p(pindx{:},find(whichcomps==q),2,i).title(plotlabels{i});
                %t.FontWeight = 'Bold'; t.FontSize = 20;
                
            case 'embody_top'
                
                atlas = EasyParse(argsin,'embodyatlas');
                
                theseloads = reshape(plsmdl.(coeffplots{i}{1})(:,q),sum(unique(atlas)>0),[]);
                theseloads_p = reshape(plsmdl.([coeffplots{i}{1} '_bootp'])(:,q),sum(unique(atlas)>0),[]);
                theseloads = theseloads.*(theseloads_p<0.05);
                
                emos = {'Anger','Contempt','Disgust','Envy','Fear','Happiness','Pride','Sadness','Shame','Surprise'};
                
                emolabs = [strcat(emos,' Activation') strcat(emos,' Deactivation')];
                
                %loadsbyemo = reshape(theseloads,sum(unique(atlas)>0),length(emos),[]); % activation
                %loadsbyemo = loadsbyemo.*(reshape(theseloads_p,sum(unique(atlas)>0),length(emos),[])<0.05);
                
                bestemos = sum(abs(theseloads),1);
                
                [~,sortbest] = sort(bestemos); sortbest = fliplr(sortbest);
                
                ntop = 6;
                
                if ncomps == 1
                    dims = numSubplots(ntop,1.5); dims = fliplr(dims);
                else
                    dims = numSubplots(ntop,2.5); %dims = fliplr(dims);
                end
                dims = fliplr(dims);
                
                %emos = repmat(emos,1,size(theseloads,2)/10);
                p(pindx{:},find(whichcomps==q),2,i).pack(dims(1),dims(2))
                for qq = 1:ntop
                    [i1,i2] = ind2sub(dims,qq);
                    p(pindx{:},find(whichcomps==q),2,i,i1,i2).select()
                    axs(qq) = gca;
                    map = vect2vol(theseloads(:,sortbest(qq)),atlas);
                    embody_plotmap(map)
                    set(gca,'YDir','reverse')
                    axis tight
                    title(emolabs{sortbest(qq)},'FontSize',14);
                end
                Normalize_Clim(axs,1);
                set(gca,'YDir','reverse')
                cbar = colorbar;
                cbar.Label.String = 'Loading'; cbar.FontSize = 14;
                
                %t = p(pindx{:},find(whichcomps==q),2,i).title(plotlabels{i});
                %t.FontWeight = 'Bold'; t.FontSize = 20; t.Position = t.Position+[0 0.05 0];
                
            case 'posture'
                %if CheckInput(argsin,'postureweights') % translate back to posture feature space in case the model was done on components
                %    posweights = EasyParse(varargin,'postureweights'); % this should be a table
                %else
                posweights = array2table(horz(plsmdl.(coeffplots{i}{1})(:,q)),'VariableNames',labs); % otherwise, it's assumed that the model weights are the full posture dataset
                %end
                if CheckInput(argsin,'postureweights')
                    posweights = ttimes(posweights,EasyParse(varargin,'postureweights'));
                else
                    posweights = ttimes(posweights,25);
                end
                
                p(pindx{:},find(whichcomps==q),2,i).select()
                postureplot(posweights)
                axis off
                
            case 'posture_continuum'
                if ~CheckInput(argsin,'continuum_nsteps')
                    nsteps = 6;
                else
                    nsteps = EasyParse(argsin,'continuum_nsteps');
                end
                %thesesteps = linspace(0,100,nsteps);
                
                meandat = nanmean(plsmdl.(coeffplots{i}{1}(1)),1);
                
                % construct the steps by taking the projections of the
                % percentiles onto the components
                
                thesesteps = linspace(100,0,nsteps);
                
                for ii = 1:nsteps
                    thisscr = plsmdl.([coeffplots{i}{1}(1) 'S'])(:,find(whichcomps==q));
                    thisindx = FindClosest(thisscr,prctile(thisscr,thesesteps(ii)),1);
                    
                    tmp1 = (plsmdl.(coeffplots{i}{1}(1))(thisindx,:)-nanmean(plsmdl.(coeffplots{i}{1}(1)),1));
                    tmp2 = plsmdl.([coeffplots{i}{1}(1) 'L'])(:,find(whichcomps==q));
                    thesesteps(ii) = (tmp1./norm(tmp1))*(tmp2./norm(tmp2));
                end

                thesesteps = linspace(0,100,nsteps);
                
                %thesesteps = linspace(-0.7,0.7,nsteps);

                
                p(pindx{:},find(whichcomps==q),2,i).pack('h',repmat({1/nsteps},1,nsteps))
                
                for ii = 1:nsteps
                    p(pindx{:},find(whichcomps==q),2,i,ii).select()
                    
                    %thisscr = plsmdl.XS(:,find(whichcomps==q));
                    thisscr = plsmdl.([coeffplots{i}{1}(1) 'S'])(:,find(whichcomps==q));
                    thisindx = FindClosest(thisscr,prctile(thisscr,thesesteps(ii)),1);
                    postureplot(plsmdl.(coeffplots{i}{1}(1))(thisindx,:))
                    
                    %postureplot(meandat+thesesteps(ii)*plsmdl.([coeffplots{i}{1}(1) 'L'])(:,find(whichcomps==q))')
                    
                    %title([num2str(thesesteps(ii)) 'th' newline 'percentile'],'FontSize',16)
                    view(-90,0)
                    ax(ii) = gca;
                    axis off
                end
                
                Normalize_Xlim(ax,0); Normalize_Ylim(ax,0); Normalize_Zlim(ax,0);
                p(pindx{:},find(whichcomps==q),2,i).de.margin = [10 10 5 5];
        end
        
    end
    
end

p(pindx{:}).margin = [20 30 5 5];
p.margintop = 10;
p.marginleft = 25;
%p(pindx{:},2).marginleft = 22;



