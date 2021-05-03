function p = pca_plot(pcamdl,xlabels,varargin)
% this function makes a standard plot for a pca model defined by pca_parallel

argsin = varargin;

if CheckInput(argsin,'panel')
    p = EasyParse(argsin,'panel');
else
    figure
    p = panel('no-manage-font');
end

if CheckInput(argsin,'panelindx')
    pindx = EasyParse(argsin,'panelindx');
else
    pindx = {};
end

if CheckInput(argsin,'plottypes')
    plottypes = EasyParse(argsin,'plottypes');
else
    plottypes = {'bar'};
end

if CheckInput(argsin,'numplots')
    numplots = EasyParse(argsin,'numplots');
else
    numplots = pcamdl.ncomps;
end

if ~exist('xlabels','var') || isempty(xlabels)
    xlabels = cellcat('Predictor ',cellstr(num2str([1:size(pcamdl.XL,1)]')),'',0);
end

if CheckInput(argsin,'fontsize')
    fsize = EasyParse(argsin,'fontsize');
else
    fsize = 14;
end


if isempty(pindx)
    p.pack('h',{0.2 0.8});
else
    p(pindx{:}).pack('v',repmat({1/ncomps},1,ncomps))
    for i = 1:ncomps
        p(pindx{:},i).pack('h',{1/3 2/3})
    end
end


for i = 1:length(plottypes)
   if strcmpi(plottypes{i},'posture')
      plotwidth(i) = 0.75;
   else
       plotwidth(i) = 1;
   end
end
plotwidth = (100./sum(plotwidth))*plotwidth;
plotwidth = mat2cell(plotwidth,1,ones(1,length(plotwidth)));



p(pindx{:},1).pack('v',{1/2 1/2})

p(pindx{:},1,1).select()
l = lines;
plot(1:pcamdl.ncomps,pcamdl.explained(1:pcamdl.ncomps),'LineWidth',5,'color',l(1,:));
hold on
plot(1:length(pcamdl.explained),pcamdl.explained(1:end),'LineWidth',2,'color',l(1,:),'linestyle','--','HandleVisibility','off')
stdshade(1:length(pcamdl.explained),pcamdl.expl_perm','k',0.3,2,'prctileci');
xl = get(gca,'XLim');
line(xl,[100/length(xlabels) 100/length(xlabels)],'LineWidth',2,'color','r')
set(gca,'XLim',xl);
legend({'Eigenvalues','Parallel analysis eigenvalues (95% CI)','Kaiser threshold'})
FixAxes(gca,fsize*1.2)
xlabel('Component')
ylabel('Explained variance (%)')

p(pindx{:},1,1).marginbottom = 25;

p(pindx{:},1,2).pack()
p(pindx{:},1,2).pack({[0 0 1 1]})

p(pindx{:},1,2,1).select()
l = lines;
plot(cumsum(pcamdl.explained),'LineWidth',3,'color',l(2,:));
hold on
yl = get(gca,'YLim');
line([pcamdl.ncomps pcamdl.ncomps],yl,'color','k','LineWidth',1.5)
line([pcamdl.kaiser pcamdl.kaiser],yl,'color','k','LineWidth',1.5)
set(gca,'YLim',yl)
FixAxes(gca,fsize*1.2)
xlabel('Component')
ylabel('Cumulative explained variance')
ax1 = gca;
p(pindx{:},1,2,2).select()
set(gca,'color','none','YColor','none','FontSize',fsize*1.2);
if pcamdl.ncomps < pcamdl.kaiser
    set(gca,'XTick',[pcamdl.ncomps pcamdl.kaiser],...
        'XTickLabel',{'Parallel analysis','Kaiser criterion'})
elseif pcamdl.ncomps == pcamdl.kaiser
    set(gca,'XTick',[pcamdl.ncomps],...
        'XTickLabel',{'Parallel & kaiser'})
else
    set(gca,'XTick',[pcamdl.kaiser pcamdl.ncomps],...
        'XTickLabel',{'Kaiser criterion','Parallel analysis'})
end
xtickangle(45)
set(gca,'XLim',get(p(pindx{:},1,2,1).axis,'XLim'))

p(pindx{:},2).pack('h',plotwidth);

for q = 1:length(plottypes)
    
    plottype = plottypes{q};
    plotdim = numSubplots(numplots);
    
    p(pindx{:},2,q).pack(plotdim(1),plotdim(2));
    
    for i = 1:numplots
        [i1,i2] = ind2sub(plotdim,i);
        switch plottype
            case 'bar'
                p(pindx{:},2,q,i1,i2).select()
                b = bar(pcamdl.loads(:,i),'facecolor',palecol(l(1,:)));
                hold on
                if isfield(pcamdl,'loads_boot')
                    e = errorbar(pcamdl.loads(:,i),1.96*std(pcamdl.loads_boot(:,i,:),[],3),...
                        'LineStyle','none','LineWidth',2,'Color','k');
                end
                set(gca,'XTick',[1:size(pcamdl.loads,1)],'XTickLabel',xlabels)
                ylabel('Loading');
                xtickangle(90)
                FixAxes(gca,fsize)
            case 'barh'
                p(pindx{:},2,q,i1,i2).select()
                b = barh(pcamdl.loads,'facecolor',palecol(l(1,:)));
                hold on
                if isfield(pcamdl,'loads_boot')
                    e = errorbar(pcamdl.loads(:,i),1:length(pcamdl.loads(:,i)),1.96*std(pcamdl.loads_boot(:,i,:),[],3),...
                        'horizontal','LineStyle','none','LineWidth',2,'Color','k');
                end
                set(gca,'XTick',[1:size(pcamdl.loads,1)],'XTickLabel',xlabels)
                ylabel('Loading');
                FixAxes(gca,fsize)
            case 'wordcloud'
                f = figure;
                wordcloud(xlabels,abs(pcamdl.loads))
                ax = findobj('parent','f','type','axes');
                p(pindx{:},2,q,i1,i2).select(ax)
                close(f)
            case 'posture'
                p(pindx{:},2,q,i1,i2).select()
                if CheckInput(argsin,'posturenames') % translate back to posture feature space in case the model was done on components
                    posnames = EasyParse(varargin,'posturenames'); % this should be a table
                else
                    posnames = xlabels;
                end
                
                if CheckInput(argsin,'posturescale')
                    posscale = EasyParse(varargin,'posturescale');
                else
                    posscale = 45;
                end
                
                posweights = array2table(pcamdl.loads(:,i)'.*posscale,'VariableNames',posnames); % otherwise, it's assumed that the lasso model weights are the full posture weights
                
                posturescreen_plot_ang(posweights)
        end
        title(['Component ' num2str(i) ' loadings'])
    end
    
end

p(pindx{:}).margin = [20 30 5 10];
p(pindx{:},2).marginleft = 22;
p(pindx{:},2).de.margin = [20 20 5 5];



