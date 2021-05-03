function p = pls_plot(plsmdl,ncomps,plotlabels,xlabels,ylabels,varargin)
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

if ~exist('xlabels','var') || isempty(xlabels)
    xlabels = cellcat('Predictor ',cellstr(num2str([1:size(plsmdl.XL,1)]')),'',0);
end

if ~exist('ylabels','var') || isempty(ylabels)
    ylabels = cellcat('Response ',cellstr(num2str([1:size(plsmdl.YL,1)]')),'',0);
    
end

for i = 1:length(coeffplots)
   if strcmpi(coeffplots{i}{2},'posture') || strcmpi(coeffplots{i}{2},'embody')
      plotwidth(i) = 0.5;
   else
       plotwidth(i) = 1;
   end
end
plotwidth = (100./sum(plotwidth))*plotwidth;
plotwidth = mat2cell(plotwidth,1,ones(1,length(plotwidth)));

if isempty(pindx)
    p.pack('v',repmat({1/ncomps},1,ncomps))
    for i = 1:ncomps
        p(i).pack('h',{1/3 2/3})
    end
else
    p(pindx{:}).pack('v',repmat({1/ncomps},1,ncomps))
    for i = 1:ncomps
        p(pindx{:},i).pack('h',{1/3 2/3})
    end
end


for q = 1:ncomps
    p(pindx{:},q,1).select()
    nicecorrplot(plsmdl.XS(:,q),plsmdl.YS(:,q),{['Latent component ' num2str(q) ' - ' plotlabels{1}],...
        ['Latent component ' num2str(q) ' - ' plotlabels{2}]},'type','pearson','externalp',plsmdl.pperm(q));
    
    p(pindx{:},q,2).pack('h',plotwidth);
    
    for i = 1:length(coeffplots)
        if contains(coeffplots{i}{1},'X')
            labs = xlabels;
        else
            labs = ylabels;
        end
        switch coeffplots{i}{2}
            case 'bar'
                p(pindx{:},q,2,i).select()
                l = lines;
                b = bar(plsmdl.(coeffplots{i}{1})(:,q),'facecolor',palecol(l(i,:)));
                hold on
                if isfield(plsmdl,[coeffplots{i}{1} '_boot'])
                    e = errorbar(plsmdl.(coeffplots{i}{1})(:,q),1.96*std(plsmdl.([coeffplots{i}{1} '_boot'])(:,q,:),[],3),...
                        'LineStyle','none','LineWidth',2,'Color','k');
                end
                set(gca,'XTick',[1:size(plsmdl.(coeffplots{i}{1}),1)],'XTickLabel',labs)
                ylabel('Loading');
                xtickangle(90)
                FixAxes(gca,14)
            case 'barh'
                p(pindx{:},q,2,i).select()
                l = lines;
                b = barh(plsmdl.(coeffplots{i}{1}),'facecolor',palecol(l(i,:)));
                hold on
                if isfield(plsmdl,[coeffplots{i}{1} '_boot'])
                    e = errorbar(plsmdl.(coeffplots{i}{1})(:,q),1:length(plsmdl.(coeffplots{i}{1})(:,q)),1.96*std(plsmdl.([coeffplots{i}{1} '_boot'])(:,q,:),[],3),...
                        'horizontal','LineStyle','none','LineWidth',2,'Color','k');
                end
                set(gca,'XTick',[1:size(plsmdl.(coeffplots{i}{1}),1)],'XTickLabel',labs)
                ylabel('Loading');
                FixAxes(gca,14)
            case 'wordcloud'
                f = figure;
                wordcloud(labs,abs(plsmdl.(coeffplots{i}{1})))
                ax = findobj('parent','f','type','axes');
                p(pindx{:},q,2,i).select(ax)
                close(f)
            case 'embody'
                atlas = EasyParse(argsin,'embodyatlas');
                p(pindx{:},q,2,i).select()
                
                map = vect2vol(plsmdl.(coeffplots{i}{1})(:,q),atlas);
                embody_plotmap(map)
                set(gca,'YDir','reverse')
                cbar = colorbar;
                cbar.Label.String = 'Loading'; cbar.FontSize = 14;
                
            case 'posture'
                %if CheckInput(argsin,'postureweights') % translate back to posture feature space in case the model was done on components
                %    posweights = EasyParse(varargin,'postureweights'); % this should be a table
                %else
                    posweights = array2table(horz(plsmdl.(coeffplots{i}{1})(:,q)),'VariableNames',labs); % otherwise, it's assumed that the model weights are the full posture dataset
                %end
                if CheckInput(argsin,'postureweights')
                   posweights = ttimes(posweights,EasyParse(varargin,'postureweights')); 
                end
                
                p(pindx{:},q,2,i).select()
                posturescreen_plot_ang(posweights)
        end
    end
    
end

p(pindx{:}).margin = [20 30 5 5];
%p(pindx{:},2).marginleft = 22;



