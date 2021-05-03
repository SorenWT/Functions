function p = lasso_plot(lassomdl,Y,xlabels,ylabels,varargin)
% this function makes a standardized plot for a lasso model defined by
% lasso_perm

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
    coeffplots = {'bar'};
end

if ~exist('xlabels','var') || isempty(xlabels)
    xlabels = cellcat('Predictor ',cellstr(num2str([1:size(Y,1)]')),'',0);
end

if ~exist('ylabels','var') || isempty(ylabels)
    ylabels = 'Response';
end

for i = 1:length(coeffplots)
   if strcmpi(coeffplots{i},'posture')
      plotwidth(i) = 0.5;
   else
       plotwidth(i) = 1;
   end
end
plotwidth = (100./sum(plotwidth))*plotwidth;
plotwidth = mat2cell(plotwidth,1,ones(1,length(plotwidth)));

if isempty(pindx)
    p.pack('h',{1/3 2/3})
else
    p(pindx{:}).pack('h',{1/3 2/3}) 
end

p(pindx{:},1).select()
nicecorrplot(lassomdl.pred,Y,{['Predicted ' ylabels],ylabels},'type','pearson','externalp',lassomdl.pperm);

p(pindx{:},2).pack('h',plotwidth);

for i = 1:length(coeffplots)
    switch coeffplots{i}
        case 'bar'
            p(pindx{:},2,i).select()
            l = lines;
            b = bar(lassomdl.coeffs,'facecolor',palecol(l(1,:)));
            set(gca,'XTick',[1:length(lassomdl.coeffs)],'XTickLabel',xlabels)
            ylabel('Coefficient');
            xtickangle(90)
            FixAxes(gca,14)
        case 'barh'
            p(pindx{:},2,i).select()
            l = lines;
            b = barh(lassomdl.coeffs,'facecolor',palecol(l(1,:)));
            set(gca,'XTick',[1:length(lassomdl.coeffs)],'XTickLabel',xlabels)
            ylabel('Coefficient');
            FixAxes(gca,14)
        case 'wordcloud'
            f = figure;
            wordcloud(xlabels,abs(lassomdl.coeffs))
            ax = findobj('parent','f','type','axes');
            p(pindx{:},2,i).select(ax)
            close(f)
        case 'posture'
            if CheckInput(argsin,'postureweights') % translate back to posture feature space in case the model was done on components
                posweights = EasyParse(varargin,'postureweights'); % this should be a table
            else
                posweights = array2table(lassomdl.coeffs,'VariableNames',xlabels); % otherwise, it's assumed that the lasso model weights are the full posture weights
            end
%                     posweights = array2table(horz(lassomdl.coeffs),'VariableNames',labs); % otherwise, it's assumed that the model weights are the full posture dataset
%                 %end
%                 if CheckInput(argsin,'postureweights')
%                    posweights = ttimes(posweights,EasyParse(varargin,'postureweights')); 
%                 end

            p(pindx{:},2,i).select()
            posturescreen_plot_ang(posweights)
    end
end

p(pindx{:}).margin = [20 30 5 5];
p(pindx{:},2).marginleft = 22;



