function hout = plotAdded_swt(mdl,coef,opts)

if nargin < 3
opts = struct;    
end

opts = setdefault(opts,'DummyVarCoding','full');

if isa(mdl,'LinearMixedModel')
    oldmdl = mdl;
    
    newformula = char(oldmdl.Formula); 
    newformula = erase(newformula,' ');
    % remove the random effects
    newformula = regexprep(newformula, '\(.*?\)', '');
    newformula(end) = []; % only will work with one random effect, fix later
    
    mdl = fitlm(oldmdl.Variables,newformula,'DummyVarCoding',opts.DummyVarCoding);
    %mdl.Coefficients = oldmdl.Coefficients;
    %mdl = oldmdl; 
%    mdl.Coefficients = oldmdl.Coefficients;
end


hold on
hout = plotAdded(mdl,coef);


x = hout(1).XData; y = hout(1).YData; 
xlab = get(gca,'XLabel'); xlab = xlab.String;
ylab = get(gca,'YLabel'); ylab = ylab.String;


clf
ci = coefCI(mdl); ci = ci(coef,:);
nicecorrplot_v2(x,y,{xlab,ylab},'Plot','off','regci',ci)


% hout(1).Color = palecol(hout(1).Color,0.3);
% hout(1).Marker = 'o'; hout(1).MarkerFaceColor = hout(1).Color;
% hout(2).LineWidth = 3;
% hout(3).LineWidth = 1.5;
% FixAxes(gca,16)
legend({'Adjusted data','95% conf. bounds','Fit'})

if ischar(coef)
    [corrrho,corrp] = partialcorr(mdl.Variables.(coef),mdl.Variables.(mdl.ResponseName),...
        mdl.Variables{:,setdiff(mdl.VariableNames,{coef mdl.ResponseName})},'rows','pairwise');
    corrp = mdl.Coefficients.pValue(find(strcmpi(mdl.PredictorNames,coef))+1);
end

% ax = gca;
% pos = ax.Position;
% tb = annotation('textbox','String',{['r_{partial} = ' num2str(round(corrrho,3))];['p = ' num2str(round(corrp,3,'significant'))]},...
%     'FitBoxToText','on','LineStyle','none','FontSize',14);
% tbsize = get(tb,'Position');
% delete(tb)
% edges = [pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)]; %left bottom right top
% tb = annotation('textbox','Position',[edges(3)-tbsize(3)-pos(3)*0.05 edges(4)-tbsize(4)-0.05*pos(4) tbsize(3) tbsize(4)],...
%     'String',{['r_{partial} = ' num2str(round(corrrho,3))];['p = ' num2str(round(corrp,3,'significant'))]},'FitBoxToText','on','LineStyle','none','FontSize',14);


end
