function [p,result] = corrplot_hist(X,Y,varargin)

argsin = varargin;

if ~CheckInput(argsin,'varnames')
    varnames = {'X','Y'};
else
    varnames = EasyParse(argsin,'varnames');
end

if CheckInput(argsin,'panel')
    p = EasyParse(argsin,'panel');
    pindx = EasyParse(argsin,'panelindx');
else
    f = figure;
    p = panel('no-manage-font');
    p.pack(1)
    pindx = {1};
    
    p.margin = [20 20 10 10];
    p.de.margin = [15 10 5 5];
end

argsin = setdefault(argsin,'removeoutliers','off');

if EasyParse(argsin,'removeoutliers','on')
   outsindx = isoutlier(X) | isoutlier(Y);
   X = X.*nanmask(~outsindx); Y = Y.*nanmask(~outsindx);
end

argsin = removeargs(argsin,{'varnames','panel','panelindx','removeoutliers'});

[result.r,result.p,result.ci,result.n] = corrswt(X,Y,argsin{:});

% 
% result.mdl{1} = fitlm([X,vert(grp),X.*vert(grp)],Y,'VarNames',{varnames{1},varnames{3},[varnames{1} '*' varnames{3}],varnames{2}});
% 
% result.mdl{2} = fitlm([Y,vert(grp),Y.*vert(grp)],X,'VarNames',{varnames{2},varnames{3},[varnames{2} '*' varnames{3}],varnames{1}});

% do the actual plotting

l = lines;

tmptbl = array2table([X,Y],'VariableNames',varnames);

%s = scatterhistogram(tmptbl,varnames{1},varnames{2},'GroupVariable',varnames{3},...
%    'MarkerSize',36,'Color',palecol(l(1:2,:)),'HistogramDisplayStyle','smooth')
scatterhist_panel(X,Y,p,pindx,'Kernel','on','MarkerSize',24,'Marker','.',...
    'color',palecol(l(1,:)),'Location','NorthEast','Direction','out');
p(pindx{:},2,2).marginleft = 5;
p(pindx{:},2,1).select()
hold on
plotregline(X,Y,'color','r','linewidth',2)
%plotregline(X(grp==ungrp(1)),Y(grp==ungrp(1)),'color',l(1,:),'linewidth',2)
%plotregline(X(grp==ungrp(2)),Y(grp==ungrp(2)),'color',l(2,:),'linewidth',2)
%manlegend({'Overall estimate',grpnames{1},grpnames{2}},[0 0 0; l(1:2,:)])
FixAxes(gca,16)
set(gca,'XAxisLocation','bottom','YAxisLocation','left');
xlabel(varnames{1})
ylabel(varnames{2})

if ~CheckInput(argsin,'Type')
    rstring = '\rho';
else
    typ = EasyParse(argsin,'Type');
    
    if strcmpi(typ,'spearman')
        rstring = '\rho';
    elseif strcmpi(typ,'pearson')
        rstring = 'r';
    end
end

easyannotate({[rstring ' = ' num2str(round(result.r,3))];['p = ' num2str(round(result.p,3,'significant'))]},[],'Color','k')
%easyannotate({[rstring ' = ' num2str(round(result.(grpnames{1}).r,3))];['p = ' num2str(round(result.(grpnames{1}).p,3,'significant'))]},[],'Color',l(1,:))
%easyannotate({[rstring ' = ' num2str(round(result.(grpnames{2}).r,3))];['p = ' num2str(round(result.(grpnames{2}).p,3,'significant'))]},[],'Color',l(2,:))

%title(['Interaction p = ' num2str(round(result.mdl{1}.Coefficients.pValue(4),3,'significant'))])

%p(pindx{:},1,1).select()
%hold on
%xlim = get(gca,'XLim');
%if result.mdl{2}.Coefficients.pValue(3) < 0.05
%    ss = sigstar({[min(xlim) max(xlim)]},result.mdl{2}.Coefficients.pValue(3),0,24);
%    delete(ss(1))
%end
% 
% p(pindx{:},2,2).select()
% hold on
% xlim = get(gca,'XLim');
% if result.mdl{1}.Coefficients.pValue(3) < 0.05
%     ss = sigstar({[min(xlim) max(xlim)]},result.mdl{1}.Coefficients.pValue(3),0,24);
%     delete(ss(1));
% end

%title(['Difference z = ' num2str(round(moder.rdiff.z,3,'significant')) ', p = ' num2str(round(moder.rdiff.p,3,'significant'))])

%title(['Difference z = ' num2str(round(moder.rdiff.z,3,'significant')) ', p = ' num2str(round(moder.rdiff.p,3,'significant'))])

