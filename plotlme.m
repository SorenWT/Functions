function plotlme(mdl_lme,var,labels,varargin)

l = lines;

argsin = varargin;
argsin = setdefault(argsin,'scatterclr',palecol(l(1,:),0.33));
scatterclr = EasyParse(argsin,'scatterclr');

scatter(mdl_lme.Variables.(var),mdl_lme.Variables.(mdl_lme.ResponseName),384./round(log(mdl_lme.NumObservations)),scatterclr,'filled');

plotregline(mdl_lme.Variables.(var),mdl_lme.Variables.(mdl_lme.ResponseName),'mdl',mdl_lme,'plotCI','on')

FixAxes(gca,20)
xlabel(labels{1}); ylabel(labels{2})