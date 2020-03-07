function [delR,R] = ErrorProp(formula,varargin)

nvars = length(varargin)/2;
allvars = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p'};
vars = allvars(1:nvars);
varcommas = cellcat(',',vars,'',1);
varcommas{end} = varcommas{end}(1:end-1);
delvars = cellcat('del',vars,'',0);
delvarcommas = cellcat(',',delvars,'',1);
delvarcommas{end} = delvarcommas{end}(1:end-1);

eval(['syms y(' varcommas{:} ')'])
eval(['syms z(' varcommas{:} ',' delvarcommas{:} ')'])

eval(['y(' varcommas{:} ') = eval(formula);']);

sqrtstring = [];
for q = 1:length(vars)
    sqrtstring = cat(2,sqrtstring,['(diff(y,' vars{q} ')^2)*(' delvars{q} '^2)+']);
end

sqrtstring(end) = [];

eval(['z(' varcommas{:} ',' delvarcommas{:} ') = sqrt(' sqrtstring ');']);

R = double(y(varargin{1:2:end}));
delR = double(z(varargin{1:2:end},varargin{2:2:end}));