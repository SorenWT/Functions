function [delR,R] = ErrorProp3(formula,varargin)

nvars = length(varargin)/2;
allvars = {'a','b','c','d','e','f','g','h','i','j','k'};
vars = allvars(1:nvars);
varcommas = cellcat(',',vars,'',1);
varcommas{end} = varcommas{end}(1:end-1);
delvars = cellcat('del',vars,'',0);
delvarcommas = cellcat(',',delvars,'',1);
delvarcommas{end} = delvarcommas{end}(1:end-1);

eval(['syms f(' varcommas{:} ')'])
eval(['syms g(' varcommas{:} ',' delvarcommas{:} ')'])

eval(['f(' varcommas{:} ') = eval(formula);']);

sqrtstring = [];
for z = 1:length(vars)
    sqrtstring = cat(2,sqrtstring,['(diff(f,' vars{z} ')^2)*(' delvars{z} '^2)+']);
end

sqrtstring(end) = [];

eval(['g(' varcommas{:} ',' delvarcommas{:} ') = sqrt(' sqrtstring ');']);

R = double(f(varargin{1:2:end}));
delR = double(g(varargin{1:2:end},varargin{2:2:end}));




%R = form(x,y);

% if strcmpi(formType,'add')
%     %syms delform(A,delx,B,dely)
%     %delform(A,delx,B,dely) = ;
%     delR = sqrt((A^2)*(delx^2) + (B^2)*(dely^2));
% elseif strcmpi(formType,'multiply')
%     %    syms delform(A,delx,B,dely)
%     %delform(A,x,delx,B,y,dely,R) = abs(R)*sqrt((A^2)*((delx/x)^2) + (B^2)*((dely/y)^2));
%     delR = abs(R)*sqrt((A^2)*((delx/x)^2) + (B^2)*((dely/y)^2));
% end
