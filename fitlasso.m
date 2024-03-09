function lassomdl = fitlasso(X,Y,distr,varargin)

if nargin < 3
   distr = 'normal';
end

[coeffs,fitinfo] = lassoglm(X,Y,distr,varargin{:});

lassomdl.coeffs = [fitinfo.Intercept(fitinfo.IndexMinDeviance); coeffs(:,fitinfo.IndexMinDeviance)];
lassomdl.fitinfo = fitinfo;