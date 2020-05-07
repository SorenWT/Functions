function s = ft_statfun_partialcorrT(cfg, dat,design)

if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end
if ~isfield(cfg, 'type'),              cfg.type           = 'Pearson'; end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
  ft_error('P-values can only be calculated if the test statistics are calculated.');
end
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
  ft_error('cfg.uvar should not exist for a correlation statistic');
end

[nsmpl,nrepl] = size(dat);
df = nrepl - 1;
if df<1
  ft_error('Insufficient error degrees of freedom for this analysis.')
end

if strcmp(cfg.computestat,'yes') % compute the statistic
  % calculate the correlation coefficient between the dependent variable and the predictor
  rho1 = corr(dat', design', 'type', cfg.type);
  rho2 = corr(cfg.partial', design', 'type', cfg.type);
  [z,p] = diffcorr(rho1,rho2,length(cfg.design),length(cfg.design));
  clear dat
  
  % convert correlation coefficient to t-statistic (for MCP correction): t^2 = DF*R^2 / (1-R^2)
  %tstat = rho*(sqrt(nrepl-2))./sqrt((1-rho.^2));
  
  s.stat = z; % store t values in s.stat variable for use with ft_statistics_montecarlo.m
  s.rho = rho1-rho2; % store r values in s.rho variable (these are the actual correlation coefficients)
  clear rho tstat
end

if strcmp(cfg.computecritval,'yes')
  % also compute the critical values
  s.df      = df-2;
%   if cfg.tail==-1
%     s.critval = tinv(cfg.alpha,df);
%   elseif  cfg.tail==0
%     s.critval = [tinv(cfg.alpha/2,df),tinv(1-cfg.alpha/2,df)];
%   elseif cfg.tail==1
%     s.critval = tinv(1-cfg.alpha,df);
%   end
s.critval = [-1.96 1.96];
end

if strcmp(cfg.computeprob,'yes')
  % also compute the p-values
  s.df      = df-2;
  s.prob = p;
%   if cfg.tail==-1
%     s.prob = tcdf(s.stat,s.df);
%   elseif  cfg.tail==0
%     s.prob = 2*tcdf(-abs(s.stat),s.df);
%   elseif cfg.tail==1
%     s.prob = 1-tcdf(s.stat,s.df);
%   end
end
