function [s, cfg] = ft_statfun_reverse_indepsamplesF(cfg, dat, design)
% SWT custom function for doing an F-test where the imaging data gives the
% independent variable, while the design gives the dependent variable
%
% FT_STATFUN_INDEPSAMPLESF calculates the independent samples F-statistic on the
% biological data in dat (the dependent variable), using the information on the
% independent variable (ivar) in design.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_indepsamplesF'
%
% Configuration options
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification
%   cfg.ivar  = row number of the design that contains the labels of the conditions that must be
%               compared (default=1). The labels range from 1 to the number of conditions.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2006, Eric Maris
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% set the defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    ft_error('P-values can only be calculated if the test statistics are calculated.');
end
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    ft_error('cfg.uvar should not exist for an independent samples statistic');
end

% this assumes that dat has the same number of conditions at each channel
ncond = length(unique(dat(1,:)));
nrepl=0;
for condindx=1:ncond
    nrepl=nrepl+length(find(dat(1,:)==condindx));
end
if nrepl<size(design,2)
    ft_error('Invalid specification of the independent variable in the design array.');
end
if nrepl<=ncond
    ft_error('The must be more trials/subjects than levels of the independent variable.');
end
dfnum = ncond - 1;
dfdenom = nrepl - ncond;

nsmpls = size(dat,1);

design = vert(design);

if strcmp(cfg.computestat, 'yes')
    for i = 1:nsmpls
        % compute the statistic
        nobspercell = zeros(1,ncond);
        avgs = zeros(1,ncond);
        pooledvar = zeros(1);
        for condindx=1:ncond
            sel=find(dat(i,:)==condindx);
            nobspercell(condindx)=length(sel);
            avgs(condindx)=mean(design(sel));
            pooledvar = pooledvar + nobspercell(condindx)*var(design(sel),1);
        end
        pooledvar = pooledvar/dfdenom;
        globalavg = mean(design);
        mseffect = ((avgs-repmat(globalavg,1,ncond)).^2)*nobspercell'./dfnum;
        s.stat(i) = mseffect./pooledvar;
    end
    s.stat = vert(s.stat);
end


if strcmp(cfg.computecritval,'yes')
    % also compute the critical values
    s.dfnum   = dfnum;
    s.dfdenom = dfdenom;
    if cfg.tail==-1
        ft_error('For an independent samples F-statistic, it does not make sense to calculate a left tail critical value.');
    end
    if cfg.tail==0
        ft_error('For an independent samples F-statistic, it does not make sense to calculate a two-sided critical value.');
    end
    if cfg.tail==1
        s.critval = finv(1-cfg.alpha,s.dfnum,s.dfdenom);
    end
end

if strcmp(cfg.computeprob,'yes')
    % also compute the p-values
    s.dfnum   = dfnum;
    s.dfdenom = dfdenom;
    if cfg.tail==-1
        ft_error('For an independent samples F-statistic, it does not make sense to calculate a left tail p-value.');
    end
    if cfg.tail==0
        ft_error('For an independent samples F-statistic, it does not make sense to calculate a two-sided p-value.');
    end
    if cfg.tail==1
        s.prob = 1-fcdf(s.stat,s.dfnum,s.dfdenom);
    end
end

