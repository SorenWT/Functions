function [teststat, h, stats] = fast_signrank(x,y,varargin)
%SIGNRANK Wilcoxon signed rank test for zero median.
%   P = SIGNRANK(X) performs a two-sided signed rank test of the hypothesis
%   that the data in the vector X come from a distribution whose median
%   (and mean, if it exists) is zero, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("median is zero") is true.
%   Small values of P cast doubt on the validity of the null hypothesis.
%   The data are assumed to come from a continuous distribution, symmetric
%   about its median.
%
%   P = SIGNRANK(X,M) performs a two-sided test of the hypothesis that the
%   data in the vector X come from a distribution whose median is M.  M
%   must be a scalar.
%
%   P = SIGNRANK(X,Y) performs a paired, two-sided test of the hypothesis
%   that the difference between the matched samples in the vectors X and Y
%   comes from a distribution whose median is zero.  The differences X-Y
%   are assumed to come from a continuous distribution, symmetric about its
%   median.  X and Y must be the same length.  The two-sided p-value is
%   computed by doubling the most significant one-sided value.
%
%   SIGNRANK treats NaNs in X or Y as missing values, and removes them.
%
%   [P,H] = SIGNRANK(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("median is zero") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = SIGNRANK(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = SIGNRANK(...,'method',METHOD) computes the p-value using an
%   exact algorithm if METHOD is 'exact', or a normal approximation if
%   METHOD is 'approximate'.  The default is to use an exact method for
%   small samples.
%
%   [P,H] = SIGNRANK(...,'tail',TAIL) performs the test against the
%   alternative hypothesis specified by TAIL:
%       'both'  -- "median is not zero (or M)" (two-tailed test, default)
%       'right' -- "median is greater than zero (or M)" (right-tailed test)
%       'left'  -- "median is less than zero (or M)" (left-tailed test)
%   TAIL must be a single string.
%
%   [P,H,STATS] = SIGNRANK(...) returns STATS, a structure with one or two
%   fields.  The field 'signedrank' contains the value of the signed rank
%   statistic for positive values in X, X-M, or X-Y. If P is calculated
%   using a normal approximation, then the field 'zval' contains the value
%   of the normal (Z) statistic.
%
%   See also SIGNTEST, RANKSUM, TTEST, ZTEST.

%   Copyright 2010-2012 The MathWorks, Inc.

%   For the two-sample case, SIGNRANK uses a tolerance based on the
%   values EPSD=EPS(X)+EPS(Y). Any pair of values of D=X-Y that differ by
%   no more than the sum of their two EPSD values are treated as ties.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.



% Check most of the inputs now
% if nargin > 2
%     [varargin{:}] = convertStringsToChars(varargin{:});
% end
% 
% alpha = 0.05;
% if nargin>2 && isnumeric(varargin{1})
%    % Grandfathered syntax:  signrank(x,y,alpha)
%    alpha = varargin{1};
%    varargin(1) = [];
% end
% oknames = {'alpha' 'method' 'tail'};
% dflts   = {alpha   ''  'both'};
% [alpha,method,tail] = internal.stats.parseArgs(oknames,dflts,varargin{:});
% 
% if ~isscalar(alpha)
%    error(message('stats:signrank:BadAlpha'));
% end
% if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
%    error(message('stats:signrank:BadAlpha'));
% end
% 
% onesample = false;
% if nargin < 2 || isempty(y)
%     y = zeros(size(x));
%     onesample = true;
% elseif isscalar(y)
% 	y = repmat(y, size(x));
% end
% 
% if ~isvector(x) || ~isvector(y)
% 	error(message('stats:signrank:InvalidData'));
% elseif numel(x) ~= numel(y)
% 	error(message('stats:signrank:InputSizeMismatch'));
% end

method = 'approximate';
alpha = 0.05;
tail = 'both';
onesample = 0;

diffxy = x(:) - y(:);
if onesample
	epsdiff = zeros(size(x(:)));
else
	epsdiff = eps(x(:)) + eps(y(:));
end

% Remove missing data
t = isnan(diffxy);
diffxy(t) = [];
epsdiff(t) = [];
if isempty(diffxy)
	error(message('stats:signrank:NotEnoughData'));
end

t = (abs(diffxy) <= epsdiff);
diffxy(t) = [];
epsdiff(t) = [];

n = length(diffxy);

if (n == 0)         % degenerate case, all ties
	p = 1;
    teststat = 0;
	if (nargout > 1)
		h = 0;
		if (nargout > 2)
			stats.signedrank = 0;
		end
	end
	return
end

% Now deal with the method argument
if isempty(method)
	if n<=15
		method = 'exact';
	else
		method = 'approximate';
	end
elseif strcmpi(method, 'oldexact')
	method = 'oldexact';
else   % method not recognized, throw error
	method = internal.stats.getParamVal(method,{'exact' 'approximate'},'''method''');
end

% Check the tail argument
tail = internal.stats.getParamVal(tail,{'both' 'right' 'left'},'''tail''');


% Calculations for Sign Rank Test

% Find positive differences and ranks of absolute differences
iPos = (diffxy>0);
[tie_rank, tieadj] = tiedrank(abs(diffxy),0,0,epsdiff);

% Compute signed rank statistic (most extreme version)
w = sum(tie_rank(iPos));
w2 = (n*(n+1)/2) - w;


if strcmpi(method, 'approximate')

	switch tail
		case 'both'
			z = (w-n*(n+1)/4) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
			%p = 2*normcdf(-abs(z),0,1);
			
		case 'right'
			z = (w-n*(n+1)/4 - 0.5) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
			%p = normcdf(-z, 0, 1);
			
		case 'left'
			z = (w-n*(n+1)/4 + 0.5) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
			%p = normcdf(z, 0, 1);
	end
	
	if (nargout > 2)
		stats.zval = z;
	end
	
elseif strcmpi(method, 'oldexact')
	% Enumerates all possibilities and does not adjust for ties
	allposs = (ff2n(n))';
	idx = (1:n)';
	idx = idx(:,ones(2.^n,1));
	pranks = sum(allposs.*idx,1);
	
	switch tail
		case 'both'
			in_tails =  sum(pranks <= (min(w,w2)))  +  ...
				sum(pranks >= max(w,w2));
			% Avoid p>1 if w is in the middle and middle is double-counted
			p = min(1, in_tails./(2.^n));
			
		case 'right'
			in_tail =  sum(pranks >= w);
			p =  in_tail ./ (2.^n);
			
		case 'left'
			in_tail = sum(pranks <= w);
			p =  in_tail ./ (2.^n);
	end
	
else   % strcmpi(method, 'exact')
	[p, many_w_p] = statsrexact(tie_rank, w);   % probability of smaller tail

	switch tail
		case 'both'
			p = min(1, 2*p);   % two-sided, don't double-count the middle value
			
		case 'right'
			if w < w2   % right tail is larger
				p =  1 - p + many_w_p(end,2);
			end
			
		case 'left'
			if w > w2   % left tail is larger
				p =  1 - p + many_w_p(end,2);
			end
	end
	
end   % end of conditional on method


if nargout > 1
    h = (p<=alpha);
    if (nargout > 2)
        stats.signedrank = w;
    end
end

if strcmpi(method,'approximate')
   teststat = z;
else
    teststat = w;
end