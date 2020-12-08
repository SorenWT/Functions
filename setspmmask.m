function setspmmask(th)
% SETSPMMASK sets the threshold mask to the specified value or -Inf
%     if no value supplied. Particularly needed for perfusion data.

global defaults

% set the global value
if nargin == 1
     defaults.mask.thresh = th;
else
     defaults.mask.thresh = -Inf;
end

% clear the local copy of the global defaults.
clear defaults

% reload the global
global defaults

% check if it's been changed.
if defaults.mask.thresh == -Inf;
     fprintf('\nAnalysis masking threshold set to -Inf.\n')
else
     fprintf('\nProblem setting analysis masking threshold.\n')
end
end