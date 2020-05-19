function [PLEout] = PLE_JF_EEG_wrapper(EEG,frange,varargin)

if nargin < 2
   frange = [0.5 50];
end

if nargin < 3
    argsin = {[]};
else
    argsin = varargin;
end

PLEout = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing power law exponent...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    PLEout(c) = JF_power_law(EEG.data(c,:),1/EEG.srate,frange(1),frange(2),argsin{:});
end