function [bp] = Bandpower_EEG_wrapper(EEG,frange,norm_bandpass)

if nargin < 3
   norm_bandpass = 'no'; 
end

if any(isnan(EEG.data))
   warning('Data contains NaNs - discontinuities may introduce errors in frequency analysis')
end

for c = 1:EEG.nbchan
   
    % handle NaNs if they were included to mark artefacts
   thisdat = EEG.data(c,:); thisdat = thisdat(~isnan(thisdat));
  
   [pxx,f] = pwelch(thisdat,[],[],2^nextpow2((3/2)*EEG.srate),EEG.srate);
   findx = intersect(find(f > frange(1)),find(f < frange(2)));
   
   %bp(c) = trapz(f(findx),pxx(findx))/numel(findx);

   bp(c) = (norm(pxx(findx))^2)./numel(findx);
    %bp(c) = bandpower(EEG.data(c,:),EEG.srate,frange); 
   if ~strcmpi(norm_bandpass,'no')
       allfindx = intersect(find(f > norm_bandpass(1)),find(f < norm_bandpass(2)));
       %bp(c) = bp(c)/(trapz(f(allfindx),pxx(allfindx))/numel(allfindx));
       bp(c) = bp(c)/((norm(pxx(allfindx))^2)./numel(allfindx));
        %bp(c) = bp(c)/bandpower(EEG.data(c,:),EEG.srate,norm_bandpass); 
   end
end
