function [MFOut] = MF_EEG_wrapper(EEG,frange,oscifrac)

if nargin < 2
    frange = [0.5 50];
end

if isfield(EEG,'nbchan')
    
    MFOut = zeros(1,EEG.nbchan);
    
    disp(' ')
    disp('Computing median frequency...')
    
    for c = 1:EEG.nbchan
        fprintf([num2str(c) ' ']);
        MFOut(c) = medfreq(EEG.data(c,:),EEG.srate,frange);
    end
    
else
    
    if nargin < 3
        oscifrac = 'osci';
    end
    
    for c = 1:size(EEG.(oscifrac),2)
        fprintf([num2str(c) ' ']);
        MFOut(c) = medfreq(EEG.(oscifrac),EEG.freq,frange);
    end
end