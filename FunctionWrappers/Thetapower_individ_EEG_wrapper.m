function [bpout] = Thetapower_individ_EEG_wrapper(EEG,lo_cutoff,norm_bandpass)

if nargin < 2
   norm_bandpass = 'no'; 
end

disp('Computing individualized theta power...')
disp('')

[psum,~,f] = restingIAF(EEG.data,EEG.nbchan,3,[1 30],EEG.srate,[5 15],11,5);

if ~isnan(psum.paf)
    pf = psum.paf;
else
    pf = psum.cog;
end

try
    iarange = [f(psum.iaw(1)) f(psum.iaw(2))];
    bpout = Bandpower_EEG_wrapper(EEG,[lo_cutoff iarange(1)],norm_bandpass);
catch
    bpout = NaN(1,EEG.nbchan);
end
