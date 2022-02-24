function [bpout] = Alphapower_individ_EEG_wrapper(EEG,norm_bandpass,irasa)

if nargin < 2
    norm_bandpass = 'no';
end

if nargin < 3
    irasa = 'no';
end

disp('Computing individualized alpha power...')
disp('')

[psum,~,f] = restingIAF(EEG.data,EEG.nbchan,3,[1 30],round(EEG.srate),[5 15],11,5);

if ~isnan(psum.paf)
    pf = psum.paf;
else
    pf = psum.cog;
end

%try
    iarange = [f(psum.iaw(1)) f(psum.iaw(2))];
    if strcmpi(irasa,'no')
        bpout = Bandpower_EEG_wrapper(EEG,iarange,norm_bandpass);
    else
        EEG = parload(fullfile(EEG.filepath,[EEG.filename '_IRASA_specs.mat']),'EEG');
        bpout = IRASAPower_EEG_wrapper(EEG,'osci',iarange,norm_bandpass);
    end
%catch
%    bpout = NaN(1,EEG.nbchan);
%end

