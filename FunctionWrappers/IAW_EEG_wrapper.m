function [iawout] = IAW_EEG_wrapper(EEG)

disp('Computing alpha peak frequency...')
disp('')

[psum,~,f] = restingIAF(EEG.data,EEG.nbchan,3,[1 40],EEG.srate,[5 15],11,5);

if ~isnan(psum.paf)
    pf = psum.paf;
else
    pf = psum.cog;
end

try
iaw = f(psum.iaw(2))-f(psum.iaw(1));

iawout = repmat(iaw,1,EEG.nbchan);
catch
    iawout = repmat(NaN,1,EEG.nbchan);
end