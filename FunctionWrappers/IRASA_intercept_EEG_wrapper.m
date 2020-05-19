function [IrasaOut] = IRASA_intercept_EEG_wrapper(spec,frange)

disp(' ')
disp('Computing PLE with IRASA...')

if nargin < 2
   frange = [0.5 50]; 
end

tmp = amri_sig_plawfit(spec,frange);
IrasaOut = tmp.Cons;
IrasaOut = IrasaOut';
