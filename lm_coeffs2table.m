function [tblout]=lm_coeffs2table(lm,varnames)

if nargin < 2
   varnames = lm.Coefficients.Properties.RowNames(2:end);
end

tblout = lm.Coefficients;

tblout{:,:} = niceround(tblout{:,:});

tblout.rnames = [{'Intercept'} horz(varnames)]';

tblout = tblout(:,[end 1:end-1]);
tblout.Properties.RowNames = {};
%tblout.Properties.RowNames = [{'Intercept'} horz(varnames)];

end
