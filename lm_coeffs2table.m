function [tblout]=lm_coeffs2table(lm,varnames)

tblout = lm.Coefficients;

tblout{:,:} = niceround(tblout{:,:});

tblout.Properties.RowNames = [{'Intercept'} varnames];

end

function rounddat = niceround(dat)

if any(size(dat)>1)
    rounddat = arrayfun(@niceround,dat);
else
    
    if abs(dat)>10
        rounddat = round(dat,1,'decimal');
    elseif abs(dat)>1
        rounddat = round(dat,2,'decimal');
    else
        rounddat = round(dat,2,'significant');
    end
end

end