function [p,F,contrval] = coefTest_swt(mdl,H)

[p,F] = coefTest(mdl,H);

contrval = H*mdl.Coefficients.Estimate;