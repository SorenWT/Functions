function [F,p]=ftest_modelcompare(mdl1,mdl2)
% first model must be the simpler model

if mdl1.DFE<mdl2.DFE
    error('1st model should be the simpler one!')
end

F = ((mdl1.SSE-mdl2.SSE)/(mdl1.DFE-mdl2.DFE))/(mdl2.SSE/mdl2.DFE);

p = 1-fcdf(F,mdl1.DFE,mdl2.DFE);