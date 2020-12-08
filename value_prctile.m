function [p] = value_prctile(historicalData,exogenousVariable)
% if testing multiple values, permutations should be rows


nless = sum(historicalData < exogenousVariable);
nequal = sum(historicalData == exogenousVariable);
p = (nless + 0.5*nequal) / length(historicalData);