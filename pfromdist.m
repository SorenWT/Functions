function pvalue = pfromdist(statobs,nulldist)
% two-tailed p value from a given null distribution

pvalue = sum(statobs>nulldist)/length(nulldist);

if pvalue > 0.5
    pvalue = 1-pvalue;
end

pvalue = 2*pvalue;

if pvalue > 1
    pvalue = 1;
end