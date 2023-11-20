function ICC = ICCbin_anova(data)

% Formula from Wu et al., 2012: Comparison of Methods for Estimating the Intraclass Correlation Coefficient for Binary Responses in Cancer Prevention Cluster Randomized Trials

data(any(isnan(data),2),:) = [];

k = size(data,1);
ni = repmat(size(data,2),size(data,1),1); 

Z = sum(data,2); 
Ztotsq = sum(Z).^2;

N = numel(data);

na = (1/(k-1))*(N-sum(ni.^2)/N);

MSB  = (1/(k-1))*(sum((Z.^2)./ni) - Ztotsq/N);
MSW = (1/(N-k))*(sum(Z)-sum((Z.^2)./ni));

ICC = (MSB-MSW)./(MSB+(na-1)*MSW);