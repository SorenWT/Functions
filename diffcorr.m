function [z, p] = diffcorr(r1, r2, n1, n2)

% Modified from Voytek et al.'s erpac software 

r1 = vert(r1);
r2 = vert(r2);
n1 = vert(n1); 
n2 = vert(n2);

if length(r1) ~= length(n1)
   n1 = repmat(n1,length(r1),1); 
end

if length(r2) ~= length(n2)
   n2 = repmat(n2,length(r2),1); 
end


% Fisher's z-transform to normalize correlation coefficients
z1 = 0.5 * log((1 + r1)./(1 - r1));
z2 = 0.5 * log((1 + r2)./(1 - r2));
zdiff = z1 - z2;

% Calculate sigma given the number of trials
se = sqrt((1./(n1 - 3)) + (1./(n2 - 3)));

% Get z-score of the differences between correlation coefficients
z = zdiff./se;

% Calculate p-value
p = 1 - normcdf(abs(z));

p = 2*p; %two-tailed test

clear r1 r2 n1 n2 z1 z2 zdiff se

