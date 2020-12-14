function [ H, p ] = lemeshow( pred, target, D )

if nargin<3
    D=10; % Number of deciles to use.
    % Note that the HL test is only defined for D=10 (or D=8 if applying to a training set)
end
if min(size(pred))~=1 || min(size(target))~=1
    error('Function only accepts vector inputs.');
end

% ensure they are column vectors
target = target(:); pred = pred(:);
N = numel(target);

[pred,idxSort] = sort(pred,1,'ascend');
target = target(idxSort);

% create an index of the outcomes into each of the D deciles
idxDecile = ceil( (1:N)/N*D )';

Htbl = zeros(D,2);
Htbl(:,1) = accumarray(idxDecile,target);
Htbl(:,2) = accumarray(idxDecile,pred);

% calculate expected probability in each decile
n = hist(idxDecile,1:D);
for i = 1:size(Htbl,1)
   H(i) = (((Htbl(i,1)-Htbl(i,2)).^2)./Htbl(i,2))+((((n(i)-Htbl(i,1))-(n(i)-Htbl(i,2))).^2)./(n(i)-Htbl(i,2)));
end

H = sum(H);

p = 1-chi2cdf(H,D-2);

%p = H(:,2)./n;
%HL = (H(:,1)-H(:,2)).^2 ./ (p.*n.*(1-p)+eps); % eps prevents division by 0
%HL = sum(HL);


end