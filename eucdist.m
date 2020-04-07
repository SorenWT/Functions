function dist = eucdist(X,Y)
% should have same syntax as corr (ie. X, Y need to have same number of
% rows)

if nargin < 2
   Y = X; 
end

dist = zeros(size(X,2),size(Y,2));
for i = 1:size(X,2)
    dist(i,:) = vecnorm(repmat(X(:,i),1,size(Y,2))-Y,2,1);
end