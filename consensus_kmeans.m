function [cl] = consensus_clust(clfunc,datin,k,nreps,varargin)

% func can be any clustering function that requires you to specify a number
% of clusters

argsin = varargin;

for i = 1:nreps
    cl(:,i) = clfunc(datin,k,argsin{:});
    adj(:,:,i) = cl(:,i)==cl(:,i)';
end
meanadj = mean(adj,3);
cl = clfunc(meanadj,k,argsin{:});

