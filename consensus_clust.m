function [cl,ctrs] = consensus_clust(clfunc,datin,k,nreps,varargin)

% func can be any clustering function that requires you to specify a number
% of clusters

argsin = varargin;

for i = 1:nreps
    %if nargout > 1
    %    [cl(:,i),ctrs(:,:,i)] = clfunc(datin,k,argsin{:});
    %else
    cl(:,i) = clfunc(datin,k,argsin{:});
    %end
    adj(:,:,i) = cl(:,i)==cl(:,i)';
end
meanadj = mean(adj,3);
cl = clfunc(meanadj,k,argsin{:});

if CheckInput(argsin,'Distance')
    dist = EasyParse(argsin,'Distance');
else
    dist = 'sqeuclidean';
end

ctrs = getcentroids(cl,datin,dist);

% if nargout > 1
%     if any(all(clout==cl,1))
%        ctrs = ctrs(:,:,find(all(clout==cl,1),1));
%        cl = clout;
%     else
%         while 1
%             [cl,ctrs] = clfunc(datin,k,argsin{:});
%             if all(cl==clout)
%                cl = clout;
%                break
%             end
%         end
%     end
% else
%     cl = clout;
% end

