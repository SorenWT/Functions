function [stats_mat,sigstar_pairs,sigstar_p]=signrank_multcompare(datin,corrmethod)

if nargin < 2
   corrmethod = 'bonf_holm'; 
end

stats_mat = ones(size(datin,2));
for i = 1:size(datin,2)
    for ii = 1:i
        stats_mat(i,ii) = signrank(datin(:,i),datin(:,ii));
    end
end

switch corrmethod
    case 'bonf_holm'
        stats_mat = bdiagtomat(bonf_holm(belowDiag(stats_mat),0.05),size(datin,2));
    case 'fdr'
        stats_mat = bdiagtomat(fdr(belowDiag(stats_mat)),size(datin,2));
end

pairs = find(stats_mat<0.05);
[indx(:,1),indx(:,2)] = ind2sub(pairs,size(stats_mat));

sigstar_pairs = mat2cell(indx,2,ones(size(indx,1),1));
sigstar_p = stats_mat(find(stats_mat<0.05));
