function [stats_mat,sigstar_pairs,sigstar_p]=multcompare_swt(datin,statmethod,corrmethod,xvals)

if nargin < 2
    statmethod = 'signrank';
end

if nargin < 3
    corrmethod = 'bonf_holm';
end

if nargin < 4
    xvals = 1:size(datin,2);
end

stats_mat = ones(size(datin,2));
for i = 1:size(datin,2)
    for ii = 1:i
        switch statmethod
            case 'signrank'
                stats_mat(i,ii) = signrank(datin(:,i),datin(:,ii));
            case 'ttest'
                [~,stats_mat(i,ii)] = ttest(datin(:,i)-datin(:,ii));
        end
    end
end

switch corrmethod
    case 'bonf_holm'
        stats_mat = bdiagtomat(bonf_holm(belowDiag(stats_mat),0.05),size(datin,2));
    case 'fdr'
        stats_mat = bdiagtomat(fdr(belowDiag(stats_mat)),size(datin,2));
end

pairs = find(stats_mat<0.1);
if ~isempty(pairs)
    [indx(:,1),indx(:,2)] = ind2sub(size(stats_mat),pairs);
    
    indx(:,1) = xvals(indx(:,1)); indx(:,2) = xvals(indx(:,2));
    
    sigstar_pairs = mat2cell(indx,ones(size(indx,1),1),2);
    sigstar_p = stats_mat(find(stats_mat<0.1));
else
    sigstar_pairs = {}; sigstar_p = [];
end
