function [p,orig_stat] = Permtest_ISCD(data,nperm,type,varargin)
% Permtest_ISCD uses a permutation test to test for difference in means
% between two intersubject correlation or intersubject distance matrices
%
% Input arguments:
%
% data are the data from which the ISC matrices are defined.
%    This should be a cell array of matrices in the form
%    observations x subjects
% nperm is the number of permutations to use
%
% type defines what is tested. If testing correlation, input 'pearson',
% 'spearman', or 'kendall'. If testing euclidean distance, input 'eucdist'
%
% this method assumes the same number of observations in each group
%
% Output arguments:
%
% p is the p-value from the permutation test
% orig_stat is the observed value of the test statistic

if ~CheckInput(varargin,'stattype')
    if length(data) == 2
        stattype = 'mediandiff';
    else
        stattype = 'kruskalwallis';
    end
else
    stattype = EasyParse(varargin,'stattype');
end

if ~exist('type','var')
    type = 'pearson';
end

for i = 1:length(data)
    switch type
        case {'spearman','pearson','kendall'}
            corrmat{i} = corr(data{i},'Type',type);
        case 'eucdist'
            corrmat{i} = eucdist(data{i});
    end
    corrvals{i} = vert(rtoz(belowDiag(corrmat{i})));
end

switch stattype
    case 't'
        [~,~,~,stat] = ttest2(corrvals{1},corrvals{2});
        orig_stat = stat.tstat;
    case 'meandiff'
        orig_stat = mean(corrvals{1}) - mean(corrvals{2});
    case 'mediandiff'
        orig_stat = median(corrvals{1}) - median(corrvals{2});
    case 'anova'
        grp = Make_designVect(cellfun(@length,corrvals,'UniformOutput',true));
        [~,stat] = anova1(cat(1,corrvals{:}),grp,'off');
        orig_stat = stat{2,5};
    case 'kruskal'
        grp = Make_designVect(cellfun(@length,corrvals,'UniformOutput',true));
        [~,stat] = kruskalwallis(cat(1,corrvals{:}),grp,'off');
        orig_stat = stat{2,5};
end


alldata = cat(2,data{:});
origdesign = Make_designVect(cellfun(@(dat)size(dat,2),data,'UniformOutput',true));


for i = 1:nperm
    design = origdesign(randperm(length(origdesign)));
    for ii = 1:length(data)
        switch type
            case {'spearman','pearson','kendall'}
                newcorrmat{ii} = corr(alldata(:,design == ii),'Type',type);
            case 'eucdist'
                newcorrmat{ii} = eucdist(alldata(:,design==ii));
        end
        newcorrvals{ii} = rtoz(belowDiag(newcorrmat{ii}));
    end
    
    switch stattype
        case 't'
            [~,~,~,stat] = ttest2(newcorrvals{1},newcorrvals{2});
            perm_stat(i) = stat.tstat;
        case 'meandiff'
            perm_stat(i) = mean(newcorrvals{1}) - mean(newcorrvals{2});
        case 'mediandiff'
            perm_stat(i) = median(newcorrvals{1}) - median(newcorrvals{2});
        case 'anova'
            grp = Make_designVect(cellfun(@length,newcorrvals,'UniformOutput',true));
            [~,stat] = anova1(cat(1,newcorrvals{:}),grp,'off');
            perm_stat(i) = stat{2,5};
        case 'kruskal'
            grp = Make_designVect(cellfun(@length,newcorrvals,'UniformOutput',true));
            [~,stat] = kruskalwallis(cat(1,newcorrvals{:}),grp,'off');
            perm_stat(i) = stat{2,5};
    end
end

p = value_prctile(perm_stat,orig_stat);

if p > 0.5
    p = 1-p;
end

p = 2*p; %two-sided test
