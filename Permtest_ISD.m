function [p, orig_tstat] = Permtest_ISD(data1,data2,nperm,varargin)
% Permtest_ISD uses a permutation test to test for difference in means 
% between two intersubject distance matrices
%
% Input arguments: 
%
% data1 and data2 are the data from which the ISC matrices are defined. 
%    They should be in the form observations x subjects
% nperm is the number of permutations to use
%
% this method assumes the same number of observations in each group
%
% Output arguments:
%
% p is the p-value from the permutation test

if ~CheckInput(varargin,'stattype')
   stattype = 'meandiff';
else
    stattype = EasyParse(varargin,'stattype');
end

corrmat1 = EucDist(data1);
%corrmat1 = corr(data1,'Type',type);
corrmat2 = EucDist(data2);
%corrmat2 = corr(data2,'Type',type);
corrvals1 = belowDiag(corrmat1);
corrvals2 = belowDiag(corrmat2);

switch stattype
    case 't'
        [~,~,~,stat] = ttest2(corrvals1,corrvals2);
        orig_stat = stat.tstat;
    case 'meandiff'
        orig_stat = mean(corrvals1) - mean(corrvals2);
end

alldata = cat(2,data1,data2);
origdesign = Make_designVect([size(data1,2) size(data2,2)]);


for i = 1:nperm
    design = origdesign(randperm(length(origdesign)));
    newcorrmat1 = EucDist(alldata(:,find(design == 1)));
    newcorrmat2 = EucDist(alldata(:,find(design == 2)));
    %newcorrmat1 = corr(alldata(:,find(design == 1)),'Type',type);
    %newcorrmat2 = corr(alldata(:,find(design == 2)),'Type',type);
    newcorrvals1 = belowDiag(newcorrmat1);
    newcorrvals2 = belowDiag(newcorrmat2);
    switch stattype
        case 't'
            [~,~,~,stat] = ttest2(newcorrvals1,newcorrvals2);
            perm_stat(i) = stat.tstat;
        case 'meandiff'
            perm_stat(i) = mean(newcorrvals1) - mean(newcorrvals2);
    end
end

p = value_prctile(perm_stat,orig_stat);

if p > 0.5
   p = 1-p; 
end

p = 2*p; %two-sided test

end

function [dist] = EucDist(data)
    for c = 1:size(data,2)
        for cc = 1:size(data,2)
            dist(c,cc) = norm(data(:,c) - data(:,cc));
        end
    end
end