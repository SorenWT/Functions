function [lm,permlm,restbl] = vartest_control(dat,ctrl,type,names)
% Compares variability of the two samples in dat, while controlling for the
% variables in ctrl
%
% Input arguments: 
%
% dat: the samples for which you wish to compare variability. This should
%   be an ngroups-element cell array, with each cell containing a vector or
%   table (width 1) with one group's data
% ctrl: the variables you want to control for when comparing variability.
%   This should again be an ngroups-element cell array, with each element
%   containing either a nsubjects x nvariables matrix or table. If
%   inputting tables, the variable names must be the same in each table
%
% Optional inputs:
%
% type: defines what type of variability you're comparing. 'CV' compares
%   relative variability using the procedure defined by Schultz (1985).
%   'std' compares absolute variability using the Brown-Forsythe type
%   transformation. 'std_lev' compares absolute variability using the
%   Levene type transformation (default = 'CV')
% names: if inputting matrices and vectors, this is a cell array containing
%   the names of the variables of interest. The order should be independent
%   variable, dependent variable, control 1, control 2, etc. Otherwise
%   default names will be assigned. If this argument is supplied with
%   tables input for dat and ctrl, the values in this argument will be
%   overridden
%
% Outputs:
% 
% lm: the linear model from fitlm corresponding to the variables input
% permlm: if you request more than one output, this will call permtest_lm
%   to test each variable's significance using a permutation test in R. If 
%   this output is not requested, the R function will not be called
% restbl: a table summarizing the coefficient estimates and p-values,
%   including permutation p-values 



if ~exist('type','var') || isempty(type)
   type = 'CV'; 
end

if istable(dat{1})
    depvarname = char(dat{1}.Properties.VariableNames);
    for i = 1:length(dat)
        dat{i} = dat{i}{:,:}; % convert to vector
    end
elseif exist('names','var')
    depvarname = names{2};
else
    depvarname = 'dependent_var';
end

if strcmpi(type,'CV')
    for i = 1:length(dat)
        dat{i} = abs(dat{i}-nanmedian(dat{i}))./nanmedian(dat{i});
    end
elseif strcmpi(type,'std')
    for i = 1:length(dat)
        dat{i} = abs(dat{i}-nanmedian(dat{i}));
    end
elseif strcmpi(type,'std_lev')
    for i = 1:length(dat)
        dat{i} = abs(dat{i}-nanmean(dat{i}));
    end
end

if ~isempty(ctrl)
    if istable(ctrl{1})
        ctrlnames = ctrl{1}.Properties.VariableNames;
    else
        if exist('names','var')
            ctrlnames = names(3:end);
        else
            ctrlnames = cellcat('ctrl',cellstr(num2str([1:size(ctrl,2)]')),'',0);
        end
        for i = 1:length(ctrl)
            ctrl{i} = array2table(ctrl{i},'VariableNames',ctrlnames);
        end
    end
    ctrltbl = cat(1,ctrl{:});
end

if exist('names','var')
   indvarname = names{1}; 
else
    indvarname = 'independent_var';
end

tbl = table;
tbl.(indvarname) = Make_designVect(cellfun(@length,dat,'UniformOutput',true))';
dat = cellfun(@vert,dat,'UniformOutput',false);
tbl.(depvarname) = cat(1,dat{:});
tbl = cat(2,tbl,ctrltbl);

tmp = cellcat('+',ctrlnames,'',1);
tmp = cat(2,tmp{:});
mdl = [depvarname '~' indvarname '+' tmp(1:end-1)];


if length(dat) > 2
    [~,lm] = anovan(tbl.(depvarname),[{tbl.(indvarname)} mat2cell(tbl{:,3:end},height(tbl),ones(1,width(tbl)-2))],'Continuous',[2 4],...
        'varnames',[{indvarname} ctrlnames]); %broken right now
else
    lm = fitlm(tbl,mdl);
end


if nargout > 1
    permlm = permtest_lm(tbl,mdl);
    
    restbl = lm.Coefficients;
    restbl.permutation_p = permlm.table.permutationPr___t__;
end

