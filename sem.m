function [summary,paramscov] = sem(data,syntax,clustervar)

if nargin < 3
    clustervar = 'none'; 
end

currdir = pwd;

% create all two-way interaction terms, just in case



useddata = zeros(1,width(data));
for i = 1:length(data.Properties.VariableNames)
    if any(contains(syntax,data.Properties.VariableNames{i}))
       useddata(i) = 1;
    end
end

numericdata = varfun(@isnumeric,data);
numericdata = numericdata{:,:};

useddata = useddata & numericdata; % get rid of any non-numerics that snuck in there because of being subsets of names

%useddata = find(useddata);

names = data.Properties.VariableNames(useddata);

design = x2fx(data{:,useddata},'interaction');
design(:,1) = [];
for i = 1:length(names)
    for ii = 1:length(names)
    rxnnames{ii,i} = [names{i} '.' names{ii}];
    end
end
allnames = [names reshape(rxnnames(find(tril(ones(length(names)))-eye(length(names)))),1,[])];
allnames(cellfun(@length,allnames)>63) = cellfun(@(d)indexme(d,1:63),allnames(cellfun(@length,allnames)>63),'UniformOutput',false);
design = array2table(design,'VariableNames',allnames);
design{:,(length(names)+1):end} = nancenter(design{:,(length(names)+1):end},1);
design = cat(2,design,data(:,~useddata)); % keep non-numeric data (e.g. ids) in but don't do interactions

writetable(design,fullfile(currdir,'datatbl.csv'))

path = which('sem');

path = erase(path,'sem.m');

funcname = 'sem_matlab.R';
filein = 'datatbl.csv';
fileout = 'output';

[~,sysname] = system('hostname');

if contains(sysname,'Mac-the-Knife')
    setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
end

%syntax = replace(syntax,':','.');

%system('ulimit -s 16384')

system(['R -e ''source("' fullfile(path,funcname) '"); sem_matlab("'...
    fullfile(currdir,filein) '","' fullfile(currdir,fileout) '","' syntax '","' clustervar '")'''])

%mdl = jsonread(fullfile(currdir,[fileout '_mdl.json']));
summary = jsonread(fullfile(currdir,[fileout '_summary.json']));
%paramscov = jsonread(fullfile(currdir,[fileout '_cov.json']));

%system(['rm ' filein]);
system(['rm ' fileout '_mdl.json']);
system(['rm ' fileout '_summary.json']);
%system(['rm ' fileout '_cov.json']);



summary.tbl = table;
summary.tbl.frmla = strcat(summary.pe.lhs,summary.pe.op,summary.pe.rhs);
summary.tbl.est = summary.pe.est; summary.tbl.se = summary.pe.se; 
summary.tbl.z = summary.pe.z; summary.tbl.pvalue = summary.pe.pvalue;
for i = 1:height(summary.tbl)
    summary.tbl.ci{i} = [summary.tbl.est(i)-1.96*summary.tbl.se(i) summary.tbl.est(i)+1.96*summary.tbl.se(i)];
end








