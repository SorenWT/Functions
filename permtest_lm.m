function [mdl] = permtest_lm(data,frmla)
% data should be a table
% frmla is the relevant anova formula

currdir = pwd;

rmindx = zeros(1,height(data));
for i = 1:height(data)
    dat = table2cell(data(i,:));
    if any(cellfun(@isnan,dat))
       rmindx(i) = 1; 
    end
end
data(find(rmindx),:) = [];

writetable(data,fullfile(currdir,'datatbl.csv'))

path = which('permtest_lm');

path = erase(path,'permtest_lm.m');

funcname = 'permtest_lm.R';
filein = 'datatbl.csv';

[~,sysname] = system('hostname');

if contains(sysname,'Mac-the-Knife')
    setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
end

%system('ulimit -s 16384')

system(['R -e ''source("' fullfile(path,funcname) '"); permtest_lm("'...
    fullfile(currdir,filein) '",' frmla ')'''])

mdl = jsonread(fullfile(currdir,[filein '_results.json']));
%multcompare = jsonread(fullfile(currdir,[fileout '_multcompare.json']));

system(['rm ' filein]);
system(['rm ' filein '_results.json']);





