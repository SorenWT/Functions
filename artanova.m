function [mdl] = artanova(data,frmla)
% data should be a table
% frmla is the relevant anova formula

currdir = pwd;

writetable(data,fullfile(currdir,'datatbl.csv'))

path = which('artanova');

path = erase(path,'artanova.m');

funcname = 'artanova.R';
filein = 'datatbl.csv';
fileout = 'output';

[~,sysname] = system('hostname');

if contains(sysname,'Mac-the-Knife')
    setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
end

%system('ulimit -s 16384')

system(['R -e ''source("' fullfile(path,funcname) '"); artanova("'...
    fullfile(currdir,filein) '",' frmla ',"' ... 
    fullfile(currdir,fileout) '")'''])

mdl = jsonread(fullfile(currdir,[fileout '_results.json']));
%multcompare = jsonread(fullfile(currdir,[fileout '_multcompare.json']));

system(['rm ' filein]);
system(['rm ' fileout '_results.json']);
system(['rm ' fileout '_multcompare.json']);





