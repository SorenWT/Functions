function [ace,smry,compare,reduc,reducsmry] = heritability_ACE(mzdata,dzdata)
% mzdata and dzdata should be tables
% frmla is the relevant anova formula

currdir = pwd;

writetable(mzdata,fullfile(currdir,'mztbl.csv'))
writetable(dzdata,fullfile(currdir,'dztbl.csv'))

path = which('heritability_ACE');

path = erase(path,'heritability_ACE.m');

funcname = 'heritability_ACE.R';
mzin = 'mztbl.csv';
dzin = 'dztbl.csv';
fileout = 'output';

[~,sysname] = system('hostname');

if contains(sysname,'Mac-the-Knife')
    setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
end

%system('ulimit -s 16384')

system(['R -e ''source("' fullfile(path,funcname) '"); heritability_ACE("'...
    fullfile(currdir,mzin) '","' fullfile(currdir,dzin) '","' ... 
    fullfile(currdir,fileout) '")'''])

ace = jsonread(fullfile(currdir,[fileout '_acesmry.json']));
smry = jsonread(fullfile(currdir,[fileout '_fullsmry.json']));
compare = jsonread(fullfile(currdir,[fileout '_mdlcompare.json']));
reduc = jsonread(fullfile(currdir,[fileout '_reducsmry.json']));
reducsmry = jsonread(fullfile(currdir,[fileout '_fullsmry_reduc.json']));


ace.ci = struct2table(ace.ci);
ace.ci.Properties.RowNames = {'Intercept','a','c','e'};
ace.ci.Properties.VariableNames = {'l_ci','u_ci'};


reduc.ci = struct2table(reduc.ci);
%reduc.ci.Properties.RowNames = {'Intercept','a','c','e'};
reduc.ci.Properties.VariableNames = {'l_ci','u_ci'};

%compare = compare.NULL;
compare = struct2table(compare);

system(['rm ' mzin]);
system(['rm ' dzin]);
system(['rm ' fileout '_acesmry.json']);
system(['rm ' fileout '_fullsmry.json']);





