function [p,teststat] = diptest(xmat)

tbl = array2table(xmat,'VariableNames',cellstr(num2str([1:size(xmat,2)]')));

currdir = pwd;

writetable(tbl,fullfile(currdir,'datatbl.csv'))


path = which('diptest');

path = erase(path,'diptest.m');

funcname = 'diptest_matlab.R';
filein = 'datatbl.csv';
fileout = 'output.json';

setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')

system(['R -e ''source("' fullfile(path,funcname) '"); diptest_matlab("'...
    fullfile(currdir,filein) '","' fullfile(currdir,fileout) '")'''])

output = jsonread(fullfile(currdir,fileout),0);

fields = fieldnames(output);

system(['rm ' filein])
system(['rm ' fileout])

for i = 1:length(fields)
    p(i) = output.(fields{i}).p_value;
    teststat(i) = output.(fields{i}).statistic.D;
end
