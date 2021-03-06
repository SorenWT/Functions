function [p,teststat] = CVtest(data)
% data input can be a N x k array, or a k element cell array
% this wrapper function assumes that the original R file is located in the
% same folder as this wrapper function

if iscell(data)
    vertdata = cellfun(@vert,data,'UniformOutput',false);
    newdat = cat(1,vertdata{:});
    grp = Make_designVect(cellfun(@length,data,'UniformOutput',true))';
else
    newdat = reshape(data,[],1);
    grp = Make_designVect([size(data,1) size(data,1)])';
end

datatbl = array2table([newdat grp],'VariableNames',{'X','G'});

nanindx = find(isnan(newdat));
if ~isempty(nanindx)
    datatbl(nanindx,:) = [];
end

currdir = pwd;

writetable(datatbl,fullfile(currdir,'datatbl.csv'))

path = which('CVtest');

path = erase(path,'CVtest.m');

funcname = 'cvtest_matlab.R';
filein = 'datatbl.csv';
fileout = 'output.json';

setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')

system(['R -e ''source("' fullfile(path,funcname) '"); cvtest_matlab("'...
    fullfile(currdir,filein) '","' fullfile(currdir,fileout) '")'''])

output = jsonread(fullfile(currdir,fileout),0);

p = output.p_value;
teststat = output.D_AD;

system(['rm ' filein])
system(['rm ' fileout])





