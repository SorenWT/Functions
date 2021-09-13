function [p,teststat] = CVtest(data)
% data input can be a N x k array, or a k element cell array
% this wrapper function assumes that the original R file is located in the
% same folder as this wrapper function

if iscell(data)
    vertdata = cellfun(@vert,data,'UniformOutput',false);
    newdat = cat(1,vertdata{:});
    grp = Make_designVect(cellfun(@(d)size(d,1),vertdata,'UniformOutput',true))';
else
    newdat = reshape(data,[],1);
    grp = Make_designVect([size(data,1) size(data,1)])';
end

currdir = pwd;

for i = 1:size(newdat,2)
    datatbl{i} = array2table([newdat(:,i) grp],'VariableNames',{'X','G'});
    nanindx = find(isnan(newdat(:,i)));
    if ~isempty(nanindx)
        datatbl{i}(nanindx,:) = [];
    end
    writetable(datatbl{i},fullfile(currdir,'datatbl.xlsx'),'Sheet',i)
end

path = which('CVtest');

path = erase(path,'CVtest.m');

funcname = 'cvtest_matlab.R';
filein = 'datatbl.xlsx';
fileout = 'output.json';

setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')

system(['R -e ''source("' fullfile(path,funcname) '"); cvtest_matlab("'...
    fullfile(currdir,filein) '","' fullfile(currdir,fileout) '")'''])

output = jsonread(fullfile(currdir,fileout),0);

p = extractfield(output,'p_value');
teststat = extractfield(output,'D_AD');

system(['rm ' filein])
system(['rm ' fileout])





