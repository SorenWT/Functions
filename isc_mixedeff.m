function stats = isc_mixedeff(datamats,type,covars,mdl)

currdir = pwd;

templatedir = which('isc_mixedeff.m');
templatedir = erase(templatedir,'isc_mixedeff.m');
templatedir = fullfile(templatedir,'templates','afni_demo');
cd(templatedir)

func = 'ARzs_CW_avvr.DEL+orig';

[err, Vfunc, Infofunc, ErrMessage] = BrikLoad (func);

dimns = [4 4 3];
maskcoord = {2 2 2};
%Vfunc = Vfunc(:,:,:,1);
%mask = zeros(size(Vfunc));
mask = zeros(dimns(1),dimns(2),dimns(3));
mask(maskcoord{1}-1:maskcoord{1}+1,maskcoord{2}-1:maskcoord{2}+1,maskcoord{3}-1:maskcoord{3}+1) = 1; % 3x3 voxel mask
%mask = logical(mask);

InfoDelOut = Infofunc;
InfoDelOut.RootName = ''; %that'll get set by WriteBrik
InfoDelOut.DATASET_RANK(1) = 1; %one sub-brick
InfoDelOut.DATASET_DIMENSIONS = dimns;
InfoDelOut.BRICK_TYPES = [1]; %store data as shorts
InfoDelOut.BRICK_STATS = []; %automatically set
InfoDelOut.BRICK_FLOAT_FACS = [];%automatically set
InfoDelOut.BRICK_LABS = 'Mask';
InfoDelOut.IDCODE_STRING = '';%automatically set

OptDelOut.Scale = 1;
OptDelOut.Prefix = 'Mask_3x3';
OptDelOut.verbose = 0;
[err, ErrMessage, InfoDelOut] = WriteBrik (mask, InfoDelOut, OptDelOut);

alldata = cat(2,datamats{:});
alldesign = Make_designVect(cellfun(@(d)size(d,2),datamats,'UniformOutput',true));

if ~isempty(covars)
    if length(covars) > 1
        for i = 1:2
            varnames{i} = covars{i}.Properties.VariableNames;
        end
        commonvars = intersect(varnames{1},varnames{2});
        for i = 1:2
            covars{i} = covars{i}(:,contains(covars{i}.Properties.VariableNames,commonvars));
        end
    end
end

allcovars = cat(1,covars{:}); % assume covariates are in table format
n = fieldnames(allcovars);
% center quantitative variables
for i = 1:size(allcovars,2)
    if length(unique(allcovars.(n{i}))) > 3
        allcovars.(n{i}) = allcovars.(n{i})-nanmean(allcovars.(n{i}));
    end
end
covars = [];

mdlvars = tokenize(mdl,'+');
mdlcovars = mdlvars;
mdlcovars(find(strcmpi(mdlcovars,'grp'))) = [];
mdlcovars(find(contains(mdlcovars,'Subj'))) = [];
mdlcovars(find(strcmpi(mdlcovars,'1'))) = [];

covars = struct;
covars.name = cell(0);
covars.data = cell(0);
covars.type = cell(0);

for i = 1:length(mdlcovars)
    if contains(mdlcovars{i},'sum')
        varname = char(extractBefore(mdlcovars{i},'sum'));
        func = @plus;
    elseif contains(mdlcovars{i},'dif')
        varname = char(extractBefore(mdlcovars{i},'dif'));
        func = @minus;
    end
    tmp = allcovars.(varname);
    covars.data{i} = bsxfun(func,tmp,tmp');
    if contains(mdlcovars{i},'dif')
        covars.data{i} = abs(covars.data{i});
    end
    covars.name{i} = mdlcovars{i};
    if length(unique(tmp)) <= 3
        covars.type{i} = 'factor';
        if ~isempty(find(covars.data{i}==0))
            covars.data{i} = covars.data{i}+1;
        end
        % can only handle covariate factors with two levels right now
        covars.factnames{i} = {[varname '11'],[varname '12'],[varname '22']};
    else
        covars.type{i} = 'quant';
    end
end


switch type
    case {'spearman','pearson','kendall'}
        allisc = rtoz(corr(alldata,'Type',type));
    case 'eucdist'
        allisc = eucdist(alldata);
        allisc = sqrt(allisc);
        %         alliscvals = boxcox(belowDiag(allisc));
        %         allisc = bdiagtomat(alliscvals,size(allisc,1));
        %         allisc = allisc+allisc';
        
end

InfoDelOut.BRICK_LABS = 'Correlation';
tblcommand = [];
tblcommand.Subj1 = zeros(0); tblcommand.Subj2 = zeros(0);
if contains(mdl,'grp')
    tblcommand.grp = zeros(0); 
end
tblcommand.InputFile = zeros(0);



% add more columns for demographics
count = 1;

for i = 1:size(allisc,1)
    for ii = 1:(i-1)
        if i ~= ii
            %Vfunc(maskcoord(1),maskcoord(2),maskcoord(3)) = allisc(i,ii);
            OptDelOut.Prefix = ['s' num2str(i) '_' num2str(ii) '_corr'];
            if i == 2 && ii == 1
                WriteBrik (repmat(allisc(i,ii),dimns(1),dimns(2),dimns(3)), InfoDelOut, OptDelOut);
            end
            grp1 = alldesign(i); grp2 = alldesign(ii);
            isc(count) = allisc(i,ii);
            
            %only two groups possible right now
            if contains(mdl,'grp')
                if grp1 == grp2 && grp1 == 1
                    tblcommand.grp(count) = 1;
                elseif grp1 == grp2 && grp1 == 2
                    tblcommand.grp(count) = -1;
                elseif grp1 ~= grp2
                    tblcommand.grp(count) = 0;
                end
            end
            
            for q = 1:length(covars.name)
                if strcmpi(covars.type{q},'quant')
                    tblcommand.(covars.name{q})(count) = covars.data{q}(i,ii);
                elseif strcmpi(covars.type{q},'factor')
                    tblcommand.(covars.name{q})(count) = {covars.factnames{q}{covars.data{q}(i,ii)}};
                end
            end
            
            tblcommand.Subj1{count} = ['s' num2str(i)];
            tblcommand.Subj2{count} = ['s' num2str(ii)];
            tblcommand.InputFile{count} = ['s' num2str(i) '_' num2str(ii) '_corr+orig'];
            count = count+1;
        end
    end
end
tblcommand = structfun(@vert,tblcommand,'UniformOutput',false);
tblcommand = struct2table(tblcommand);
tmp = tblcommand.Properties.VariableNames;
indx = find(strcmpi(tmp,'InputFile'));
% put inputfile at the end
tblcommand = tblcommand(:,[except(1:size(tblcommand,2),indx) indx]);
writetable(tblcommand,fullfile(pwd,'filestbl.txt'),'Delimiter',' ')
!awk 'NR > 1{print line" \\"}{line=$0;}END{print $0" "}' filestbl.txt > filestbl2.txt

isc = array2table(vert(isc),'VariableNames',{'ISC'});
isc.ISC(find(isinf(isc.ISC))) = NaN;
writetable(isc,fullfile(pwd,'isc.csv'))

tblfacts = tblcommand;
tblfacts.Subj1 = []; tblfacts.Subj2 = []; tblfacts.InputFile = [];
tblvarnames = tblfacts.Properties.VariableNames;
nfacts = size(tblfacts,2);
mdlfacts = [];
for i = 1:size(tblfacts,2)
    if iscell(tblfacts.(tblvarnames{i}))
        nfacts = nfacts+length(unique(tblfacts.(tblvarnames{i})))-2;
        tmp = repmat({tblvarnames{i}},1,length(unique(tblfacts.(tblvarnames{i})))-1);
        mdlfacts = cat(2,mdlfacts,horz(cellcat('dum',tmp,'',0)));
    else
        mdlfacts = cat(2,mdlfacts,{tblvarnames{i}});
    end
end
if contains(mdl,'grp')
    types = [{'quant'} horz(covars.type)];
else
    types = horz(covars.type);
end


gltcode = struct;
gltcode.ave = zeros(nfacts+1,1); gltcode.ave(1) = 1;
for i = 1:size(tblfacts,2)
    if strcmpi(tblvarnames{i},'grp')
        gltcode.g11vg22 = zeros(nfacts+1,1); gltcode.g11vg22(find(contains(mdlfacts,'grp'))+1) = 1;
        gltcode.g11 = zeros(nfacts+1,1); gltcode.g11(1) = 1; gltcode.g11(find(contains(mdlfacts,'grp'))+1) = 0.5;
        gltcode.g22 = zeros(nfacts+1,1); gltcode.g22(1) = 1; gltcode.g22(find(contains(mdlfacts,'grp'))+1) = -0.5;
    else
        if strcmpi(types{i},'quant')
            gltcode.(tblvarnames{i}) = zeros(nfacts+1,1); gltcode.(tblvarnames{i})(find(contains(mdlfacts,tblvarnames{i}))+1) = 1;
        elseif strcmpi(types{i},'factor')
            vn = tblvarnames{i}; dumindxs = find(contains(mdlfacts,vn))+1;
            comps = {[vn '11'],[vn '12'],[vn '22'],[vn '11v' vn '22'],[vn '11v' vn '12'],[vn '12v' vn '22']};
            for ii = 1:length(comps)
                gltcode.(comps{ii}) = zeros(nfacts+1,1);
            end
            gltcode.(comps{1})(1) = 1; gltcode.(comps{1})(dumindxs(1)) = 1;
            gltcode.(comps{2})(1) = 1; gltcode.(comps{2})(dumindxs(2)) = 1;
            gltcode.(comps{3})(1) = 1; gltcode.(comps{3})(dumindxs(1)) = -1; gltcode.(comps{3})(dumindxs(2)) = -1;
            gltcode.(comps{4})(dumindxs(1)) = 2; gltcode.(comps{4})(dumindxs(2)) = 1;
            gltcode.(comps{5})(dumindxs(1)) = 1; gltcode.(comps{5})(dumindxs(2)) = -2;
            gltcode.(comps{6})(dumindxs(1)) = 1; gltcode.(comps{6})(dumindxs(2)) = 2;
            gltcode.ave(dumindxs(2)) = -0.5;
        end
    end
end

gltfields = fieldnames(gltcode);
gltstr = [];
for i = 1:length(gltfields)
    gltnumstr = cellstr(num2str(vert(gltcode.(gltfields{i}))));
    gltnumstr = erase(gltnumstr,' ');
    gltnumstr = cellcat(' ',gltnumstr,'',1);
    gltnumstr = cat(2,gltnumstr{:}); gltnumstr(end) = [];
    gltstr = cat(2,gltstr,['-gltCode ' gltfields{i} ' ' char(39) gltnumstr char(39) ' \' newline]);
end
qvars = find(strcmpi(types,'quant'));
qvarsstr = tblvarnames(qvars);
qvarsstr = cellcat(',',qvarsstr,'',1);
qvarsstr = cat(2,qvarsstr{:}); qvarsstr(end) = [];
qvarsstr = [char(39) qvarsstr char(39)];

% command = ['3dISC -prefix ISCout -jobs 2 \' newline ...
%     '-mask Mask_3x3+orig \' newline ...
%     '-qVars grp \' newline ...
%     '-model ''grp+(1|Subj1)+(1|Subj2)'' \' newline ...
%     '-gltCode ave ''1 0'' \' newline ...
%     '-gltCode G11vG22 ''0 1'' \' newline ...
%     '-gltCode G11 ''1 0.5'' \' newline ...
%     '-gltCode G22 ''1 -0.5'' \' newline ...
%     '-dataTable @' fullfile(pwd,'filestbl2.txt') ' \'];

command = ['export R_MAX_VSIZE=8GB' newline ...
    '3dISC -prefix ISCout -jobs 2 \' newline ...
    '-mask Mask_3x3+orig \' newline ...
    '-qVars ' qvarsstr ' \' newline ...
    '-model ' char(39) mdl char(39) ' \' newline ...
    gltstr ...
    '-dataTable @' fullfile(pwd,'filestbl2.txt') ' \'];

setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
path = getenv('PATH');
setenv('PATH',[path ':/Users/Soren/abinexport PATH=$PATH:/Users/Soren/abin']);
system(['afni' newline newline command])



[~,iscin,iscinfo] = BrikLoad('ISCout+orig');

try
    dimord = tokenize(iscinfo.BRICK_LABS,'~');
    
    stats = struct;
    for i = 1:length(gltfields)
        stats.estimate(i) = iscin(maskcoord{:},find(strcmpi(dimord,gltfields{i})));
        stats.t(i) = iscin(maskcoord{:},find(strcmpi(dimord,[gltfields{i} ' t'])));
        if stats.t(i) > 0
            stats.p(i) = 2*tcdf(stats.t(i),size(allisc,1)-1,'upper');
        else
            stats.p(i) = 2*tcdf(stats.t(i),size(allisc,1)-1);
        end
    end
    
    stats = structfun(@vert,stats,'UniformOutput',false);
    stats = struct2table(stats);
    stats.Properties.RowNames = gltfields;
catch
    system('for f in s*_corr+orig*; do rm "$f"; done')
    !rm Mask_3x3+orig*
    !rm ISCout+orig*
    !rm filestbl*
    !rm isc.csv
    error(lasterror)
end

%!rm s*_corr+orig*
system('for f in s*_corr+orig*; do rm "$f"; done')
!rm Mask_3x3+orig*
!rm ISCout+orig*
!rm filestbl*
!rm isc.csv

cd(currdir)

