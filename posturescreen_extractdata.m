function [posdat_all,posdat_avg,alltbl,avgtbl] = posturescreen_extractdata(report)

reptext = extractFileText(report);
reptext = char(reptext);

reptext = tokenize(reptext,newline);
empties = cellfun(@isempty,reptext,'UniformOutput',true);
reptext(empties) = [];

indx1 = find(strcmpi(reptext,'Anterior View Right View Posture Displacements'));
tbl1 = array2table(zeros(5,4),'VariableNames',{reptext{indx1+1}(13:end) [reptext{indx1+2} reptext{indx1+3}] ...
    [reptext{indx1+4} reptext{indx1+5}] [reptext{indx1+6} reptext{indx1+7}]},'RowNames',...
    reptext(indx1+8:indx1+12));

% get the first table
dat = zeros(20,1);
count = 1;
for i = 1:20
   tmp = reptext{indx1+12+count};
   if strcmpi(tmp,'0')
       dat(i) = 0;
       count = count+1;
   elseif strcmpi(tmp,'N/A')
       dat(i) = NaN;
       count = count+1;
   else
       atLeastOneDigit = '(?=.*?\d).*';
       count2= 1;
       while 1
           tmp2 = reptext{indx1+12+count+count2};
           if isempty(regexp(tmp2,atLeastOneDigit,'once')) && ~strcmpi(tmp2,'N/A')
               tmp = [tmp tmp2];
               count2 = count2+1;
           else
               break
           end
       end
       dat(i) = posturestring_decode(tmp);
       count = count+count2;
   end
end

tbl1{:,:} = reshape(dat,5,4);

indx2 = find(strcmpi(reptext,'Posture Displacements'));
tbl2 = array2table(zeros(10,4),'VariableNames',{reptext{indx2+1}(13:end) [reptext{indx2+2} reptext{indx2+3}] ...
    [reptext{indx2+4} reptext{indx2+5}] [reptext{indx2+6} reptext{indx2+7}]},'RowNames',...
    reptext(indx2+8:indx2+17));

dat = zeros(40,1);
count = 1;
for i = 1:40
   tmp = reptext{indx2+17+count};
   if strcmpi(tmp,'0')
       dat(i) = 0;
       count = count+1;
   elseif strcmpi(tmp,'N/A')
       dat(i) = NaN;
       count = count+1;
   else
       atLeastOneDigit = '(?=.*?\d).*';
       count2= 1;
       while 1
           tmp2 = reptext{indx2+17+count+count2};
           if isempty(regexp(tmp2,atLeastOneDigit,'once')) && ~strcmpi(tmp2,'N/A')
               tmp = [tmp tmp2];
               count2 = count2+1;
           else
               break
           end
       end
       dat(i) = posturestring_decode(tmp);
       count = count+count2;
   end
end

tbl2{:,:} = reshape(dat,10,4);

tbl1 = cat(1,tbl1,array2table(NaN(5,4),'VariableNames',tbl1.Properties.VariableNames,...
    'RowNames',tbl2.Properties.RowNames(6:end)));

tbl2.Properties.VariableNames{3} = ['post_' tbl2.Properties.VariableNames{3}];
tbl2.Properties.VariableNames{4} = ['post_' tbl2.Properties.VariableNames{4}];


alltbl = cat(2,tbl1,tbl2);
alltbl.Properties.RowNames = erase(alltbl.Properties.RowNames,' ');

posdat_all = flatten_table(alltbl);
posdat_all(:,isnan(posdat_all{:,:})) = [];

if contains(report,'Anterior')
   badvars = find(contains(posdat_all.Properties.VariableNames,'Anterior'));
   posdat_all{:,badvars} = NaN(size(posdat_all{:,badvars}));
end

if contains(report,'Posterior')
   badvars = find(contains(posdat_all.Properties.VariableNames,'Posterior'));
   posdat_all{:,badvars} = NaN(size(posdat_all{:,badvars}));
end

if contains(report,'Lateral')
   badvars = find(contains(posdat_all.Properties.VariableNames,'Lateral'));
   posdat_all{:,badvars} = NaN(size(posdat_all{:,badvars}));
end

if contains(report,'post_Lateral')
   badvars = find(contains(posdat_all.Properties.VariableNames,'post_Lateral'));
   posdat_all{:,badvars} = NaN(size(posdat_all{:,badvars}));
end

avgtbl = alltbl;
avgtbl(6:end,:) = [];

avgtbl.LateralTranslations = nanmean([avgtbl.LateralTranslations,avgtbl.post_LateralTranslations],2);
avgtbl.LateralAngulations = nanmean([avgtbl.LateralAngulations,avgtbl.post_LateralAngulations],2);
avgtbl = avgtbl([1 2 4 5],3:4);

posdat_avg = flatten_table(avgtbl);




end

function num = posturestring_decode(str)

keywords = {'left','right','anterior','posterior','flexed','extended'};
signs = [1 -1 1 -1 1 -1];

if contains(str,'º')
    num = str2num(char(extractBefore(str,'º')));
else
    num = str2num(char(extractBefore(str,' ')));
end

num = num.*signs(find(cellfun(@(k)contains(str,k),keywords,'UniformOutput',true)));
end

function tblout = flatten_table(tblin)

tblout = reshape(tblin{:,:},[],1);
varnamesout = cell(size(tblout));

for i = 1:length(varnamesout)
   [i1,i2] = ind2sub(size(tblin),i);
   varnamesout{i} = [tblin.Properties.RowNames{i1} '_' tblin.Properties.VariableNames{i2}];
end

tblout = array2table(tblout','VariableNames',varnamesout);

end