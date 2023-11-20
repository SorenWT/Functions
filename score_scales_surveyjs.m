function [scales]=score_scales_surveyjs(jatosdata,scaleinfo)
% the order should match the order that scaleinfo is in (since one was made
% from the other)

sclnames = getfield_list(scaleinfo,'shortname');
if ~isempty(jatosdata)
    fnames = fieldnames(jatosdata);
end
for i = 1:length(scaleinfo)
    if ~isempty(jatosdata)
        jatosdata.(fnames{i}) = structfun(@(d)str2num(erase(d,'Column ')),jatosdata.(fnames{i}),'UniformOutput',false);
        jatosdata.(fnames{i}) = struct2table(jatosdata.(fnames{i}));
        %tmptbl = cell2table(scaleinfo{i}.respcoding(jatosdata.(fnames{i}){:,:}),'VariableNames',scaleinfo{i}.items);
        %scl = calc_scales(tmptbl,scaleinfo{i});
        
        rawscale = jatosdata.(fnames{i});
        
        tmprows = rawscale.Properties.VariableNames;
        tmprows = cellfun(@(d)str2num(erase(d,'Row')),tmprows,'UniformOutput',true);
        
        rawscale.Properties.VariableNames = scaleinfo{i}.items(tmprows);
        if length(scaleinfo{i}.items)~=size(rawscale,2)
            fprintf(['Mismatched number of items on scale ' sclnames{i} ': expected ' ...
                num2str(length(scaleinfo{i}.items)) ' items but only found ' num2str(size(rawscale,2))])
        end
        
        tmpnans = array2table(NaN(1,length(scaleinfo{i}.items)),'VariableNames',scaleinfo{i}.items);
        tmpnans{:,tmprows} = rawscale{:,:};
        rawscale = tmpnans;
        
        if any(scaleinfo{i}.reverse==-1)
            rawscale{1,find(scaleinfo{i}.reverse==-1)} = length(scaleinfo{i}.respcoding)+1-rawscale{1,find(scaleinfo{i}.reverse==-1)};
        end
    else
        rawscale = array2table(NaN(1,length(scaleinfo{i}.items)),'VariableNames',scaleinfo{i}.items);
    end
    
    
    
    info = scaleinfo{i};
    scales.raw.(sclnames{i}) = rawscale;
    factnames = fieldnames(scaleinfo{i}.factors);
    
    
    for ii = 1:length(factnames)
        % use mean and not sum so that people who missed a response are
        % comparable to the full-scale people
        
        indx = match_str(rawscale.Properties.VariableNames,info.items(info.factors.(factnames{ii})));
        scales.(sclnames{i}).(factnames{ii}) = nanmean(rawscale{1,indx})*length(info.factors.(factnames{ii}));
        scales.allscales.([sclnames{i} '_' factnames{ii}]) = nanmean(rawscale{1,indx})*length(info.factors.(factnames{ii}));
    end
    
    
end

scales.allscales = struct2table(scales.allscales);


