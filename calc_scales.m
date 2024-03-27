function scales = calc_scales(tbl,scaleinfo)

sclnames = fieldnames(scaleinfo);

varnames = tbl.Properties.VariableNames;

if size(tbl,1) > 1
    for i = 1:size(tbl,1)
       scalesin{i} = calc_scales(tbl(i,:),scaleinfo);
       scales = struct;
       fields = fieldnames(scalesin{1});
       for i = 1:length(fields)
           tmp = getfield_list(scalesin,fields{i});
           scales.(fields{i}) = cat(1,tmp{:});
       end
    end
else
    
    for i = 1:length(sclnames)
        info = scaleinfo.(sclnames{i});
        %fprintf(sclnames{i})
        [i1,i2] = match_str(info.items,varnames);
        rawresps = [];
        rawscale = tbl(:,i2);
        for q = 1:width(rawscale)
            tmp = find(strcmpi(info.respcoding,rawscale{1,q}));
            if isempty(tmp)
                tmp = NaN;
            end
            rawresps(q) = tmp;
            if info.reverse(q) == -1
                rawresps(q) = length(info.respcoding)+1-rawresps(q);
            end
        end
        rawscale = array2table(NaN(1,length(info.items)),'VariableNames',info.items);
        rawscale{:,i1} = rawresps;
        %rawscale = array2table(rawresps,'VariableNames',rawscale.Properties.VariableNames);
        %rawscale{1,:} = mat2cell(rawresps,1,ones(length(rawresps),1));
        scales.raw.(sclnames{i}) = rawscale;
        factnames = fieldnames(info.factors);
        if length(info.items)==size(rawscale,2)
            for ii = 1:length(factnames)
                % use mean and not sum so that people who missed a response are
                % comparable to the full-scale people
                scales.(sclnames{i}).(factnames{ii}) = nanmean(rawscale{1,info.factors.(factnames{ii})})*length(info.factors.(factnames{ii}));
                scales.allscales.([sclnames{i} '_' factnames{ii}]) = nanmean(rawscale{1,info.factors.(factnames{ii})})*length(info.factors.(factnames{ii}));
            end
        else
            for ii = 1:length(factnames)
                % use mean and not sum so that people who missed a response are
                % comparable to the full-scale people
                warning(['Mismatched number of items on scale ' sclnames{i} ': expected ' ...
                    num2str(length(info.items)) ' items but only found ' num2str(size(rawscale,2))])
                scales.(sclnames{i}).(factnames{ii}) = nanmean(rawresps(info.factors.(factnames{ii})))*length(info.factors.(factnames{ii}));
                scales.allscales.([sclnames{i} '_' factnames{ii}]) = nanmean(rawresps(info.factors.(factnames{ii})))*length(info.factors.(factnames{ii}));
                %scales.(sclnames{i}).(factnames{ii}) = NaN;
                %scales.allscales.([sclnames{i} '_' factnames{ii}]) = NaN;
            end
        end
        scales.(sclnames{i}) = struct2table(scales.(sclnames{i}));
    end
    scales.allscales = struct2table(scales.allscales);
end