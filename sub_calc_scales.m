function sub = sub_calc_scales(sub,tbl,scaleinfo)

if ~isfield(sub.raw,'scales')
   sub.raw.scales = table;
end

id = sub.emailid;

sclnames = fieldnames(scaleinfo);

varnames = tbl.Properties.VariableNames;
thisrow = tbl{:,find(contains(varnames,'PleaseEnterYour'),1)};
thisrow = find(strcmpi(thisrow,id),1);

tbl = tbl(thisrow,:);
startqs = find(contains(varnames,'PleaseEnterYour'),1)+1;
sub.raw.scales = [sub.raw.scales tbl(:,startqs:end)];

for i = 1:length(sclnames)
    info = scaleinfo.(sclnames{i});
    [i1,i2] = match_str(info.items,varnames);
    rawscale = tbl(:,i2);
    rawresps = [];
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
    rawscale = array2table(rawresps,'VariableNames',rawscale.Properties.VariableNames);
    %rawscale{1,:} = mat2cell(rawresps,1,ones(length(rawresps),1));
    sub.scales.raw.(sclnames{i}) = rawscale;
    factnames = fieldnames(info.factors);
    if length(info.items)==size(rawscale,2)
    for ii = 1:length(factnames)
        % use mean and not sum so that people who missed a response are
        % comparable to the full-scale people
        sub.scales.(sclnames{i}).(factnames{ii}) = nanmean(rawresps(info.factors.(factnames{ii})))*length(info.factors.(factnames{ii})); 
        sub.scales.allscales.([sclnames{i} '_' factnames{ii}]) = nanmean(rawresps(info.factors.(factnames{ii})))*length(info.factors.(factnames{ii}));
    end
    else
        for ii = 1:length(factnames)
            % use mean and not sum so that people who missed a response are
            % comparable to the full-scale people
            warning(['Mismatched number of items on scale ' sclnames{i} ': expected ' ...
                num2str(length(info.items)) ' items but only found ' num2str(size(rawscale,2))])
            sub.scales.(sclnames{i}).(factnames{ii}) = NaN;
            sub.scales.allscales.([sclnames{i} '_' factnames{ii}]) = NaN;
        end
    end
end

if isstruct(sub.scales.allscales)
    sub.scales.allscales = struct2table(sub.scales.allscales);
end


