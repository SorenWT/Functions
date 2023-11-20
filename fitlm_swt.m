function lm = fitlm_swt(tbl,modelspec,varargin)

if ischar(modelspec) && istable(tbl)
    % check all the variables and make sure they're centered
    for i = 1:size(tbl,1)
        if isnumeric(tbl{:,i})
            tbl{:,i} = nancenter(tbl{:,i},1);
        elseif iscell(tbl{:,i})
            tbl{:,i} = categorical(tbl{:,i});
            dummy = dummyvar(tbl{:,i}); 
            dummy = nancenter(dummy,1);
            vname = tbl.Properties.VariableNames{i};
            for q = 1:size(dummy,2)
                tbl.([vname '_' num2str(q)]) = dummy(:,q);
            end
            tmp = cellcat('+',vname,'',1); tmp = cat(1,tmp{:});
            modelspec = replace(modelspec,vname,tmp(1:end-1)); % doesn't work atm for interactions with a factor with more than 2 levels
        end
    end
    
    lm = fitlm(tbl,modelspec,varargin{:});
    
    
else
    vnames = {cellcat('Var',cellstr(num2str([1:size(x,2)]')),'',0) 'Depvar'};
    tmptbl = array2table([x y],'VariableNames',vnames);
    
    tmp = cellcat('+',vnames(1:end-1),'',1);
    modelspec = cat(1,tmp{:}); modelspec = modelspec(1:end-1); modelspec = ['Depvar~' modelspec];
    fitlm_swt(tmptbl,modelspec);
end