function res = rsa(data,nperm,type,varargin)

if CheckInput(varargin,'mask')
    mask = EasyParse(varargin,'mask');
else
    mask = ones(size(corr(data{1})));
end

if ~exist('type','var')
    type = repmat({'pearson'},1,length(data));
end

if ischar(type)
    type = repmat({type},1,length(data));
end

for i = 1:length(data)
    switch type{i}
        case {'spearman','pearson','kendall'}
            corrmat{i} = corr(data{i},'Type',type{i},'rows','pairwise');
            corrmat{i} = corrmat{i}.*mask;
            corrvals{i} = vert(rtoz(belowDiag(corrmat{i})));
        case 'eucdist'
            corrmat{i} = eucdist(data{i});
            corrmat{i} = corrmat{i}.*mask;
            corrvals{i} = vert(belowDiag(corrmat{i}));
            corrvals{i} = log(corrvals{i});
        otherwise
            corrmat{i} = squareform(pdist(data{i}',type{i}));
            corrmat{i} = corrmat{i}.*mask;
            corrvals{i} = vert(belowDiag(corrmat{i}'));
    end
    
    if any(any(isnan(data{i})))
        warning(['Data matrix ' num2str(i) ' contains NaNs - these observations will be ignored. Ignore these at your own risk!'])
    end
end

orig_stat = corr(corrvals{1},corrvals{2},'rows','pairwise');


%alldata = cat(2,data{:});
origdesign = Make_designVect(cellfun(@(dat)size(dat,2),data,'UniformOutput',true));


for i = 1:nperm
    %design = origdesign(randperm(length(origdesign)));
    for ii = 1:length(data)
        switch type{ii}
            case {'spearman','pearson','kendall'}
                newcorrmat{ii} = corr(data{ii}(:,randperm(size(data{ii},2))),'Type',type{ii},'rows','pairwise');
                newcorrmat{ii} = newcorrmat{ii}.*mask;
                newcorrvals{ii} = rtoz(belowDiag(newcorrmat{ii}));
            case 'eucdist'
                newcorrmat{ii} = eucdist(data{ii}(:,randperm(size(data{ii},2))));
                newcorrmat{ii} = newcorrmat{ii}.*mask;
                newcorrvals{ii} = belowDiag(newcorrmat{ii});
                %newcorrvals{ii} = log(newcorrvals{ii});
            otherwise
                newcorrmat{ii} = squareform(pdist(data{ii}(:,randperm(size(data{ii},2)))',type{ii}));
                newcorrmat{ii} = newcorrmat{ii}.*mask;
                newcorrvals{ii} = vert(belowDiag(newcorrmat{ii}));
        end
    end
    
    perm_stat(i) = corr(newcorrvals{1},newcorrvals{2},'rows','pairwise');
end

p = value_prctile(perm_stat,orig_stat);

if p > 0.5
    p = 1-p;
end

p = 2*p; %two-sided test

res.pperm = p; res.statobs = orig_stat; res.type = type;
res.permstat = perm_stat; res.corrmats = corrmat; res.corrvals = corrvals;

% jackknifing to get pseudo-leverage
if CheckInput(varargin,'jackknife') && EasyParse(varargin,'jackknife','on')
    for ii = 1:2
        for i = 1:size(data{ii},1)
            switch type{ii}
                case {'spearman','pearson','kendall'}
                    newcorrmat{ii} = corr(data{ii}(except(1:size(data{ii},1),i),:),'Type',type{ii},'rows','pairwise');
                    newcorrmat{ii} = newcorrmat{ii}.*mask;
                    newcorrvals{ii} = rtoz(belowDiag(newcorrmat{ii}));
                case 'eucdist'
                    newcorrmat{ii} = eucdist(data{ii}(except(1:size(data{ii},1),i),:));
                    newcorrmat{ii} = newcorrmat{ii}.*mask;
                    newcorrvals{ii} = belowDiag(newcorrmat{ii});
                    %newcorrvals{ii} = log(newcorrvals{ii});
                otherwise
                    newcorrmat{ii} = squareform(pdist(data{ii}(except(1:size(data{ii},1),i),:)',type{ii}));
                    newcorrmat{ii} = newcorrmat{ii}.*mask;
                    newcorrvals{ii} = vert(belowDiag(newcorrmat{ii}));
            end
            newcorrvals{except(1:2,ii)} = corrvals{except(1:2,ii)};
            
            pseudolev{ii}(i) = corr(newcorrvals{1},newcorrvals{2},'rows','pairwise');
        end
        pseudolev{ii} = orig_stat-pseudolev{ii};
    end
    
    res.pseudolev = pseudolev; % positive values of pseudolev indicate positive contribution to RSA
end



