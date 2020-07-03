function [s, cfg] = ft_statfun_partspearman(cfg, dat, design)

switch cfg.tail
    case 0
        tail = 'both';
        s.critval = [-1+cfg.alpha 1-cfg.alpha];
    case 1
        tail = 'right';
        s.critval = 1-cfg.alpha;
    case -1
        tail = 'left';
        s.critval = -1+cfg.alpha;
end

if iscell(design)
    design = design{1};
    tmp = unique(design(find(~isnan(design))));
    
    for c = 1:size(design,1)
        indices{c,1} = find(design(c,:) == tmp(1));
        indices{c,2} = find(design(c,:) == tmp(2));
    end
    
    p = ones(size(dat,1),1);
    for c = 1:size(dat,1)
        [r,p(c)] = partcorr_pairwise(find(~isnan(dat(c,indices{c,1})))',find(~isnan(dat(c,indices{c,2})))',t3d(cfg.partial),'Type','Spearman','tail',tail);
        if r > 0
            p(c) = 1-p(c);
        else
            p(c) = -1+p(c);
        end
    end
else
    poslabel1        = find(design(cfg.ivar,:)==1);
    poslabel2        = find(design(cfg.ivar,:)==2);
    [dum,i]          = sort(design(cfg.uvar,poslabel1), 'ascend');
    poslabelsperunit(:,1) = poslabel1(i);
    [dum,i]          = sort(design(cfg.uvar,poslabel2), 'ascend');
    poslabelsperunit(:,2) = poslabel2(i);
    
    indices{1} = poslabelsperunit(:,1);
    indices{2} = poslabelsperunit(:,2);
    
    p = ones(size(dat,1),1);
    for c = 1:size(dat,1)
        [r,p(c)] = partialcorr(dat(c,indices{1})',dat(c,indices{2})',t3d(cfg.partial(c,:,:)),'Type','Spearman','tail',tail);
        if r > 0
            p(c) = 1-p(c);
        else
            p(c) = -1+p(c);
        end
    end
end
s.stat = p;

s.df = size(dat,2)-2;
end
