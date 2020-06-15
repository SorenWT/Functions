function [s, cfg] = ft_statfun_fast_signrank(cfg, dat, design)

if iscell(design)
    design = design{1};
    
    tmp = unique(design(find(~isnan(design))));
    for c = 1:size(design,1)
        indices{c,1} = find(design(c,:) == tmp(1));
        indices{c,2} = find(design(c,:) == tmp(2));
    end
    
    teststat = ones(size(dat,1),1);
    for c = 1:size(dat,1)
        teststat(c) = signrank(find(~isnan(dat(c,indices{c,1}))),find(~isnan(dat(c,indices{c,2}))));
%         if median(dat(c,indices{c,1})) > median(dat(c,indices{c,2}))
%             p(c) = 1-p(c);
%         else
%             p(c) = -1+p(c);
%         end
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
    
    teststat = ones(size(dat,1),1);
    for c = 1:size(dat,1)
        teststat(c) = fast_signrank(dat(c,indices{1}),dat(c,indices{2}));
%         if median(dat(c,indices{1})) > median(dat(c,indices{2}))
%             teststat(c) = 1-teststat(c);
%         else
%             teststat(c) = -1+teststat(c);
%         end
    end
end
s.stat = teststat;

% Only works when using approximate method! Use ft_statfun_signrank with
% fewer than 15 samples
switch cfg.tail
%     case 0
%         s.critval = [-0.95 0.95];
%     case 1
%         s.critval = 0.95;
%     case -1
%         s.critval = -0.95;
    case 0
        s.critval = [-1.96 1.96];
    case 1
        s.critval = 1.645;
    case -1
        s.critval = -1.645;
end

s.df = size(dat,2)-2;

end