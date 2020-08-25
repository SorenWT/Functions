function [stats] = ft_iscdiff(cfg,data)

if isempty(cfg)
    cfg = struct;
end

cfg = setdefault(cfg,'cond',{'Condition1' 'Condition2'});
cfg = setdefault(cfg,'thresh',0.1);
cfg = setdefault(cfg,'nrand',1000);
cfg = setdefault(cfg,'corrtype','pearson');

for c = 1:length(data{1}.meas)
    origdesign = Make_designVect([length(data{1}.sub),length(data{2}.sub)]);
    alldata = cat(2,data{1}.data(:,:,c)',data{2}.data(:,:,c)');
    allisc = rtoz(corr(alldata,'type',cfg.corrtype));
    
    isc1{c} = allisc(1:length(data{1}.sub),1:length(data{1}.sub));
    isc2{c} = allisc(length(data{1}.sub)+1:end,length(data{1}.sub)+1:end);
    mixisc{c} = allisc(length(data{1}.sub)+1:end,1:length(data{1}.sub));
    
    medisc1(c) = median(belowDiag(isc1{c})); medisc2(c) = median(belowDiag(isc2{c})); 
    mediscmix(c) = median(reshape(mixisc{c},[],1));
    
    isc1vsisc2obs = medisc1(c)-medisc2(c);
    isc1vsmixobs = medisc1(c) - mediscmix(c);
    isc2vsmixobs = medisc2(c) - mediscmix(c);
    
    for i = 1:cfg.nrand
        assign = randperm(length(origdesign));
        allisc = rtoz(corr(alldata(:,assign),'type',cfg.corrtype));
        isc1perm = allisc(1:length(data{1}.sub),1:length(data{1}.sub));
        isc2perm = allisc(length(data{1}.sub)+1:end,length(data{1}.sub)+1:end);
        mixiscperm = allisc(length(data{1}.sub)+1:end,1:length(data{1}.sub));
        
        isc1vsisc2perm(i) = median(belowDiag(isc1perm))-median(belowDiag(isc2perm));
        isc1vsmixperm(i) = median(belowDiag(isc1perm))-median(reshape(mixiscperm,[],1));
        isc2vsmixperm(i) = median(belowDiag(isc2perm))-median(reshape(mixiscperm,[],1));
    end
    
    p = value_prctile(isc1vsisc2perm,isc1vsisc2obs);
    if p > 0.5
        p = 1-p;
    end
    p = 2*p; %two-sided test
    isc1vsisc2(c) = p;
    
    p = value_prctile(isc1vsmixperm,isc1vsmixobs);
    if p > 0.5
        p = 1-p;
    end
    p = 2*p; %two-sided test
    isc1vsmix(c) = p;
    
        p = value_prctile(isc2vsmixperm,isc2vsmixobs);
    if p > 0.5
        p = 1-p;
    end
    p = 2*p; %two-sided test
    isc2vsmix(c) = p;

   %p(c) = Permtest_ISCD({data{1}.data(:,:,c)',data{2}.data(:,:,c)'},cfg.nrand,'pearson');
%    isc1{c} = corr(data{1}.data(:,:,c)','type','pearson');
%    isc2{c} = corr(data{2}.data(:,:,c)','type','pearson');
%    meanisc1(c) = mean(mean(isc1{c}));
%    meanisc2(c) = mean(mean(isc2{c}));
%    diff(c) = meanisc1(c)-meanisc2(c);
%    if diff(c) > 0
%       dirn{c} = [cfg.cond{1} ' > ' cfg.cond{2}]; 
%    else
%        dirn{c} = [cfg.cond{2} ' > ' cfg.cond{1}];
%    end
end

stats = struct;
stats.meas = data{1}.meas;
stats.([cfg.cond{1} 'vs' cfg.cond{2}]) = isc1vsisc2;
stats.([cfg.cond{1} 'vsmix']) = isc1vsmix;
stats.([cfg.cond{2} 'vsmix']) = isc2vsmix;
stats.(cfg.cond{1}) = isc1;
stats.(cfg.cond{2}) = isc2;
stats.mix = mixisc;
stats.([cfg.cond{1} '_median']) = medisc1;
stats.([cfg.cond{2} '_median']) = medisc2;
stats.mix_median = mediscmix;



% tbl = table;
% p = cat(2,vert(isc1vsisc2),vert(isc1vsmix),vert(isc2vsmix));
% p = p < cfg.thresh;
% p = sum(p,2);
% pmask = find(p);
% tbl.meas = vert(cellfun(@func2str,data{1}.meas(pmask),'UniformOutput',false));
% tbl.([cfg.cond{1} 'vs' cfg.cond{2}]) = isc1vsisc2(pmask);
% tbl.([cfg.cond{1} 'vsmix']) = isc1vsmix(pmask);
% tbl.([cfg.cond{2} 'vsmix']) = isc2vsmix(pmask);
% tbl.([cfg.cond{1} '_median']) = vert(medisc1(pmask));
% tbl.([cfg.cond{2} '_median']) = vert(medisc2(pmask));
% tbl.mix_median = vert(mediscmix(pmask));
