function [indices,meanindices,subs] = bmaps_get_indices(subs)

nmaps = length(subs{1}.embody.bodymap);

if ~exist('mask','var') || isempty(mask)
    %mask = abs(mapmat)>prctile(reshape(abs(mapmat(abs(mapmat)>0)),[],1),1);
    mask = imread('~/Desktop/armonylab/embody-test/embody/matlab/mask.png');
    mask = [zeros(522,2) mask zeros(522,2)];
    mask = [zeros(1,175); mask; zeros(1,175)];
    inmask = find(mask > 128);
end


for i = 1:length(subs)
    subs{i}.embody.respmetrics.totaltime = NaN(1,nmaps);
    subs{i}.embody.respmetrics.drawtime = NaN(1,nmaps);
    subs{i}.embody.respmetrics.thinktime = NaN(1,nmaps);
    subs{i}.embody.respmetrics.numclicks = NaN(1,nmaps);
    if isfield(subs{i}.raw,'embody') && ~all(isnan(subs{i}.embody.origorder))
        maptrls = jspsych_result_filter(subs{i}.raw.embody,'trial_type','embody');
        %emoorder = getfield_list(maptrls,'stimulus');
        %[~,m1] = match_str(subs{1}.embody.emotion,emoorder);
        embodytottime = getfield_list(maptrls,'rt');
        try
        subs{i}.embody.respmetrics.totaltime(subs{i}.embody.origorder) = embodytottime(1:length(subs{i}.embody.origorder));
        catch
            disp('test')
            
        end
        tmp = getfield_list(maptrls,'arrTimeD');
        tmp = cellfun(@range,tmp,'UniformOutput',false);
        tmp(cellfun(@isempty,tmp)) = {NaN}; tmp = cat(2,tmp{:});
        subs{i}.embody.respmetrics.drawtime(subs{i}.embody.origorder) = tmp(1:length(subs{i}.embody.origorder));
        subs{i}.embody.respmetrics.thinktime = subs{i}.embody.respmetrics.totaltime - subs{i}.embody.respmetrics.drawtime;
        
        numclicks = getfield_list(maptrls,'arrMD');
        numclicks = cellfun(@length,numclicks,'UniformOutput',true);
        subs{i}.embody.respmetrics.numclicks(subs{i}.embody.origorder) = numclicks(1:length(subs{i}.embody.origorder));
    end
end

% make bodily map indices

for i = 1:length(subs)
    subs{i}.embody.indices.prcactarea = []; subs{i}.embody.indices.prccolored = [];
    subs{i}.embody.indices.hasdeactivation = [];
    for q = 1:length(subs{i}.embody.bodymap)
        % Percent area activated: fraction of total colored area that is
        % activation
        subs{i}.embody.indices.prcactarea(q) = (sum(sum(subs{i}.embody.bodymap{q}(inmask)>0,1),2)./(sum(sum(subs{i}.embody.bodymap{q}(inmask)<0,1),2)+nansum(nansum(subs{i}.embody.bodymap{q}(inmask)>0,1),2))).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(inmask)),1),2)));
        subs{i}.embody.indices.prccolored(q) = (sum(sum(subs{i}.embody.bodymap{q}(inmask)~=0 & ~isnan(subs{i}.embody.bodymap{q}(inmask))))./numel(inmask)).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(inmask)),1),2)));
        subs{i}.embody.indices.hasdeactivation(q) = double(any(any(subs{i}.embody.bodymap{q}(inmask)<0))).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(inmask)),1),2)));
        
        % same metrics as above but within bigatlas
        %for qq = 1:3
        %    tmpmask = find(bigatlas.atlas==qq);
        %    subs{i}.embody.parcindices.prcactarea(q,qq) = (sum(sum(subs{i}.embody.bodymap{q}(tmpmask)>0,1),2)./(sum(sum(subs{i}.embody.bodymap{q}(tmpmask)<0,1),2)+nansum(nansum(subs{i}.embody.bodymap{q}(tmpmask)>0,1),2))).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(tmpmask)),1),2)));
        %    subs{i}.embody.parcindices.prccolored(q,qq) = (sum(sum(subs{i}.embody.bodymap{q}(tmpmask)~=0 & ~isnan(subs{i}.embody.bodymap{q}(tmpmask))))./numel(tmpmask)).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(tmpmask)),1),2)));
        %    subs{i}.embody.parcindices.hasdeactivation(q,qq) = double(any(any(subs{i}.embody.bodymap{q}(tmpmask)<0))).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(tmpmask)),1),2)));
        %end
        
        ccpos = bwconncomp(subs{i}.embody.bodymap{q}>0);
        ccneg = bwconncomp(subs{i}.embody.bodymap{q}<0);
        subs{i}.embody.clusts{q} = [ccpos ccneg];
        subs{i}.embody.indices.numclusters(q) = (ccpos.NumObjects+ccneg.NumObjects).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(inmask)),1),2)));
        if ~isempty([cellfun(@length,ccpos.PixelIdxList) cellfun(@length,ccneg.PixelIdxList)])
            subs{i}.embody.indices.avgclustsize(q) = mean([cellfun(@length,ccpos.PixelIdxList) cellfun(@length,ccneg.PixelIdxList)]).*nanmask(double(all(all(~isnan(subs{i}.embody.bodymap{q}(inmask)),1),2)));
        else
            subs{i}.embody.indices.avgclustsize(q) = NaN;
        end
    end
end

for i = 1:length(subs)
    embodyexcluded(i) = isfield_nest(subs{i},'excluded.embody');
end
embodyexcluded = vert(embodyexcluded);

indices = struct;
f = fieldnames(subs{1}.embody.indices);
for i = 1:length(f)
    indices.(f{i}) = getfield_list(subs,['embody.indices.' f{i}]);
    indices.(f{i}) = cat(1,indices.(f{i}){:});
end

f = fieldnames(subs{1}.embody.respmetrics);
for i = 1:length(f)
    indices.(f{i}) = getfield_list(subs,['embody.respmetrics.' f{i}]);
    indices.(f{i}) = cat(1,indices.(f{i}){:});
    indices.(f{i}) = indices.(f{i}).*vert(nanmask(~embodyexcluded));
end

indices.nclustsclean = min(cat(3,indices.numclusters,indices.numclicks),[],3);
indices.clustsizeclean = indices.prccolored./indices.nclustsclean;

meanindices = structfun(@(d)mean(d,2),indices,'UniformOutput',false);
meanindices = struct2table(meanindices);