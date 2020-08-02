function grpdata = mmp_reassign(roidata,mmpgroup)

if ~exist('mmpgroup','var')
   load('mmp_atlas_grouped.mat') 
end

grpdata = roidata;
grpdata.label = mmpgroup.parcellationlabel;

for i = 1:length(grpdata.trial)
    grpdata.trial{i} = NaN(length(grpdata.label),size(roidata.trial{i},2));
    for ii = 1:length(grpdata.label)
        roicontrib = find(mmpgroup.roiassign == ii);
        nvertex = mmpgroup.nvertex(roicontrib);
        grpdata.trial{i}(ii,:) = vert(nvertex).*roidata.trial{i}(roicontrib,:)./sum(nvertex);
    end
end