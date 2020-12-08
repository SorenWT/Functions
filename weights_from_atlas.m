function W = weights_from_atlas(atlas,info)
% computes spatial weight matrix for a volumetric atlas for Moran spectral
% randomization
% atlas is a 3D image
% info is a nifti header containing the spatial transform matrix

% first, make the neighbour matrix
atlas(find(atlas == 0)) = NaN;
%label = cellstr(num2str([1:length(unique(atlas(~isnan(atlas))))]'));
label = 1:length(unique(atlas(~isnan(atlas))));

neighbs = struct;
for c = 1:length(label)
    neighbs(c).label = label(c);
    neighbs(c).neighblabel = [];
end

for c = 1:3
    edges{c} = diff(atlas,1,c);
    edges{c} = edges{c} > 0;
    boundaryindx = find(edges{c});
    for cc = 1:length(boundaryindx)
        [i1,i2,i3] = ind2sub(size(edges{c}),boundaryindx(cc));
        indx = {i1 i2 i3};
        indx{c} = indx{c}+1;
        roi1val = atlas(i1,i2,i3);
        roi2val = atlas(indx{:});
        neighbs(roi1val).neighblabel = [neighbs(roi1val).neighblabel neighbs(roi2val).label];
        neighbs(roi2val).neighblabel = [neighbs(roi2val).neighblabel neighbs(roi1val).label];
    end
end

for c = 1:length(neighbs)
    %neighbs(c).neighblabel(1) = [];
    neighbs(c).neighblabel = unique(neighbs(c).neighblabel);
end

atlas(isnan(atlas))=0;

W = eye(length(unique(atlas(atlas>0))));

for i = 1:length(neighbs)
    W(i,neighbs(i).neighblabel) = 1;
    W(neighbs(i).neighblabel,i) = 1;
end

inside = atlas>0; uniquevals = unique(atlas(atlas>0));
% now figure out the spatial distances
dim = size(atlas); 
transform = info.Transform.T;
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos_from  = ft_warp_apply(sparse(transform), [X(:) Y(:) Z(:)]);
pos_from  = pos_from(inside,:);

% get mean positions for each voxel
atlas1d = reshape(atlas,[],1);
atlas1d = atlas1d(inside);
for i = 1:length(uniquevals)
    meanpos(i,:) = mean(pos_from(atlas1d==uniquevals(i),:),1);
end
distmat = zeros(size(W));
for i = 1:length(uniquevals)
    for ii = 1:(i-1)
        distmat(i,ii) = sum((meanpos(i,:)-meanpos(ii,:)).^2);
    end
end
distmat = distmat+distmat'; % make symmetric

% now weight the adjacency matrix by distance
W = W./distmat; 
W(isinf(W)) = 0;

