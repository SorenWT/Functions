function [symdata]=symmetrize(data,dir,pixonly)

if nargin < 3
    pixonly = 0;
end

if pixonly
    comp = bwconncomp(data); 
    [~,maxpix] = max(cellfun(@length,comp.PixelIdxList));
    comp.PixelIdxList = comp.PixelIdxList(maxpix); comp.NumObjects = 1;
    theseprops = regionprops(comp,{'BoundingBox'});
    theseprops.BoundingBox = ceil(theseprops.BoundingBox);
    theseprops.bbindx = {theseprops.BoundingBox(1):theseprops.BoundingBox(1)+theseprops.BoundingBox(3),...
        theseprops.BoundingBox(2):theseprops.BoundingBox(2)+theseprops.BoundingBox(4)}
    
    tmp = symmetrize(data(theseprops.bbindx{2},theseprops.bbindx{1}),dir);
    data(theseprops.bbindx{2},theseprops.bbindx{1}) = tmp;
    symdata = data;
else
    symdata = (data+flip(data,dir))/2;
end
