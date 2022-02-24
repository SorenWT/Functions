function clust_sigmasks(ax,clust,reverse)
% reverse flips the order so that blue is positive and red is negative

monored = customcolormap_preset('red-white-blue');
monoblue = monored(1:32,:);
monored = flipud(monored(33:end,:));

if nargin < 3
   reverse = 0; 
end

if reverse
   tmp = monoblue;
   monoblue = monored;
   monored = tmp;
   offsets = [0.08 0.04];
else
    offsets = [0.04 0.08];
end

if ~isempty(clust.posclusters)
    signums = extractfield(clust.posclusters,'prob');
    signums = find(signums < 0.05);
    if ~isempty(signums)
        Plot_sigmask(ax,ismember(clust.posclusterslabelmat,signums),'cmapline','LineWidth',5,'cmap',monored,'yoffset',offsets(1));
    end
end

if ~isempty(clust.negclusters)
    signums = extractfield(clust.negclusters,'prob');
    signums = find(signums < 0.05);
    if ~isempty(signums)
        Plot_sigmask(ax,ismember(clust.negclusterslabelmat,signums),'cmapline','LineWidth',5,'cmap',monoblue,'yoffset',offsets(2));
    end
end