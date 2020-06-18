function p = clusterp(clust,minp)

if isfield(clust,'negclusters') && ~isempty(clust.negclusters)
    p.neg = clust.negclusters(1).prob;
else
    p.neg = NaN;
end

if isfield(clust,'posclusters') && ~isempty(clust.posclusters)
   p.pos = clust.posclusters(1).prob;
else
   p.pos = NaN;
end

if nargin > 1 && strcmpi(minp,'min')
	p = min(p.pos,p.neg);
end
