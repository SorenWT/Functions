function [torf,dirn] = issig(cluster,alpha)

if nargin < 2
    alpha = cluster.cfg.alpha;
end

torf = 0;
dirn = 0;

if isfield(cluster,'posclusters')
    for c = 1:length(cluster.posclusters)
        if cluster.posclusters(c).prob <= alpha
           torf = 1;
           dirn = 1;
        end
    end
end

if isfield(cluster,'negclusters')
    for c = 1:length(cluster.negclusters)
        if cluster.negclusters(c).prob <= alpha
           torf = 1;
           dirn = dirn-1;
        end
    end
end

if torf == 0
   dirn = NaN; 
end