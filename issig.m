function [torf,dirn,p] = issig(cluster,alpha)

if nargin < 2
    alpha = cluster.cfg.alpha;
end

torf = 0;
dirn = 0;

p = [];
count = 1;
if isfield(cluster,'posclusters')
    for c = 1:length(cluster.posclusters)
        if cluster.posclusters(c).prob <= alpha
           torf = 1;
           dirn = 1;
           p(count) = cluster.posclusters(c).prob*0.05/alpha;
           count = count+1;
        end
    end
end

if isfield(cluster,'negclusters')
    for c = 1:length(cluster.negclusters)
        if cluster.negclusters(c).prob <= alpha
           torf = 1;
           dirn = dirn-1;
           p(count) = cluster.negclusters(c).prob*0.05/alpha;
           count = count+1;
        end
    end
end

if torf == 0
   dirn = NaN; 
end