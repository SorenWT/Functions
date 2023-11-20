function [compsout] = applypca(pcamdl,tblin,dorot)

if nargin < 3
    if isfield(pcamdl,'rotated')
        dorot = 1;
    else
        dorot = 0;
    end
end

if istable(tblin)
    [m1,m2] = match_str(pcamdl.varnames,tblin.Properties.VariableNames);
    
    if dorot
        weights = pcamdl.rotated.weights(m1,:);
    else
        weights = pcamdl.weights(m1,:);
    end
    tblin = tblin(:,m2);
    
    compsout = nancenter(tblin{:,:},1)*weights;
elseif ismatrix(tblin)
    
end