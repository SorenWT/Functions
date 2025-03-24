function [thisplot] = plotInteraction_swt(mdl,variables,type,p,pindx)

if nargin < 3
    type = 'predictions';
end

if nargin < 4
   p = panel('no-manage-font');
   p.pack('h',1);
end

if isa(mdl,'LinearMixedModel')
    oldmdl = mdl;
    
    newformula = char(oldmdl.Formula); 
    newformula = erase(newformula,' ');
    newformula = regexprep(newformula, '\(.*?\)', '');
    
    % remove the random effects
    
    mdl = fitlm(oldmdl.Variables,newformula);
end

if length(variables)==2
    thisplot = plotInteraction(mdl,variables{1},variables{2},type); 
else
    p = panel('no-manage-font');
    p(pindx{:}).pack('h',{1/2 1/2})
    
    
end
