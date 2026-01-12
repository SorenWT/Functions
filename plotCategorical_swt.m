function [thisplot] = plotCategorical_swt(mdl,variable,opts)

if nargin < 3
    opts = struct;
end

opts = setdefault(opts,'type','predictions');
opts = setdefault(opts,'nstrat', 3);
opts = setdefault(opts,'stratprc',linspace(0,100,opts.nstrat));
opts = setdefault(opts,'labels',[{mdl.ResponseName} variable]);


if nargin < 4
   p = panel('no-manage-font');
   p.pack('h',1);
end

if isa(mdl,'LinearMixedModel')
    oldmdl = mdl;
    
    newformula = char(oldmdl.Formula); 
    newformula = erase(newformula,' ');
    % remove the random effects
    newformula = regexprep(newformula, '\(.*?\)', '');
    newformula(end) = []; % only will work with one random effect, fix later
    
    mdl = fitlm(oldmdl.Variables,newformula);
end

for i = 1:length(variable)
    varinds(i) = find(strcmpi(mdl.Variables.Properties.VariableNames,variable{i}));    
end






end


function plotit(mdl,varinds,xdat,opts)


for i = 1:length(xdat)
    y = getAdjustedResponse(mdl,varinds,xdat{i});
    plot(xdat{i}(:,end),y,'LineWidth',2)
    hold on
    legendinfo{i} = [opts.labels{end-1} ' = ' num2str(xdat{i}(1,end-1))];
end
xlabel(opts.labels{end})
ylabel(opts.labels{1})

legend(legendinfo)

if size(xdat{1},2)>2
    ttl = '';
    for i = 1:(length(varinds)-2)
        ttl = [ttl opts.labels{end-1-i}  ' = ' num2str(xdat{1}(1,end-1-i)) ' & '];
    end
    
    ttl(end-1:end) = [];
    title(ttl)
end
FixAxes(gca,18)

end
