function [thisplot] = plotInteraction_swt(mdl,variables,opts)

if nargin < 3
    opts = struct;
end

opts = setdefault(opts,'type','predictions');
opts = setdefault(opts,'nstrat', repmat(3,1,length(variables)));
for i = 1:length(variables)
    tmp{i} = linspace(0,100,opts.nstrat(i));
end
opts = setdefault(opts,'stratprc',tmp);
opts = setdefault(opts,'labels',[{mdl.ResponseName} variables]);
opts = setdefault(opts,'DummyVarCoding','full');
opts = setdefault(opts,'plotdata','off');


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
    
    mdl = fitlm(oldmdl.Variables,newformula,'DummyVarCoding',opts.DummyVarCoding);
    %mdl.Coefficients = oldmdl.Coefficients;
    mdl = oldmdl; 
end


% make a writeable copy
varstbl = mdl.Variables;

for i = 1:length(variables)
    varinds(i) = find(strcmpi(varstbl.Properties.VariableNames,variables{i})); 
    if iscategorical(varstbl.(variables{i})) || iscell(varstbl.(variables{i}))
        varstbl.(variables{i}) = factor2num(varstbl.(variables{i}))';
        opts.nstrat(i) = length(unique(varstbl.(variables{i})(~isnan(varstbl.(variables{i})))));
    end
end

% actually plot it all
p = panel('no-manage-font');

tmp = zeros(1,2);
tmp(1) = min(varstbl.(variables{end})); tmp(2) = max(varstbl.(variables{end}));

if length(variables) == 4
    p.pack(opts.nstrat,opts.nstrat);
    
    for q = 1:opts.nstrat(1)
        for qq = 1:opts.nstrat(2)
            for qqq = 1:opts.nstrat(3)
                xdat{qqq}(:,length(variables)) = make_xdat(varstbl.(variables{end}));
%                 xdat{qqq} = zeros(100,length(variables));
%                 xdat{qqq}(:,end) = linspace(tmp(1),tmp(2),100);
                xdat{qqq}(:,3) = thisprctile(varstbl.(variables{3}),opts.stratprc{3},qqq);
                xdat{qqq}(:,2) = thisprctile(varstbl.(variables{2}),opts.stratprc{2},qq);
                xdat{qqq}(:,1) = thisprctile(varstbl.(variables{1}),opts.stratprc{1},q);
            end
            
            p(q,qq).select()
            plotit(mdl,varinds,xdat,varstbl,opts)
        end
    end
elseif length(variables)==3
    p.pack(1,opts.nstrat(1));
    for q = 1:opts.nstrat(1)
        for qq = 1:opts.nstrat(2)
%             xdat{qq} = zeros(100,length(variables));
%             xdat{qq}(:,end) = linspace(tmp(1),tmp(2),100);
            xdat{qq}(:,length(variables)) = make_xdat(varstbl.(variables{end}));
            xdat{qq}(:,2) = thisprctile(varstbl.(variables{2}),opts.stratprc{2},qq);
            xdat{qq}(:,1) = thisprctile(varstbl.(variables{1}),opts.stratprc{1},q);
        end
        
        p(1,q).select()
        plotit(mdl,varinds,xdat,varstbl,opts)
    end
else
    p.pack(1,1);
    for q = 1:opts.nstrat(1)
        %xdat{q} = zeros(100,length(variables));
        %xdat{q}(:,end) = linspace(tmp(1),tmp(2),100);
        xdat{q}(:,length(variables)) = make_xdat(varstbl.(variables{end}));
        xdat{q}(:,1) = thisprctile(varstbl.(variables{1}),opts.stratprc{1},q);
    end
        
    p(1,1).select()
    plotit(mdl,varinds,xdat,varstbl,opts)
end


p.margin = [25 25 5 10];



end


function plotit(mdl,varinds,xdat,varstbl,opts)


for i = 1:length(xdat)
    %y = feval(mdl,varinds,xdat{i});
    for ii = 1:size(xdat{i},2)
        vdat{ii} = xdat{i}(:,ii);
        if iscategorical(mdl.Variables{:,varinds(ii)}) || iscell(mdl.Variables{:,varinds(ii)})
            vdat{ii} = num2factor(vdat{ii},unique(mdl.Variables{:,varinds(ii)}),1:length(unique(mdl.Variables{:,varinds(ii)})));
        end
    end
    
    variables = mdl.VariableNames(varinds);
    
    %[~,sortord] = match_str(mdl.PredictorNames,variables);
    
    vtbl = table;
    for q = 1:length(vdat)
        vtbl.(variables{q}) = vdat{q};
    end
    
    if isa(mdl,'LinearMixedModel')
        for q = 1:length(mdl.Formula.GroupingVariableNames)
            vtbl.(mdl.Formula.GroupingVariableNames{q}{1}) = repmat({'new'},height(vtbl),1);
        end
    end
    
    for q = 1:length(mdl.PredictorNames)
        if ~any(strcmpi(vtbl.Properties.VariableNames,mdl.PredictorNames{q}))
           vtbl.(mdl.PredictorNames{q}) = repmat(nanmean(mdl.Variables.(mdl.PredictorNames{q})),height(vtbl),1); 
        end
    end
    
    %y = feval(mdl,vdat{sortord});
    y = predict(mdl,vtbl);
    
    l = lines;
    plot(xdat{i}(:,end),y,'LineWidth',2,'Color',l(i,:))
    hold on
    
    if ischar(vdat{end-1}) || iscell(vdat{end-1})
        val = vdat{end-1}{1};
    else
            val = num2str(xdat{i}(1,end-1));
    end
        
    legendinfo{i} = [opts.labels{end-1} ' = ' val];
    
    if strcmpi(opts.plotdata,'on')
        % only works right now for categorical moderators - need to figure
        % out for continuous
        for q = 1:length(variables)-1
           indx(:,q) = varstbl.(variables{q}) == vtbl.(variables{q})(1);
        end
        indx = all(indx,2);
        sz = 250./round(log(sum(indx)));
        
        if isa(mdl,'LinearMixedModel')
           scatter(grpstats(varstbl.(variables{end})(indx),varstbl.(vtbl.Properties.VariableNames{end})(indx)),grpstats(varstbl.(mdl.ResponseName)(indx),varstbl.(vtbl.Properties.VariableNames{end})(indx)),...
                sz,palecol(l(i,:),0.3),'filled') 
        else
            scatter(varstbl.(variables{end})(indx),varstbl.(mdl.ResponseName)(indx),sz,palecol(l(i,:),0.3),'filled') 
        end
    end

end
xlabel(opts.labels{end})
ylabel(opts.labels{1})

manlegend(legendinfo,l(1:length(xdat),:))

if size(xdat{1},2)>2
    ttl = '';
    for i = 1:(length(varinds)-2)
        ttl = [ttl opts.labels{end-1-i}  ' = ' num2str(xdat{1}(1,end-1-i)) ' & '];
    end
    
    ttl(end-1:end) = [];
    title(ttl)
end
FixAxes(gca,18)
Normalize_Ylim(gcf,0)

end

function out = thisprctile(var,stratprc,indx)

if length(unique(var(~isnan(var)))) < 10
    unvals = unique(var(~isnan(var)));
    out = unvals(indx);
else
    out = nanprctile(var,stratprc(indx));
end
end

function xdat = make_xdat(var)

% pick an arbitrary threshold for the number of unique values for it to
% consider a factor
if length(unique(var(~isnan(var)))) < 10
xdat = unique(var(~isnan(var)));
else
    tmp(1) = min(var); tmp(2) = max(var);

    xdat = linspace(tmp(1),tmp(2),100);
end

end
