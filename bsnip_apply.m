function bsnip_apply(cfg,filesuffix,varargin)
bsnipdir = '/group/northoff/share/bsnip';

if ~CheckInput(varargin,'conds')
    conds = {'EO','EC'};
else
    conds = EasyParse(varargin,'conds');
end
if ~CheckInput(varargin,'phens')
    phens = {'Case','Control'};
else
    phens = EasyParse(varargin,'phens');
end
if ~CheckInput(varargin,'types')
    types = {{'Schizophrenia','Schizoaffective','Bipolar'},{'Healthy'}};
else
    types = EasyParse(varargin,'types');
end

for i = 1:length(conds)
    for ii = 1:length(phens)
        for iii=  1:length(types{ii})
            cfg.dir = fullfile(bsnipdir,conds{i},phens{ii},types{ii}{iii},'processed');
            cfg.outfile = fullfile(bsnipdir,conds{i},'measures',[conds{i} '_' types{ii}{iii} filesuffix]);
            ft_applymeasure(cfg)
        end
    end
end
