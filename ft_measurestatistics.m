function [stats] = ft_measurestatistics(cfg,data)
% ft_measurestatistics does stats on two (eventually or more) outputs of
% ft_applymeasure
%
% Input arguments:
%
% cfg: a config structure with the following options:
%      test: the test you want to use, passed in as a string. Inputs
%         include 'ttest' (for paired samples), 'ttest2' (for unpaired),
%         'signrank','ranksum', and 'empirical'. Use empirical for comparing
%         one outputs structure with a fourth dimension (either surrogates
%         or subsamples) with another outputs structure with only three
%         dimensions. (default = 'ranksum')
%      channel: the channels you want to test - must be an input that works
%         with ft_channelselection. (default = 'all');
%      multcompare: your method of multiple comparison correction.
%         'cluster', 'fdr', or 'none' (default = 'cluster')
%      effectsize: do effect size statistics - inputs can be any of the
%         inputs to the 'mes.m' function from the Measures of Effect Size
%         Toolbox (Hentschke and St�ttgen, 2011) (default set based on
%         test string)
%      cluster: if you chose 'cluster' as the multiple comparison
%         correction method, you can add optional inputs for the
%         permutation test
%              nrand: number of permutations (default = 10000)
%              minnbchan: minimum number of channels in a cluster (default
%              = 1)
%              statfun: the fieldtrip function to use for the permutation
%              test (default = set based on cfg.test)
%
% data: a 1 X N cell array, each cell containing one "outputs" structure
%      from ft_applymeasure. Each outputs structure must have the same
%      dimensions
%
%
% Outputs:
%
% stats: a 1 X number of measures cell array, each containing a stats
%    structure with fields depending on the input configuration

%% Set defaults

if ~cfgcheck(cfg,'test')
    cfg.test = 'ranksum';
end

if ~cfgcheck(cfg,'channel')
    cfg.channel = 'all';
end

if ~cfgcheck(cfg,'multcompare')
    cfg.multcompare = 'cluster';
end

if ~cfgcheck(cfg,'effectsize')
    switch cfg.test
        case 'ttest','ttest2'
            cfg.effectsize = 'hedgesg';
        case 'ranksum'
            cfg.effectsize = 'auroc';
        case 'signrank'
            cfg.effectsize = 'rbcorr';
    end
end

if cfgcheck(cfg,'multcompare','cluster') && ~cfgcheck(cfg.cluster,'statfun')
    switch cfg.test
        case 'ttest'
            cfg.cluster.statfun = 'ft_statfun_depsamplesT';
        case 'ttest2'
            cfg.cluster.statfun = 'ft_statfun_indepsamplesT';
        case 'ranksum'
            cfg.cluster.statfun = 'ft_statfun_ranksum';
        case 'signrank'
            cfg.cluster.statfun = 'ft_statfun_signrank';
    end
end

if cfgcheck(cfg,'multcompare','cluster') && ~cfgcheck(cfg.cluster,'nrand')
    cfg.cluster.nrand = 10000;
end

if cfgcheck(cfg,'multcompare','cluster') && ~cfgcheck(cfg.cluster,'minnbchan')
    cfg.cluster.minnbchan = 1;
end

%% Select channels
if isfield(data{1},'elec')
    chans = ft_channelselection(cfg.channel,data{1}.elec);
    chans = Subset_index(data{1}.chan,chans);
elseif isfield(data{1},'grad')
    chans = ft_channelselection(cfg.channel,data{1}.grad);
    chans = Subset_index(data{1}.chan,chans);
end

for c = 1:length(data)
    % for now, chans is always the second dimension, so no need to figure
    % out how to work with dimord
    
    %dimns = tokenize(data{1}.dimord,'_');
    %chans = find(strcmpi(dimns,'chan'));
    data{c} = data{c}.data(:,chans,:,:);
end

%% Calculate stats
for i = 1:length(data{1}.meas)
    %dimns = tokenize(data{1}.dimord,'_');
    for c = 1:length(data{1}.chan)
        switch cfg.test
            case 'ttest'
                [~,stats{i}.p(c)] =  ttest(data{1}.data(:,c,i)-data{2}.data(:,c,i));
            case 'ttest2'
                [~,stats{i}.p(c)] = ttest2(data{1}.data(:,c,i),data{2}.data(:,c,i));
            case 'ranksum'
                stats{i}.p(c) = ranksum(data{1}.data(:,c,i),data{2}.data(:,c,i));
            case 'signrank'
                stats{i}.p(c) = signrank(data{1}.data(:,c,i),data{2}.data(:,c,i));
            case 'empirical'
                % finish later
        end
        stats{i}.effsize{c} = mes(data{1}(:,c,i),data{2}(:,c,i),cfg.effectsize);
    end
    
    switch cfg.multcompare
        case 'cluster'
            if isfield(data{1},'elec')
                datasetinfo.elec = data{1}.elec;
            elseif isfield(data{1},'grad')
                datasetinfo.grad = data{1}.grad;
            end
            datasetinfo.label = data{1}.chan;
            for c = 1:length(data)
                input{c} = data{c}.data(:,:,i)';
            end
            stats{i}.cluster = EasyClusterCorrect(input,datasetinfo,cfg.cluster.statfun,cfg.cluster);
        case 'fdr'
            stats{i}.fdr = mafdr(stats{i}.p(c),'BHFDR',true);
    end
end


