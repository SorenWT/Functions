%function CamCAN_IRASA_tf_slurm(indx)
% script for getting IRASA time-frequency stuff

addpath(genpath('/home/soren/Documents/MATLAB'))
addpath('/group/northoff/share/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

cd /home/soren/Documents/camcan/Preprocessed/Task/Epoched

files = dir('*1Hz.mat');

%indices = ((indx-1)*32+1):(indx*32);
rng(0)
indx = randperm(length(files),200);

%indices(indices > length(files)) = [];

cfg.oscifrac = 'osci'; cfg.winsize = 1.5; cfg.modifywindow = 'no';
cfg.toi = [-1.25:0.025:-0.75 0:0.025:0.75];
cfg.parflag = 'yes';
%cfg.toi = [-0.75:0.025:0.75];

load('atlas_MMP1.0_4k.mat')
load('mmp_atlas_grouped.mat')

%parpool(48)
%count = 1;
badcount = 1;

for i = indx
    %i = indx(count);
    
    subid = char(extractBefore(files(i).name,'_'));
    try
        grad = parload(fullfile('/group/northoff/share/camcan/Task/Continuous/',[subid '_continuous_1Hz.mat']),'cont_data');
        grad = grad.grad;
    catch
        badsubs(badcount).message = 'No continuous task file found';
        badsubs(badcount).indx = i;
        badcount = badcount+1;
        continue
    end
    
    
    if ~exist(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmpreduced.mat'],'file')
        load(files(i).name)
        data.grad = grad;
        
        
        try
            headmodel = parload(fullfile('/home/soren/Documents/camcan/sourcemodels/',subid,'sourcemodel',...
                [subid '_headmodel.mat']),'vol');
        catch
            badsubs(badcount).message = 'No headmodel file found';
            badsubs(badcount).indx = i;
            badcount = badcount+1;
            continue
        end
        try
            sourcemodel = parload(fullfile('/home/soren/Documents/camcan/sourcemodels/',subid,'sourcemodel',...
                [subid '_sourcemodel_8k_nmg.mat']),'sourcemodel');
        catch
            badsubs(badcount).message = 'No sourcemodel file found';
            badsubs(badcount).indx = i;
            badcount = badcount+1;
            continue
        end
        
        roidata = SourceEst(data,headmodel,sourcemodel,atlas,struct('parflag',1));
        
        roidata.sampleinfo = data.sampleinfo; roidata.trialinfo = data.trialinfo;
        
        save(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmp.mat'],'roidata')
        
        roidata = mmp_reassign(roidata,mmpgroup);
        roidata.sampleinfo = data.sampleinfo; roidata.trialinfo = data.trialinfo;
        
        
        save(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmpreduced.mat'],'roidata')
        
        data = roidata;
    else
        load(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmpreduced.mat'])
        data = roidata;
    end
    
    
    if ~exist(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name '_IRASAtf_source.mat'],'file')
        %load(files(i).name)
        
        latencies = data.sampleinfo(:,1)'+1000;
        
        itis = diff([latencies(1) latencies]); % don't keep the first trial since you don't know how long the first ITI really was
        goodevents = find(itis >=2000); % keep trials where the preceeding iti was 4 seconds or longer
        %ntrl(i) = length(goodevents);
        
        data = ft_selectdata(struct('trials',goodevents),data);
        
        disp(['Data has ' num2str(length(data.trial)) ' trials'])
        
        tic;
        [osci,specs] = IRASA_tf(cfg,data);
        eltime(i) = toc;
        cfg.oscifrac = 'frac';
        frac = IRASA_tf(cfg,data,specs);
        cfg.oscifrac = 'mixd';
        mixd = IRASA_tf(cfg,data,specs);
        save(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name '_IRASAtf_source.mat'],'osci','frac','mixd','specs')
    end
    
    %catch
    %   badsubs(i) = 1;
    %end
end
