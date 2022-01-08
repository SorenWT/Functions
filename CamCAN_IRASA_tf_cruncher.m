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
for i = indx
    subid = char(extractBefore(files(i).name,'_'));
    try
        grad = parload(fullfile('/group/northoff/share/camcan/Task/Continuous/',[subid '_continuous_1Hz.mat']),'cont_data');
        grad = grad.grad;
        
        
        if ~exist(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmpreduced.mat'],'file')
            load(files(i).name)
            data.grad = grad;
            
            
            headmodel = parload(fullfile('/home/soren/Documents/camcan/sourcemodels/',subid,'sourcemodel',...
                [subid '_headmodel.mat']),'vol');
            sourcemodel = parload(fullfile('/home/soren/Documents/camcan/sourcemodels/',subid,'sourcemodel',...
                [subid '_sourcemodel_8k_nmg.mat']),'sourcemodel');
            
            roidata = SourceEst(data,headmodel,sourcemodel,atlas,struct('parflag',1));
            
            save(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmp.mat'],'roidata')
            
            roidata = mmp_reassign(roidata,mmpgroup);
            
            save(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmpreduced.mat'],'roidata')
            
            data = roidata;
        else
            load(['/home/soren/Documents/camcan/source/task/epoched/' files(i).name(1:end-4) '_roidata_mmpreduced.mat'])
            data = roidata;
        end
        
        
        if ~exist([files(i).name '_IRASAtf.mat'],'file')
            %load(files(i).name)
            
            latencies = data.sampleinfo(:,1)'+1000;
            
            itis = diff([latencies(1) latencies]); % don't keep the first trial since you don't know how long the first ITI really was
            goodevents = find(itis >=2000); % keep trials where the preceeding iti was 4 seconds or longer
            %ntrl(i) = length(goodevents);
            
            data = ft_selectdata(struct('trial',goodevents),data);
            
            tic;
            [osci,specs] = IRASA_tf(cfg,data);
            toc;
            cfg.oscifrac = 'frac';
            frac = IRASA_tf(cfg,data,specs);
            cfg.oscifrac = 'mixd';
            mixd = IRASA_tf(cfg,data,specs);
            save([files(i).name '_IRASAtf.mat'],'osci','frac','mixd','specs')
        end
    catch
        badsubs(i) = 1;
    end
end
