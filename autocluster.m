function [bestclust,bestctr,allclusts,ctrs] = autocluster(data,nrange,opts)

if nargin < 4
   opts = struct; 
end

opts = setdefault(opts,'nperms',1000);
opts = setdefault(opts,'trainprc',0.5);
opts = setdefault(opts,'distance','euclidean');


for i = nrange
    [allclusts{i},ctrs{i}] = kmeans(data,i,'Distance',opts.distance);
end


for i = 1:opts.nperms
    cvp = cvpartition(size(data,2),'Holdout',opts.trainprc);
    train = training(cvp); test = test(cvp);
    for ii = nrange
        [allclusts_train{ii}] = kmeans(data(:,train),ii,'Distance',opts.distance);
        [allclusts_test{ii}] = kmeans(data(:,test),ii,'Distance',opts.distance);
        trainadj = allclusts_train{ii}==allclusts_train{ii}';
        testadj = allclusts_test{ii}==allclusts_test{ii}';
        
        stab(i,ii) = dice(belowDiag(trainadj),belowDiag(testadj));
    end
    
end

meanstab = mean(stab,1);
[~,bestindx] = max(meanstab);

bestclust = allclusts{bestindx}; bestctr = ctrs{bestindx};
    
    
    