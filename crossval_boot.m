function [cv_p,acc] = crossval_boot(indvar,depvar,cvfunc,accfunc,opts)

if nargin < 5
   opts = struct;
end

opts = setdefault(opts,'critval',0.1);
opts = setdefault(opts,'KFold',2);
opts = setdefault(opts,'perms',1000);

clear acc
for i = 1:round(opts.perms/opts.KFold)
    cvpart = cvpartition(length(depvar)/10,'KFold',opts.KFold);
    for ii = 1:opts.KFold
        trainind = ismember(grp,find(training(cvpart,ii)));
        testind = ismember(grp,find(test(cvpart,ii)));
        
        mdl_train = cvfunc(indvar(trainind),depvar(trainind));
        
        pred_test = predict(mdl_train,indvar(testind));
        
        acc(i,ii) = accfunc(pred_test,depvar(testind));
    end
    
end

acc = reshape(acc,[],1);

cv_p = 1-sum(acc>opts.critval)./length(acc); 

