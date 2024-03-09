function [cv_p,acc,mdl_train] = crossval_boot(indvar,depvar,cvfunc,accfunc,opts)

if nargin < 5
   opts = struct;
end

opts = setdefault(opts,'critacc',0.5);
opts = setdefault(opts,'KFold',2);
opts = setdefault(opts,'perms',1000);
opts = setdefault(opts,'predfun',@(d1,d2)predict(d1,d2));

clear acc
for i = 1:round(opts.perms/opts.KFold)
    cvpart = cvpartition(depvar,'KFold',opts.KFold);
    for ii = 1:opts.KFold
        trainind = find(training(cvpart,ii));
        testind = find(test(cvpart,ii));
        
        mdl_train{i,ii} = cvfunc(indvar(trainind,:),depvar(trainind));
        
        pred_test = opts.predfun(mdl_train{i,ii},indvar(testind,:));
        
        acc(i,ii) = accfunc(pred_test,depvar(testind));
    end
    
end

acc = reshape(acc,[],1);

cv_p = 1-sum(acc>opts.critacc)./length(acc); 

