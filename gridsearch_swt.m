function [xf,fvals,gridvals,gridindx,indx] = gridsearch_swt(fun,nvars,lb,ub,opts)

lb = horz(lb); ub = horz(ub);

if nargin < 5
    opts = struct;
end

opts = setdefault(opts,'UseParallel',true);
opts = setdefault(opts,'nsteps',10);

for i = 1:nvars
    gridvals(i,:) = linspace(lb(i),ub(i),opts.nsteps);
end

derivopts = opts; derivopts.original = 0;

if nvars > 1
    for i = 1:opts.nsteps
        [~,fvals{i}] = gridsearch_swt(@(d)fun(cat(2,gridvals(1,i),d)),nvars-1,lb(2:end),ub(2:end),derivopts);
    end
else
    if opts.UseParallel
        parfor i = 1:opts.nsteps
            fvals{i} = fun(gridvals(1,i));
        end
    else
        for i = 1:opts.nsteps
            fvals{i} = fun(gridvals(1,i));
        end
    end
end

xf = 0;

if opts.original
    fvalsold = fvals;
    for i = 1:nvars
        fvals = cat(2,fvals{:});
    end
    [~,gridindx] = min(fvals);
    for i = 1:nvars
       indx(i,:) = repmat(horz(Make_designVect(repmat(opts.nsteps^(nvars-i),1,opts.nsteps))),1,opts.nsteps^(i-1));
    end
    gridindx = indx(:,gridindx);
    for i = 1:nvars
       xf(i) = gridvals(i,gridindx(i)); 
    end
end


