function seq = counterbalance_multicond(nconds,reps)
% only works with 1st order counterbalancing, omits self adjacencies
% nconds is a cell array with the different conditions
% works only with 2 conditions right now

totalconds = prod([nconds{:}]);

adjmatrix = ones(totalconds,totalconds);

condmatrix = reshape(1:totalconds,nconds{:});

for i = 1:nconds{1}
    for q = 1:nconds{2}
        for qq = 1:nconds{2}
            adjmatrix(condmatrix(i,q),condmatrix(i,qq)) = 0;
        end
    end
end

for i = 1:nconds{2}
    for q = 1:nconds{1}
        for qq = 1:nconds{1}
            adjmatrix(condmatrix(q,i),condmatrix(qq,i)) = 0;
        end
    end
end

seq = carryoverCounterbalance_SWT(totalconds,1,reps,1,adjmatrix);   