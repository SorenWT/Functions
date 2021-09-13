function [shuford,shufseq] = balance_conditions(overseq,underseqs,wincondtol,crosscondtol,balancetransitions)

if nargin < 4 || isempty(crosscondtol)
    crosscondtol = 1;
end

if nargin < 3 || isempty(wincondtol)
    wincondtol = 0.05;
end

if nargin < 5
    balancetransitions = NaN;
end

underseqs = horz(underseqs);
overseq = horz(overseq);

if size(underseqs,2) < size(overseq,2)
    underseqs = repmat(underseqs,1,ceil(size(overseq,2)./size(underseqs,2)));
end

if size(overseq,1)>1
    % if you want to balance something with respect to multiple conditions,
    % convert it to a single condition with more values
    if size(overseq,2) > 1
        overseq = overseq.*vert(10.^[1:size(overseq,1)]);
        overseq = sum(overseq,1);
    end
end

if size(underseqs,1)>1 && ~all(isnan(balancetransitions))
    if size(underseqs,2) > 1
        transbalseq = underseqs.*vert(10.^[1:size(underseqs,1)]);
        transbalseq = sum(transbalseq(balancetransitions,:),1);
    end
end

while 1
    if all(isnan(balancetransitions))
        shuford = randperm(size(underseqs,2),length(overseq));
        shufseq = underseqs(:,shuford);
    else
        uniquevals = unique(transbalseq);
        nreps = floor(size(overseq,2)./numel(ones(length(uniquevals)))); 
        if nreps > 0
        shufseq = carryoverCounterbalance_SWT(length(uniquevals),1,nreps,0);
        else
            shufseq = 0;
        end
        shufseq = shufseq(1:end-1); 
        %shufseq = shufseq(1:length(transbalseq));
        remaining = length(overseq)-length(shufseq);
        for q = 1:(remaining/length(uniquevals))
           shufseq = [shufseq randperm(length(uniquevals),length(uniquevals))]; 
        end
        
        shuford = zeros(1,length(shufseq));
        
        while 1
            for i = 1:length(uniquevals)
                tmp = find(transbalseq==uniquevals(i));
                shuford(shufseq==i) = tmp(randperm(length(tmp),sum(shufseq==i)));
            end
            newshufseq = underseqs(:,shuford);
            if ~any(~balancetransitions) || check_balance(transbalseq(shuford),newshufseq(~balancetransitions,:),wincondtol,crosscondtol)
                shufseq = newshufseq;
                break
            end
        end
    end
    
    if check_balance(overseq,shufseq,wincondtol,crosscondtol)
        break 
    end

end

end

function isbalanced = check_balance(overseq,shufseq,wincondtol,crosscondtol)

for i = 1:size(shufseq,1)
    counts{i} = count_conds(overseq,shufseq(i,:));
end

for i = 1:length(counts)
    for ii = 1:size(counts{i},1)
        %crit(i,ii) = (sum(abs(diff(counts{i})),2)./sum(counts{i},2))<tolerance;
        nexpected = ceil(sum(counts{i}(ii,:))./length(counts{i}(ii,:)));
        critval(i,ii) = (abs(max(counts{i}(ii,:))-nexpected))./nexpected;
        crit(i,ii) = critval(i,ii)<=wincondtol;
    end
end

if sum(sum(crit))./numel(crit) >= crosscondtol
    isbalanced = 1;
else
    isbalanced = 0;
end


end


function counts = count_conds(overseq,underseq)
overseqels = unique(overseq);
underseqels = unique(underseq);

for q = 1:length(overseqels)
    for qq = 1:length(underseqels)
        counts(q,qq) = sum(overseq==overseqels(q) & underseq==underseqels(qq));
    end
end
end