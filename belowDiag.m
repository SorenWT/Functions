function outmat = belowDiag(inmat)

if size(inmat,3) == 1
    mask = tril(ones(size(inmat)))-eye(size(inmat));
    %if ~iscell(inmat)
    %    maskedmat = inmat.*mask;
    %end
    %outmat = reshape(inmat(sub1,sub2,i),[],1);
    outmat = reshape(inmat(find(mask)),[],1);
else
    for i = 1:size(inmat,3)
        outmat(:,i) = belowDiag(inmat(:,:,i));
    end
end