function out = reprandperm(n,siz)

out = zeros(siz(1),siz(2),siz(3));
for ii = 1:siz(3)
    for i = 1:siz(2)
        tmp = randperm(n,siz(1));
        out(:,i,ii) = sub2ind(siz,tmp,repmat(i,1,length(tmp)),repmat(ii,1,length(tmp)));
    end
end