function ctrs = getcentroids(cl,dat,dist)

k = length(unique(cl));

switch dist
    case 'cosine'
        tmpdat = dat./vecnorm(dat,2,2);
        for i = 1:k
            ctrs(i,:) = sum(tmpdat(cl==i,:),1);
            ctrs(i,:) = ctrs(i,:)/norm(ctrs(i,:));
            ctrs(i,:) = ctrs(i,:)*mean(vecnorm(dat(cl==i,:),2,2));
        end
    case 'sqeuclidean'
        for i = 1:k
            ctrs(i,:) = mean(dat(cl==i,:),1);
        end
        
        
end

