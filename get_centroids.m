function centroids = get_centroids(data,cl)

uniquecl = unique(cl);

for i = 1:length(uniquecl)
    X = data(cl==i,:);
    X = bsxfun(@minus, X, mean(X,2));
    Xnorm = sqrt(sum(X.^2, 2));
    if any(min(Xnorm) <= eps(max(Xnorm)))
        error(message('stats:kmeans:ConstantDataForCorr'));
    end
    X =  bsxfun(@rdivide,X,Xnorm);
end