function X = nancenter(X,dim)


if istable(X)
    X{:,:} = X{:,:}-nanmean(X{:,:},dim);
else
X = X-nanmean(X,dim);
end